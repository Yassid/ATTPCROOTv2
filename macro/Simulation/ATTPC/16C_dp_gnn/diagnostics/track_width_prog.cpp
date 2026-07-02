// track_width_prog.cpp  (standalone, links libAtData -> no rootcling/ACLiC)
// Transverse track thickness vs drift z on a raw AtEventH point cloud.
// Build: see build_width.sh.  Run: ./track_width <file.root> <tag> [vDrift] [maxEvt] [fwdDeg]
#include "AtEvent.h"
#include "AtHit.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>

#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

struct WPt {
   double x, y, z, q;
};

int main(int argc, char **argv)
{
   if (argc < 3) {
      printf("usage: %s <file.root> <tag> [vDrift=1.15] [maxEvt=-1] [fwdDeg=40] [R=8] [Nmin=6] [linCut=3]\n", argv[0]);
      return 1;
   }
   std::string fname = argv[1];
   std::string tag = argv[2];
   double vDrift = argc > 3 ? atof(argv[3]) : 1.15;
   long maxEvt = argc > 4 ? atol(argv[4]) : -1;
   double fwdDeg = argc > 5 ? atof(argv[5]) : 40.0;
   double R = argc > 6 ? atof(argv[6]) : 8.0;
   int Nmin = argc > 7 ? atoi(argv[7]) : 6;
   double linCut = argc > 8 ? atof(argv[8]) : 3.0;

   TFile *f = TFile::Open(fname.c_str());
   if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", fname.c_str()); return 1; }
   TTree *tree = (TTree *)f->Get("cbmsim");
   if (!tree) { printf("ERROR: no cbmsim tree\n"); return 1; }
   TClonesArray *evtArr = nullptr;
   tree->SetBranchAddress("AtEventH", &evtArr);

   double zmin = 1e9, zmax = -1e9;
   // pad-plane transverse spread (-> CoefT) and drift-aligned transverse spread (-> CoefL), vs drift z
   TProfile pXY(Form("p_xy_%s", tag.c_str()), Form("%s;z [mm];#sigma_{xy}^{2} [mm^{2}]", tag.c_str()), 50, 0, 1000);
   TProfile pZd(Form("p_zd_%s", tag.c_str()), Form("%s;z [mm];#sigma_{z}^{2} [mm^{2}]", tag.c_str()), 50, 0, 1000);
   TH1D hQ(Form("h_q_%s", tag.c_str()), Form("%s per-hit charge;charge;hits", tag.c_str()), 200, 0, 8000);
   TH1D hWxy(Form("h_wxy_%s", tag.c_str()), Form("%s #sigma_{xy};#sigma_{xy} [mm];hits", tag.c_str()), 100, 0, 8);
   TH1D hTh(Form("h_th_%s", tag.c_str()), Form("%s local polar angle;#theta [deg];hits", tag.c_str()), 90, 0, 90);

   long n = tree->GetEntries();
   if (maxEvt > 0 && maxEvt < n) n = maxEvt;
   long nSeg = 0, nHits = 0, nEvtNonEmpty = 0;
   (void)fwdDeg;

   for (long ie = 0; ie < n; ++ie) {
      tree->GetEntry(ie);
      if (!evtArr) continue;
      for (int ke = 0; ke < evtArr->GetEntriesFast(); ++ke) {
         auto *ev = (AtEvent *)evtArr->At(ke);
         if (!ev) continue;
         const auto &hits = ev->GetHits();
         if (!hits.empty()) nEvtNonEmpty++;
         nHits += hits.size();
         std::vector<WPt> pts;
         pts.reserve(hits.size());
         for (const auto &h : hits) {
            auto pos = h->GetPosition();
            pts.push_back({pos.X(), pos.Y(), pos.Z(), h->GetCharge()});
            if (pos.Z() < zmin) zmin = pos.Z();
            if (pos.Z() > zmax) zmax = pos.Z();
            hQ.Fill(h->GetCharge());
         }
         const int N = pts.size();
         if (N < Nmin + 1) continue;

         for (int i = 0; i < N; ++i) {
            std::vector<int> nb;
            for (int j = 0; j < N; ++j) {
               double dx = pts[j].x - pts[i].x, dy = pts[j].y - pts[i].y, dz = pts[j].z - pts[i].z;
               if (dx * dx + dy * dy + dz * dz <= R * R) nb.push_back(j);
            }
            if ((int)nb.size() < Nmin) continue;
            double cx = 0, cy = 0, cz = 0;
            for (int j : nb) { cx += pts[j].x; cy += pts[j].y; cz += pts[j].z; }
            double m = nb.size();
            cx /= m; cy /= m; cz /= m;
            double sxx=0,syy=0,szz=0,sxy=0,sxz=0,syz=0;
            for (int j : nb) {
               double ax = pts[j].x - cx, ay = pts[j].y - cy, az = pts[j].z - cz;
               sxx+=ax*ax; syy+=ay*ay; szz+=az*az; sxy+=ax*ay; sxz+=ax*az; syz+=ay*az;
            }
            TMatrixDSym C(3);
            C(0,0)=sxx/m; C(1,1)=syy/m; C(2,2)=szz/m;
            C(0,1)=C(1,0)=sxy/m; C(0,2)=C(2,0)=sxz/m; C(1,2)=C(2,1)=syz/m;
            TMatrixDSymEigen eig(C);
            TVectorD ev3 = eig.GetEigenValues();   // descending: col0=axis, col1/col2=transverse
            TMatrixD evec = eig.GetEigenVectors();
            double l1 = ev3[0], l2 = ev3[1], l3 = ev3[2];
            if (l2 < 0) l2 = 0;
            if (l3 < 0) l3 = 0;
            if (l1 <= 0 || l1 < linCut * l2) continue;   // principal axis must dominate -> track-like
            // classify the two transverse eigenvectors by |z-component|:
            //   most in-plane (small |z|) -> pad-plane spread (transverse diffusion)
            //   most drift-aligned (large |z|) -> drift spread (longitudinal diffusion)
            double z2 = std::fabs(evec(2, 1)), z3 = std::fabs(evec(2, 2));
            double lamXY, lamZ;
            if (z2 <= z3) { lamXY = l2; lamZ = l3; } else { lamXY = l3; lamZ = l2; }
            double zc = pts[i].z;
            nSeg++;
            pXY.Fill(zc, lamXY);
            pZd.Fill(zc, lamZ);
            hWxy.Fill(std::sqrt(lamXY));
            double vz = std::fabs(evec(2, 0));
            hTh.Fill(std::acos(std::min(1.0, vz)) * TMath::RadToDeg());
         }
      }
   }

   printf("\n===== %s =====\n", tag.c_str());
   printf("events used: %ld   z-range: [%.1f, %.1f] mm\n", n, zmin, zmax);
   printf("total hits: %ld   hits/nonempty-evt: %.1f   track-like segments: %ld\n",
          nHits, nEvtNonEmpty > 0 ? (double)nHits / nEvtNonEmpty : 0.0, nSeg);

   // sigma^2 = a + b*z ; CoefT = b[mm^2/cm]*vDrift/200, with b[mm^2/mm]*10 = b[mm^2/cm]
   TF1 linXY("linXY", "[0]+[1]*x", 100, 1000);
   TF1 linZ("linZ", "[0]+[1]*x", 100, 1000);
   pXY.Fit(&linXY, "QR");
   pZd.Fit(&linZ, "QR");
   double axy = linXY.GetParameter(0), bxy = linXY.GetParameter(1);
   double az = linZ.GetParameter(0), bz = linZ.GetParameter(1);
   double coefT = std::fabs(bxy) * 10.0 * vDrift / 200.0;
   double coefL = std::fabs(bz) * 10.0 * vDrift / 200.0;
   printf("PAD-PLANE (transverse diffusion):\n");
   printf("  sigma_xy^2 = %.4f + %.5f*z   intercept sigma0=%.3f mm\n", axy, bxy, axy > 0 ? std::sqrt(axy) : 0.0);
   printf("  -> CoefT = %.6f cm^2/us   (current par = 0.0009)\n", coefT);
   printf("DRIFT (longitudinal diffusion):\n");
   printf("  sigma_z^2  = %.4f + %.5f*z   intercept sigma0=%.3f mm\n", az, bz, az > 0 ? std::sqrt(az) : 0.0);
   printf("  -> CoefL = %.6f cm^2/us   (current par = 0.0009)\n", coefL);
   printf("CHARGE: mean per-hit = %.1f   |  mean sigma_xy = %.3f mm\n", hQ.GetMean(), hWxy.GetMean());

   TFile out(Form("width_%s.root", tag.c_str()), "RECREATE");
   pXY.Write(); pZd.Write(); hQ.Write(); hWxy.Write(); hTh.Write();
   out.Close();
   printf("wrote width_%s.root\n", tag.c_str());
   f->Close();
   return 0;
}
