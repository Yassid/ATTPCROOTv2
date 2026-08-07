// match_diag_prog.cpp  (standalone, links libAtData -> no rootcling/ACLiC)
//
// Sim-vs-data matching diagnostic for a1954 14C(p,p'): track occupancy, per-hit charge (-> Gain)
// and transverse/longitudinal spread vs drift z (-> CoefT / CoefL).
//
// Adapted from ../../16C_dp_gnn/diagnostics/track_width_prog.cpp. The ONE structural change:
// it reads AtPatternEvent (PRA tracks) rather than AtEventH, because the surviving a1954 data
// reco files (~/a1954_C14_fit_300torr/in/) are the GATED ones and carry AtPatternEvent only.
// Reading tracks on both sides also makes the comparison apples-to-apples: the sim's truth-gated
// simg_reco.root is the analogue of the data's IC+PID-gated file.
//
// The diffusion numbers are MEASURED, not scanned: sigma^2(z) = sigma0^2 + 2*D*z/v_drift, so a
// straight-line fit of sigma^2 vs z gives CoefT/CoefL in the par-file units directly. Run it on
// the data and copy the answer into the sim par.
//
// Build: ./build_match.sh    Run: ./match_diag <file.root> <tag> [vDrift] [maxEvt] [R] [Nmin] [linCut]
#include "AtHit.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"

#include <TClonesArray.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

struct WPt {
   double x, y, z, q;
};

int main(int argc, char **argv)
{
   if (argc < 3) {
      printf("usage: %s <file.root> <tag> [vDrift=1.30] [maxEvt=-1] [R=8] [Nmin=6] [linCut=3]\n", argv[0]);
      return 1;
   }
   std::string fname = argv[1], tag = argv[2];
   double vDrift = argc > 3 ? atof(argv[3]) : 1.30;
   long maxEvt = argc > 4 ? atol(argv[4]) : -1;
   double R = argc > 5 ? atof(argv[5]) : 8.0;
   int Nmin = argc > 6 ? atoi(argv[6]) : 6;
   double linCut = argc > 7 ? atof(argv[7]) : 3.0;

   TFile *f = TFile::Open(fname.c_str());
   if (!f || f->IsZombie()) {
      printf("ERROR: cannot open %s\n", fname.c_str());
      return 1;
   }
   TTree *tree = (TTree *)f->Get("cbmsim");
   if (!tree || !tree->GetBranch("AtPatternEvent")) {
      printf("ERROR: no cbmsim/AtPatternEvent in %s\n", fname.c_str());
      return 1;
   }
   tree->SetBranchStatus("*", 0);
   tree->SetBranchStatus("AtPatternEvent*", 1);
   TClonesArray *peArr = nullptr;
   tree->SetBranchAddress("AtPatternEvent", &peArr);

   double zmin = 1e9, zmax = -1e9;
   TProfile pXY(Form("p_xy_%s", tag.c_str()), Form("%s;z [mm];#sigma_{xy}^{2} [mm^{2}]", tag.c_str()), 50, 0, 1000);
   TProfile pZd(Form("p_zd_%s", tag.c_str()), Form("%s;z [mm];#sigma_{z}^{2} [mm^{2}]", tag.c_str()), 50, 0, 1000);
   TH1D hQ(Form("h_q_%s", tag.c_str()), Form("%s per-hit charge;charge;hits", tag.c_str()), 300, 0, 6000);
   TH1D hWxy(Form("h_wxy_%s", tag.c_str()), Form("%s #sigma_{xy};#sigma_{xy} [mm];segments", tag.c_str()), 100, 0, 8);
   TH1D hNpt(Form("h_npt_%s", tag.c_str()), Form("%s hits per track;hits;tracks", tag.c_str()), 100, 0, 500);
   TH1D hNtr(Form("h_ntr_%s", tag.c_str()), Form("%s tracks per event;tracks;events", tag.c_str()), 20, 0, 20);
   TH1D hQtr(Form("h_qtr_%s", tag.c_str()), Form("%s total charge per track;charge;tracks", tag.c_str()), 200, 0, 4e5);

   long n = tree->GetEntries();
   if (maxEvt > 0 && maxEvt < n)
      n = maxEvt;
   long nSeg = 0, nHits = 0, nTrk = 0, nEvtNonEmpty = 0;
   std::vector<double> qAll, nptAll, lenAll, densAll, qmmAll;

   for (long ie = 0; ie < n; ++ie) {
      tree->GetEntry(ie);
      if (!peArr || peArr->GetEntriesFast() == 0)
         continue;
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe)
         continue;
      auto &cands = pe->GetTrackCand();
      hNtr.Fill(cands.size());
      if (!cands.empty())
         nEvtNonEmpty++;

      for (auto &trk : cands) {
         // GetHitArray() holds POINTERS -- never copy the AtHit objects out of it
         const auto &hitArr = trk.GetHitArray();
         if (hitArr.empty())
            continue;
         ++nTrk;
         nHits += hitArr.size();
         hNpt.Fill(hitArr.size());
         nptAll.push_back(hitArr.size());

         std::vector<WPt> pts;
         pts.reserve(hitArr.size());
         double qTrk = 0;
         for (const auto &h : hitArr) {
            if (!h)
               continue;
            auto pos = h->GetPosition();
            double q = h->GetCharge();
            pts.push_back({pos.X(), pos.Y(), pos.Z(), q});
            zmin = std::min(zmin, pos.Z());
            zmax = std::max(zmax, pos.Z());
            hQ.Fill(q);
            qAll.push_back(q);
            qTrk += q;
         }
         hQtr.Fill(qTrk);

         // track length = max pairwise separation (chord). hits/track alone is not comparable
         // between sim and data unless the tracks are the same length, so normalise by it.
         double len = 0;
         for (size_t a = 0; a < pts.size(); ++a)
            for (size_t b = a + 1; b < pts.size(); ++b) {
               double dx = pts[a].x - pts[b].x, dy = pts[a].y - pts[b].y, dz = pts[a].z - pts[b].z;
               len = std::max(len, dx * dx + dy * dy + dz * dz);
            }
         len = std::sqrt(len);
         if (len > 1) {
            lenAll.push_back(len);
            densAll.push_back(pts.size() / len);
            qmmAll.push_back(qTrk / len);
         }

         const int N = pts.size();
         if (N < Nmin + 1)
            continue;
         for (int i = 0; i < N; ++i) {
            std::vector<int> nb;
            for (int j = 0; j < N; ++j) {
               double dx = pts[j].x - pts[i].x, dy = pts[j].y - pts[i].y, dz = pts[j].z - pts[i].z;
               if (dx * dx + dy * dy + dz * dz <= R * R)
                  nb.push_back(j);
            }
            if ((int)nb.size() < Nmin)
               continue;
            double cx = 0, cy = 0, cz = 0;
            for (int j : nb) {
               cx += pts[j].x;
               cy += pts[j].y;
               cz += pts[j].z;
            }
            double m = nb.size();
            cx /= m;
            cy /= m;
            cz /= m;
            double sxx = 0, syy = 0, szz = 0, sxy = 0, sxz = 0, syz = 0;
            for (int j : nb) {
               double ax = pts[j].x - cx, ay = pts[j].y - cy, az = pts[j].z - cz;
               sxx += ax * ax; syy += ay * ay; szz += az * az;
               sxy += ax * ay; sxz += ax * az; syz += ay * az;
            }
            TMatrixDSym C(3);
            C(0, 0) = sxx / m; C(1, 1) = syy / m; C(2, 2) = szz / m;
            C(0, 1) = C(1, 0) = sxy / m;
            C(0, 2) = C(2, 0) = sxz / m;
            C(1, 2) = C(2, 1) = syz / m;
            TMatrixDSymEigen eig(C);
            TVectorD ev3 = eig.GetEigenValues();
            TMatrixD evec = eig.GetEigenVectors();
            double l1 = ev3[0], l2 = std::max(0.0, ev3[1]), l3 = std::max(0.0, ev3[2]);
            if (l1 <= 0 || l1 < linCut * l2)
               continue; // must be track-like locally
            double z2 = std::fabs(evec(2, 1)), z3 = std::fabs(evec(2, 2));
            double lamXY = (z2 <= z3) ? l2 : l3;
            double lamZ = (z2 <= z3) ? l3 : l2;
            ++nSeg;
            pXY.Fill(pts[i].z, lamXY);
            pZd.Fill(pts[i].z, lamZ);
            hWxy.Fill(std::sqrt(lamXY));
         }
      }
   }

   auto median = [](std::vector<double> &v) {
      if (v.empty()) return 0.0;
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };

   printf("\n===== %s =====\n", tag.c_str());
   printf("events read: %ld   non-empty: %ld   z-range [%.1f, %.1f] mm\n", n, nEvtNonEmpty, zmin, zmax);
   printf("tracks: %ld   tracks/evt: %.2f   hits: %ld   hits/track: %.1f (median %.0f)\n", nTrk,
          nEvtNonEmpty ? (double)nTrk / nEvtNonEmpty : 0.0, nHits, nTrk ? (double)nHits / nTrk : 0.0, median(nptAll));
   // means from the vectors, NOT from the histograms: the TH1 ranges clip the high-charge tail
   // and TH1::GetMean() drops the overflow, which biased the first pass low.
   auto vmean = [](const std::vector<double> &v) {
      double s = 0;
      for (double x : v) s += x;
      return v.empty() ? 0.0 : s / v.size();
   };
   printf("CHARGE/hit: mean %.1f  median %.1f   | charge/track: mean %.0f\n", vmean(qAll), median(qAll),
          vmean(qmmAll) * vmean(lenAll));
   printf("LENGTH    : median %.1f mm   | hits/mm: median %.3f   | charge/mm: median %.1f\n", median(lenAll),
          median(densAll), median(qmmAll));

   TF1 linXY("linXY", "[0]+[1]*x", 100, 1000);
   TF1 linZ("linZ", "[0]+[1]*x", 100, 1000);
   pXY.Fit(&linXY, "QR");
   pZd.Fit(&linZ, "QR");
   double axy = linXY.GetParameter(0), bxy = linXY.GetParameter(1);
   double az = linZ.GetParameter(0), bz = linZ.GetParameter(1);
   // sigma^2 [mm^2] vs z [mm]; CoefT[cm^2/us] = b[mm^2/mm] * 10 * vDrift / 200
   printf("PAD-PLANE : sigma_xy^2 = %.4f + %.6f*z   sigma0 = %.3f mm   -> CoefT = %.6f cm^2/us\n", axy, bxy,
          axy > 0 ? std::sqrt(axy) : 0.0, std::fabs(bxy) * 10.0 * vDrift / 200.0);
   printf("DRIFT     : sigma_z^2  = %.4f + %.6f*z   sigma0 = %.3f mm   -> CoefL = %.6f cm^2/us\n", az, bz,
          az > 0 ? std::sqrt(az) : 0.0, std::fabs(bz) * 10.0 * vDrift / 200.0);
   printf("mean sigma_xy = %.3f mm   (%ld local segments)\n", hWxy.GetMean(), nSeg);

   TFile out(Form("match_%s.root", tag.c_str()), "RECREATE");
   pXY.Write(); pZd.Write(); hQ.Write(); hWxy.Write(); hNpt.Write(); hNtr.Write(); hQtr.Write();
   out.Close();
   printf("wrote match_%s.root\n", tag.c_str());
   f->Close();
   return 0;
}
