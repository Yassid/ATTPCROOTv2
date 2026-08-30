/// @file dp_helios_C14.C
/// @brief The HELIOS plane for 14C(d,p)15C: ejectile energy against the axial return position,
///        both taken from hit geometry with no helix fit.
///
///   root -b -q 'dp_helios_C14.C("/mnt/f/a1954_C14dp_sm","/mnt/f/a1954_C14dp/sims_b285",2.85)'
///
/// In a solenoid a particle emitted from the axis returns to the axis after one full cyclotron
/// period, at an axial distance z_ret = v_parallel * T_cyc. Because T_cyc is the same for every
/// ejectile of a given species, E_lab is LINEAR in z_ret and states of different Q value fall on
/// PARALLEL LINES. That is the whole idea of the method, and it is what panel A tests.
///
/// WHY THIS IS WORTH REVISITING. On a1975 the return position measured badly -- bias 0.848 with
/// +-30 % scatter, against +-5 % for d_max -- and the conclusion drawn there was to use d_max.
/// But the CAUSE of that scatter was specific: those reconstructed tracks began a median 29.6 mm
/// OFF the beam axis, so the return position needed vertex anchoring plus an inner-arc
/// correction, and d0 varied 27-71 mm track to track. HELIOS assumes emission FROM the axis, and
/// in this simulation the vertex is on the axis by construction. The a1975 number therefore does
/// not automatically carry over, and the honest thing is to measure it here rather than inherit
/// it. Panel C is that measurement.
///
/// z_ret is taken directly: walk past the apex (the largest excursion from the axis) and find
/// where the track next comes closest to the axis. That is the return, and it is measured, not
/// inferred from the apex -- using 2*dz_apex instead would make this plot a restatement of the
/// SpecMAT one rather than an independent check.

#include <algorithm>
#include <vector>

static double hl_q(std::vector<double> v, double p)
{
   if (v.size() < 5)
      return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_helios_C14(TString smDir = "/mnt/f/a1954_C14dp_sm", TString simDir = "/mnt/f/a1954_C14dp/sims_b285",
                   Double_t bField = 2.85, Double_t thMin = 92.0, Int_t minHits = 120, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   const double u = 931.49401, mp = 1.007825 * u, ZPAD = 1000.0, R2D = 57.29577951;
   const char *levs[3] = {"gs", "ex0740", "ex3103"};
   const double levEx[3] = {0.0, 0.740, 3.103};
   const int col[3] = {kAzure + 2, kOrange + 7, kGreen + 3};

   std::vector<double> gz[3], gke[3], tz[3], tke[3], rZ[3], rKE[3];

   for (int L = 0; L < 3; ++L) {
      TString rf = smDir + "/" + levs[L] + "/reco.root";
      // GLOB for the sim rather than hardcoding a seed: the seed differs per field
      // (s800x at 2.85 T, s810x at 4 T, s820x at 7 T), and a hardcoded one silently reports
      // every sample as missing when the macro is pointed at another field.
      TString sf = gSystem->GetFromPipe("ls " + simDir + "/" + TString(levs[L]) + "_*_sim.root 2>/dev/null | head -1");
      sf = sf.Strip(TString::kBoth);
      TFile *fr = TFile::Open(rf);
      TFile *fs = TFile::Open(sf);
      if (!fr || fr->IsZombie() || !fs || fs->IsZombie()) {
         printf("  MISSING %s or %s\n", rf.Data(), sf.Data());
         continue;
      }
      TTree *tr = (TTree *)fr->Get("cbmsim");
      TTree *ts = (TTree *)fs->Get("cbmsim");
      TClonesArray *pe = nullptr, *mc = nullptr, *pts = nullptr;
      tr->SetBranchAddress("AtPatternEvent", &pe);
      ts->SetBranchAddress("MCTrack", &mc);
      ts->SetBranchAddress("AtTpcPoint", &pts); // NOTE: branch name lies, the objects are AtMCPoint

      Long64_t N = std::min(tr->GetEntries(), ts->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         if (!mc)
            continue;
         double keT = -1, thT = -1;
         int protonId = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pm = std::sqrt(px * px + py * py + pz * pz);
            if (pm <= 0)
               break;
            keT = std::sqrt(pm * pm + mp * mp) - mp;
            thT = std::acos(pz / pm) * R2D;
            protonId = k;
            break;
         }
         if (keT <= 0 || thT < thMin)
            continue;

         // TRUE return: walk the MC trajectory, find the apex, then the next closest approach to
         // the axis after it. Points come in order along the track.
         double zRetTrue = -1;
         if (pts) {
            std::vector<double> px_, py_, pz_;
            for (int k = 0; k < pts->GetEntriesFast(); ++k) {
               auto *q = (AtMCPoint *)pts->At(k);
               if (!q || q->GetTrackID() != protonId)
                  continue;
               px_.push_back(10.0 * q->GetX());
               py_.push_back(10.0 * q->GetY());
               pz_.push_back(10.0 * q->GetZ());
            }
            if (px_.size() > 8) {
               size_t ia = 0;
               double dm = -1;
               for (size_t k = 0; k < px_.size(); ++k) {
                  double d = std::hypot(px_[k], py_[k]);
                  if (d > dm) { dm = d; ia = k; }
               }
               const size_t W = 3;
               size_t ir = 0;
               for (size_t k = ia + W; k + W < px_.size(); ++k) {
                  double d = std::hypot(px_[k], py_[k]);
                  bool isMin = true;
                  for (size_t j = 1; j <= W; ++j)
                     if (std::hypot(px_[k - j], py_[k - j]) < d || std::hypot(px_[k + j], py_[k + j]) < d) {
                        isMin = false;
                        break;
                     }
                  if (isMin) { ir = k; break; }
               }
               if (ir > ia)
                  zRetTrue = std::fabs(pz_[ir] - pz_[0]);
            }
         }
         if (zRetTrue <= 0)
            continue;

         // ---- measured, hits only ----
         tr->GetEntry(i);
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPatternEvent *)pe->At(0);
         if (!ev)
            continue;
         AtTrack *best = nullptr;
         size_t nb = 0;
         for (auto &t : ev->GetTrackCand())
            if (t.GetHitArray().size() > nb) {
               nb = t.GetHitArray().size();
               best = const_cast<AtTrack *>(&t);
            }
         if (!best || (int)nb < minHits)
            continue;
         auto *cl = best->GetHitClusterArray();
         if (!cl || cl->size() < 12)
            continue;
         std::vector<double> X, Y, Z;
         for (auto &c : *cl) {
            auto p = c.GetPosition();
            X.push_back(p.X());
            Y.push_back(p.Y());
            Z.push_back(ZPAD - p.Z());
         }
         // order along the track: start from the end nearest the axis (fit-free, as for SpecMAT)
         if (std::hypot(X.back(), Y.back()) < std::hypot(X.front(), Y.front())) {
            std::reverse(X.begin(), X.end());
            std::reverse(Y.begin(), Y.end());
            std::reverse(Z.begin(), Z.end());
         }
         size_t iApex = 0;
         double dmax = -1;
         for (size_t k = 0; k < X.size(); ++k) {
            double d = std::hypot(X[k], Y[k]);
            if (d > dmax) { dmax = d; iApex = k; }
         }
         // the RETURN: the FIRST local minimum of the distance to the axis after the apex.
         // Taking the global minimum over the rest of the track is wrong on a multi-turn spiral,
         // which approaches the axis once per turn -- it silently returns a LATER turn (or the
         // track end when the spiral never comes back), and that alone put the measured/true
         // ratio at 1.758 with a tail out to 800 mm. Require the distance to rise again over a
         // short window so cluster noise cannot fake a minimum.
         const size_t WIN = 4;
         size_t iRet = 0;
         for (size_t k = iApex + WIN; k + WIN < X.size(); ++k) {
            double d = std::hypot(X[k], Y[k]);
            bool isMin = true;
            for (size_t j = 1; j <= WIN; ++j) {
               if (std::hypot(X[k - j], Y[k - j]) < d || std::hypot(X[k + j], Y[k + j]) < d) {
                  isMin = false;
                  break;
               }
            }
            if (isMin) { iRet = k; break; }
         }
         if (iRet <= iApex)
            continue; // truncated before it came back: there is no return to measure
         double zRet = std::fabs(Z[iRet] - Z[0]);
         if (zRet < 5)
            continue;

         // ejectile energy from the SAME geometry, so the plane is fit-free end to end
         double d0 = std::hypot(X[0], Y[0]);
         double R = 0.5 * dmax;
         if (!(R > 0) || d0 >= 2 * R)
            continue;
         double psi = TMath::Pi() - 2.0 * std::asin(d0 / (2 * R));
         double dz = Z[iApex] - Z[0];
         if (std::fabs(dz) < 1e-6)
            continue;
         double thG = std::atan2(psi * R, dz) * R2D;
         if (thG < 0)
            thG += 180.0;
         double sth = std::sin(thG / R2D);
         if (sth < 1e-3)
            continue;
         double pT = 0.299792458 * bField * R;   // MeV/c, R in mm
         double pTot = pT / sth;
         double keG = std::sqrt(pTot * pTot + mp * mp) - mp;

         gz[L].push_back(zRet);
         gke[L].push_back(keG);
         tz[L].push_back(zRetTrue);
         tke[L].push_back(keT);
         rZ[L].push_back(zRet / zRetTrue);
         rKE[L].push_back(keG / keT);
      }
      printf("  %-7s  %5zu tracks\n", levs[L], gz[L].size());
      fr->Close();
      fs->Close();
   }

   auto *c = new TCanvas("chl", "helios", 1500, 980);
   c->Divide(2, 2);
   auto frameOn = [](std::vector<double> *a, std::vector<double> *b, double &lo0, double &hi0, double &lo1,
                     double &hi1) {
      lo0 = 1e9; hi0 = -1e9; lo1 = 1e9; hi1 = -1e9;
      for (int L = 0; L < 3; ++L)
         for (size_t k = 0; k < a[L].size(); ++k) {
            lo0 = std::min(lo0, a[L][k]); hi0 = std::max(hi0, a[L][k]);
            lo1 = std::min(lo1, b[L][k]); hi1 = std::max(hi1, b[L][k]);
         }
      double p0 = 0.04 * (hi0 - lo0), p1 = 0.04 * (hi1 - lo1);
      lo0 -= p0; hi0 += p0; lo1 -= p1; hi1 += p1;
   };
   // robust frame: clip at the 1st/99th percentile so a handful of runaway points cannot
   // compress the structure into the bottom of the pad
   auto frameRobust = [](std::vector<double> *a, std::vector<double> *b, double &lo0, double &hi0,
                         double &lo1, double &hi1) {
      std::vector<double> A, B;
      for (int L = 0; L < 3; ++L) {
         A.insert(A.end(), a[L].begin(), a[L].end());
         B.insert(B.end(), b[L].begin(), b[L].end());
      }
      std::sort(A.begin(), A.end());
      std::sort(B.begin(), B.end());
      auto pc = [](std::vector<double> &v, double p) {
         return v.empty() ? 0.0 : v[(size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)))];
      };
      lo0 = pc(A, .01); hi0 = pc(A, .99); lo1 = pc(B, .01); hi1 = pc(B, .99);
      double p0 = 0.06 * (hi0 - lo0), p1 = 0.06 * (hi1 - lo1);
      lo0 -= p0; hi0 += p0; lo1 -= p1; hi1 += p1;
   };

   // one frame, taken from TRUTH, used for both scatter panels: they exist to be compared, and
   // a panel that auto-frames on its own outliers cannot be read against one that does not.
   double x0, x1, y0, y1;
   frameRobust(tz, tke, x0, x1, y0, y1);

   c->cd(1);
   gPad->SetGrid();
   auto *fA = gPad->DrawFrame(x0, y0, x1, y1,
                              "A  HELIOS plane, measured: E_{lab} vs z_{return}, hits only"
                              ";z_{return} [mm];proton KE, geometric [MeV]");
   fA->GetXaxis()->SetTitleSize(0.042);
   fA->GetYaxis()->SetTitleSize(0.042);
   auto *lg = new TLegend(0.60, 0.70, 0.98, 0.92);
   lg->SetBorderSize(0);
   lg->SetFillStyle(0);
   for (int L = 0; L < 3; ++L) {
      if (gz[L].empty()) continue;
      auto *g = new TGraph(gz[L].size(), gz[L].data(), gke[L].data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.35);
      g->SetMarkerColor(col[L]);
      g->Draw("P same");
      lg->AddEntry(g, TString::Format("E_{x} = %.3f MeV  (n=%zu)", levEx[L], gz[L].size()), "p");
   }
   {
      long nOut = 0, nTot = 0;
      for (int L = 0; L < 3; ++L)
         for (size_t k = 0; k < gz[L].size(); ++k) {
            ++nTot;
            if (gz[L][k] < x0 || gz[L][k] > x1 || gke[L][k] < y0 || gke[L][k] > y1) ++nOut;
         }
      auto *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.030);
      tx->SetTextColor(kGray + 2);
      tx->DrawLatex(0.14, 0.16,
                    TString::Format("%.1f %% of measured points fall outside this frame", 100.0 * nOut / std::max(1L, nTot)));
   }
   lg->Draw();

   c->cd(2);
   gPad->SetGrid();
   auto *fB = gPad->DrawFrame(x0, y0, x1, y1,
                              "B  the same plane from TRUTH"
                              ";z_{return} true [mm];proton KE true [MeV]");
   fB->GetXaxis()->SetTitleSize(0.042);
   fB->GetYaxis()->SetTitleSize(0.042);
   for (int L = 0; L < 3; ++L) {
      if (tz[L].empty()) continue;
      auto *g = new TGraph(tz[L].size(), tz[L].data(), tke[L].data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.35);
      g->SetMarkerColor(col[L]);
      g->Draw("P same");
   }

   c->cd(3);
   gPad->SetGrid();
   auto *hZ = new TH1D("hZ", "C  closure: z_{return} measured / true;z_{ret}/z_{ret}^{true};tracks", 100, 0.4, 1.6);
   for (int L = 0; L < 3; ++L)
      for (double v : rZ[L]) hZ->Fill(v);
   hZ->SetLineColor(kAzure + 2);
   hZ->SetLineWidth(2);
   hZ->Draw("hist");

   c->cd(4);
   gPad->SetGrid();
   auto *hK = new TH1D("hK", "D  closure: KE geometric / true;KE_{geo}/KE_{true};tracks", 100, 0.4, 1.6);
   for (int L = 0; L < 3; ++L)
      for (double v : rKE[L]) hK->Fill(v);
   hK->SetLineColor(kOrange + 7);
   hK->SetLineWidth(2);
   hK->Draw("hist");

   TString out = outDir + "dp_helios.png";
   c->SaveAs(out);

   std::vector<double> aZ, aK;
   for (int L = 0; L < 3; ++L) {
      aZ.insert(aZ.end(), rZ[L].begin(), rZ[L].end());
      aK.insert(aK.end(), rKE[L].begin(), rKE[L].end());
   }
   printf("\n  z_return meas/true : median %.3f   IQR/1.349 %.3f\n", hl_q(aZ, .5),
          (hl_q(aZ, .75) - hl_q(aZ, .25)) / 1.349);
   printf("  KE geometric/true  : median %.3f   IQR/1.349 %.3f\n", hl_q(aK, .5),
          (hl_q(aK, .75) - hl_q(aK, .25)) / 1.349);
   printf("  (a1975 measured the return position at 0.848 +- 30 %%, but with tracks starting\n"
          "   29.6 mm OFF axis; here the vertex is on the axis by construction)\n");
   printf("  wrote %s\n\n", out.Data());
}
