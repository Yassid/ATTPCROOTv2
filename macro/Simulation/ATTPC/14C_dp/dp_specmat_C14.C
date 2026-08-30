/// @file dp_specmat_C14.C
/// @brief The SpecMAT (full+partial)-track plane for 14C(d,p)15C: E_x from hit-cloud geometry,
///        no helix fit. Replaces the HELIOS/d_max validation plot, which measured the wrong thing.
///
///   root -b -q 'dp_specmat_C14.C("/mnt/f/a1954_C14dp_sm","/mnt/f/a1954_C14dp/sims_b285",2.85)'
///
/// WHY THIS AND NOT A HELIOS RETURN PLOT. On a1975 both observables were measured against the
/// trajectory integrator on the same hit clouds and they are six times apart:
///     d_max (max excursion from the beam axis)   bias 1.007   spread +-5 %
///     z at the first return (the HELIOS quantity) bias 0.848   spread +-30 %
/// That is structural, not accidental: d_max is an ENVELOPE, so it does not care where the track
/// started, how much of the inner arc is missing, or where the endpoint is called. The return z
/// depends on all three. SpecMAT therefore transfers through d_max, NOT through the return
/// position -- see [[reference_specmat_track_length_method]].
///
/// THE DELIVERABLE IS PANEL A: the (d_max, theta_lab) plane with the three levels overlaid.
/// The states are supposed to separate into BANDS there. One level cannot show separation, which
/// is why all three are read here.
///
/// TWO TRAPS THIS MACRO IS BUILT TO AVOID, both paid for on a1975:
///   * Draw the plane as a SCATTER framed on the data. A COLZ map with bins derived from sqrt(N)
///     hid the bands completely and led to the plane being called featureless on 348 tracks.
///   * d_max/2 is NOT the vertex radius. The ejectile loses energy, so the trajectory is a
///     TIGHTENING spiral and by the apex the radius is already below its vertex value. Compare
///     d_max against the TRUE rMax from the MC trajectory (AtTpcPoint), not against
///     p_T/(0.3B) from the vertex momentum, and never paper over the difference with a fudge.
///
/// theta is estimator B (SpecMAT's own L_dmax,Z with the off-axis start put back): the swept angle
/// to the apex is GEOMETRY, psi = pi - 2*asin(d0/2R), not something to measure. Only quantities
/// local to the first half turn enter, so truncation beyond the apex cannot touch it. On a1975 it
/// beat the measured-winding estimator in BOTH channels (0.972 vs 0.991 on (d,t), 0.984 vs 1.262
/// on (p,d)), and the measured-winding one is biased with a channel-dependent sign.

#include <algorithm>
#include <vector>

static double sm_q(std::vector<double> v, double p)
{
   if (v.size() < 5)
      return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_specmat_C14(TString smDir = "/mnt/f/a1954_C14dp_sm", TString simDir = "/mnt/f/a1954_C14dp/sims_b285",
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

   // per level: the measured plane, the truth plane, and the two closure ratios
   std::vector<double> gdm[3], gth[3], tdm[3], tth[3];
   std::vector<double> rDm[3], rTh[3];

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
      ts->SetBranchAddress("AtTpcPoint", &pts);

      Long64_t N = std::min(tr->GetEntries(), ts->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         if (!mc)
            continue;
         // truth: the primary proton
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
         // TRUE rMax: the largest excursion of the MC trajectory from the beam axis.
         // This is the quantity d_max estimates -- NOT p_T/(0.3B), which is the vertex radius and
         // is systematically larger because the spiral tightens as the proton slows.
         double dmaxTrue = 0; // max excursion from the beam axis = the quantity d_max measures
         if (pts)
            for (int k = 0; k < pts->GetEntriesFast(); ++k) {
               // NOTE: the branch is called "AtTpcPoint" but the objects in it are AtMCPoint.
               // AtTpcPoint has no header or source anywhere in the repo -- only a stale
               // \pragma link line in AtDetectors/AtTpc/AtTpcLinkDef.h -- so casting to it does
               // not even compile. Always ask the TClonesArray what it holds.
               auto *q = (AtMCPoint *)pts->At(k);
               if (!q || q->GetTrackID() != protonId)
                  continue;
               dmaxTrue = std::max(dmaxTrue, 10.0 * std::hypot(q->GetX(), q->GetY())); // cm -> mm
            }
         if (dmaxTrue <= 0)
            continue;

         // ---- measured, from the hit cloud only ----
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
         // d_max: the apex, the largest distance from the beam axis. Complete at HALF a turn,
         // which is what makes the method work on partial tracks.
         size_t iApex = 0;
         double dmax = -1;
         for (size_t k = 0; k < X.size(); ++k) {
            double d = std::hypot(X[k], Y[k]);
            if (d > dmax) { dmax = d; iApex = k; }
         }
         // the start is the END OF THE TRACK NEAREST THE BEAM AXIS -- a fit-free definition, and
         // the one a1975 used. d0 is how far off-axis the observed track actually begins.
         double dFront = std::hypot(X.front(), Y.front()), dBack = std::hypot(X.back(), Y.back());
         size_t iStart = (dFront <= dBack) ? 0 : X.size() - 1;
         double d0 = std::min(dFront, dBack);
         double R = 0.5 * dmax;
         if (!(R > 0) || d0 >= 2 * R)
            continue;
         // estimator B: swept angle to the apex is geometry, not a measurement
         double psi = TMath::Pi() - 2.0 * std::asin(d0 / (2 * R));
         double dz = Z[iApex] - Z[iStart];
         if (std::fabs(dz) < 1e-6)
            continue;
         double thG = std::atan2(psi * R, dz) * R2D; // signed dz keeps the hemisphere
         if (thG < 0)
            thG += 180.0;
         if (thG < thMin)
            continue;

         gdm[L].push_back(dmax);
         gth[L].push_back(thG);
         tdm[L].push_back(dmaxTrue); // same quantity, same units, no factor of 2
         tth[L].push_back(thT);
         rDm[L].push_back(dmax / dmaxTrue);
         rTh[L].push_back(thG / thT);
      }
      printf("  %-7s  %5zu tracks\n", levs[L], gdm[L].size());
      fr->Close();
      fs->Close();
   }

   // ---------------- draw ----------------
   auto *c = new TCanvas("csm", "specmat", 1500, 980);
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

   // A: the measured plane -- SCATTER, framed on the data. This is the deliverable.
   c->cd(1);
   gPad->SetGrid();
   double x0, x1, y0, y1;
   frameOn(gth, gdm, x0, x1, y0, y1);
   auto *fA = gPad->DrawFrame(x0, y0, x1, y1,
                              "A  measured plane: d_{max} vs #theta_{lab}, from hits only"
                              ";#theta_{lab} geometric [deg];d_{max} [mm]");
   fA->GetXaxis()->SetTitleSize(0.042);
   fA->GetYaxis()->SetTitleSize(0.042);
   auto *lgA = new TLegend(0.62, 0.70, 0.98, 0.92);
   lgA->SetBorderSize(0);
   lgA->SetFillStyle(0);
   for (int L = 0; L < 3; ++L) {
      if (gdm[L].empty()) continue;
      auto *g = new TGraph(gdm[L].size(), gth[L].data(), gdm[L].data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.35);
      g->SetMarkerColor(col[L]);
      g->Draw("P same");
      lgA->AddEntry(g, TString::Format("E_{x} = %.3f MeV  (n=%zu)", levEx[L], gdm[L].size()), "p");
   }
   lgA->Draw();

   // B: the same plane from TRUTH -- separates "the bands do not exist" from "the method smears them"
   c->cd(2);
   gPad->SetGrid();
   frameOn(tth, tdm, x0, x1, y0, y1);
   auto *fB = gPad->DrawFrame(x0, y0, x1, y1,
                              "B  the same plane from TRUTH (d_{max} of the MC trajectory)"
                              ";#theta_{lab} true [deg];d_{max} true [mm]");
   fB->GetXaxis()->SetTitleSize(0.042);
   fB->GetYaxis()->SetTitleSize(0.042);
   for (int L = 0; L < 3; ++L) {
      if (tdm[L].empty()) continue;
      auto *g = new TGraph(tdm[L].size(), tth[L].data(), tdm[L].data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.35);
      g->SetMarkerColor(col[L]);
      g->Draw("P same");
   }

   // C: closure on d_max. a1975 got 1.007 +- 5 % against the integrator's rMax.
   c->cd(3);
   gPad->SetGrid();
   auto *hD = new TH1D("hD", "C  closure: d_{max} / d_{max}^{true};d_{max} / d_{max}^{true};tracks", 90, 0.85, 1.15);
   for (int L = 0; L < 3; ++L)
      for (double v : rDm[L]) hD->Fill(v);
   hD->SetLineColor(kAzure + 2);
   hD->SetLineWidth(2);
   hD->Draw("hist");

   // D: closure on the angle estimator
   c->cd(4);
   gPad->SetGrid();
   auto *hT = new TH1D("hT", "D  closure: #theta_{geometric} / #theta_{true};#theta_{geo}/#theta_{true};tracks",
                       90, 0.7, 1.3);
   for (int L = 0; L < 3; ++L)
      for (double v : rTh[L]) hT->Fill(v);
   hT->SetLineColor(kOrange + 7);
   hT->SetLineWidth(2);
   hT->Draw("hist");

   TString out = outDir + "dp_specmat.png";
   c->SaveAs(out);

   std::vector<double> aD, aT;
   for (int L = 0; L < 3; ++L) {
      aD.insert(aD.end(), rDm[L].begin(), rDm[L].end());
      aT.insert(aT.end(), rTh[L].begin(), rTh[L].end());
   }
   printf("\n  d_max / d_max(true) : median %.3f   IQR/1.349 %.3f\n", sm_q(aD, .5),
          (sm_q(aD, .75) - sm_q(aD, .25)) / 1.349);
   printf("  theta_geo / theta_true: median %.3f   IQR/1.349 %.3f\n", sm_q(aT, .5),
          (sm_q(aT, .75) - sm_q(aT, .25)) / 1.349);
   printf("  wrote %s\n\n", out.Data());
}
