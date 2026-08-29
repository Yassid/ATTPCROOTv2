/// @file dp_event_view_C14.C
/// @brief Look at individual BACKWARD (d,p) tracks: the pad-plane projection, the drift coordinate
/// against hit order, and the local radius along the track.
///
///   root -b -q 'dp_event_view_C14.C()'
///
/// Written after two hypotheses about the -20 % backward energy bias failed on summary statistics
/// alone. The pattern circle is 5-8 % small for these tracks, a first-arc fit does not recover it,
/// and the summed pad-plane path implies ~9 turns where the kinematics allow ~2.5 -- which says the
/// hits are not lying on a clean trajectory. That has to be looked at, not averaged.

#include <algorithm>
#include <vector>

void dp_event_view_C14(TString simFile = "/mnt/f/a1954_C14dp_hf/sims_b285/gs_s8001_sim.root",
                       TString recoFile = "/mnt/f/a1954_C14dp_matfx/ab_reco.root", Double_t bField = 2.85,
                       Int_t nShow = 6, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);
   const double u = 931.49401, mp = 1.007825 * u, ZPAD = 1000.0;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   TFile *fr = TFile::Open(recoFile);
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);

   auto *cv = new TCanvas("evv", "evv", 1600, 950);
   cv->Divide(3, 2);
   int shown = 0;
   Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
   for (Long64_t i = 0; i < N && shown < nShow; ++i) {
      ts->GetEntry(i);
      if (!mc) continue;
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1) continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pm = std::sqrt(px * px + py * py + pz * pz);
         if (pm <= 0) continue;
         keT = std::sqrt(pm * pm + mp * mp) - mp;
         thT = std::acos(pz / pm) * TMath::RadToDeg();
         break;
      }
      if (keT <= 0 || thT < 95 || thT > 125) continue;
      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      const AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = &t; }
      if (!best || nb < 100) continue;

      std::vector<std::array<double, 3>> P;
      for (const auto &h : best->GetHitArray()) {
         auto q = h->GetPosition();
         P.push_back({q.X(), q.Y(), ZPAD - q.Z()});
      }
      std::sort(P.begin(), P.end(), [](const std::array<double, 3> &a, const std::array<double, 3> &c) {
         return a[2] > c[2]; // backward: vertex at the highest z_lab
      });
      double pTrue = std::sqrt(keT * keT + 2 * keT * mp) / 1000.0;
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField);

      cv->cd(++shown);
      gPad->SetLeftMargin(0.14);
      auto *g = new TGraph();
      for (auto &p : P) g->SetPoint(g->GetN(), p[0], p[1]);
      g->SetMarkerStyle(6);
      g->SetMarkerColor(kAzure + 2);
      g->SetTitle(TString::Format("evt %lld: #theta=%.0f#circ, KE=%.1f MeV, %zu hits, R_{true}=%.0f mm;x [mm];y [mm]",
                                  i, thT, keT, P.size(), rTrue));
      g->Draw("ap");
      // the truth circle, centred on the best-fit centre of the hits for comparison of SIZE
      double cx = 0, cy = 0;
      for (auto &p : P) { cx += p[0]; cy += p[1]; }
      cx /= P.size(); cy /= P.size();
      auto *cir = new TEllipse(cx, cy, rTrue, rTrue);
      cir->SetFillStyle(0);
      cir->SetLineColor(kRed + 1);
      cir->SetLineWidth(2);
      cir->Draw();
      auto *tx = new TLatex();
      tx->SetNDC();
      tx->SetTextSize(0.038);
      tx->SetTextColor(kRed + 1);
      tx->DrawLatex(0.17, 0.84, "red: circle of the TRUE radius");
      printf("  evt %6lld  theta %5.1f  KE %6.2f  hits %5zu  R_true %6.1f mm  R_geo %6.1f  ratio %.3f\n", i, thT, keT,
             P.size(), rTrue, best->GetGeoRadius(), best->GetGeoRadius() / rTrue);
      ++shown;
      --shown; // one pad per event
   }
   cv->SaveAs(outDir + "dp_event_view.png");
   printf("\n  wrote %sdp_event_view.png\n\n", outDir.Data());
   fs->Close(); fr->Close();
}
