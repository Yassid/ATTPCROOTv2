/// @file dp_view7T_C14.C
/// @brief Individual 7 T tracks in the band where the fit fails: theta_lab 95-125 deg.
///
///   root -b -q 'dp_view7T_C14.C("/mnt/f/a1954_C14dp_hits7T/attpc","7 T, AT-TPC pads")'
///
/// Per-10-degree plots put the 7 T failure in a band: efficiency 0.09-0.15 there, implied hit
/// sigma off scale, resolution 5-20 %. Those are the 3-6 MeV backward protons whose helix radius
/// at 7 T is ~3.6 cm -- the tightest, most wound spirals in the matrix. This draws them.
///
/// Each panel is one event: the hit cloud in the pad plane (blue), the circle of the TRUE radius
/// centred on the cloud (red), and -- where the fit produced one -- the FITTED trajectory (green).
/// If the fit is diverging, the green will leave the blue. The header carries the truth, the fit
/// and chi2/ndf, so a panel can be read without cross-referencing anything.

#include <algorithm>
#include <vector>

void dp_view7T_C14(TString stem = "/mnt/f/a1954_C14dp_hits7T/attpc", TString title = "7 T, AT-TPC pads",
                   TString simFile = "/mnt/f/a1954_C14dp/sims_b700/gs_s8201_sim.root", Double_t bField = 7.0,
                   Double_t thLo = 95, Double_t thHi = 125, Int_t nShow = 6, TString outDir = "")
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
   TFile *fr = TFile::Open(stem + "_reco.root");
   if (!fr || fr->IsZombie()) { printf("\033[1;31mno reco: %s_reco.root\033[0m\n", stem.Data()); return; }
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);
   TFile *ff = TFile::Open(stem + "_genfit.root");
   TTree *tf = ff && !ff->IsZombie() ? (TTree *)ff->Get("cbmsim") : nullptr;
   TClonesArray *te = nullptr;
   if (tf) tf->SetBranchAddress("AtTrackingEvent", &te);

   auto *cv = new TCanvas("v7", "v7", 1650, 1000);
   cv->Divide(3, 2);
   int shown = 0;
   Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
   printf("\n  %s -- tracks with true theta_lab in %.0f-%.0f deg\n", title.Data(), thLo, thHi);
   printf("  %8s %8s %8s | %8s %8s %10s | %6s %8s\n", "event", "th_true", "KE_true", "th_fit", "KE_fit", "chi2/ndf",
          "hits", "R_true");
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
      if (keT <= 0 || thT < thLo || thT > thHi) continue;
      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = const_cast<AtTrack *>(&t); }
      if (!best || nb < 40) continue;

      auto *gh = new TGraph();
      for (const auto &h : best->GetHitArray()) {
         auto q = h->GetPosition();
         gh->SetPoint(gh->GetN(), q.X(), q.Y());
      }
      double pTrue = std::sqrt(keT * keT + 2 * keT * mp) / 1000.0;
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField);

      // the fit, if there is one
      double keF = -1, thF = -1, c2n = -1;
      auto *gf = new TGraph();
      if (tf) {
         tf->GetEntry(i);
         if (te && te->GetEntriesFast() > 0) {
            auto *fe = (AtTrackingEvent *)te->At(0);
            if (fe)
               for (auto &ft : fe->GetFittedTracks()) {
                  if (!ft) continue;
                  const auto &md = ft->GetTrackMetadata();
                  double ndf = md ? md->GetNdf() : 0, ch = md ? md->GetChi2() : 0;
                  c2n = ndf > 0 ? ch / ndf : -1;
                  keF = ft->GetKinematicsXtr().kineticEnergy;
                  thF = ft->GetKinematicsXtr().theta * TMath::RadToDeg();
                  for (const auto &p : ft->GetSmoothedPositions()) gf->SetPoint(gf->GetN(), p.X(), p.Y());
                  break;
               }
         }
      }
      printf("  %8lld %8.1f %8.2f | %8.1f %8.2f %10.1f | %6zu %8.1f\n", i, thT, keT, thF, keF, c2n, nb, rTrue);

      cv->cd(++shown);
      gPad->SetLeftMargin(0.14);
      gh->SetMarkerStyle(6);
      gh->SetMarkerColor(kAzure + 2);
      gh->SetTitle(TString::Format("%s  evt %lld: #theta=%.0f#circ, KE=%.1f MeV, %zu hits, #chi^{2}/ndf=%.0f;x [mm];y [mm]",
                                   title.Data(), i, thT, keT, nb, c2n));
      gh->Draw("ap");
      double cx = 0, cy = 0;
      for (int k = 0; k < gh->GetN(); ++k) { double x, y; gh->GetPoint(k, x, y); cx += x; cy += y; }
      cx /= std::max(1, gh->GetN()); cy /= std::max(1, gh->GetN());
      auto *cir = new TEllipse(cx, cy, rTrue, rTrue);
      cir->SetFillStyle(0); cir->SetLineColor(kRed + 1); cir->SetLineWidth(2); cir->Draw();
      if (gf->GetN() > 1) {
         gf->SetLineColor(kGreen + 2); gf->SetLineWidth(2); gf->Draw("l same");
      }
      auto *tx = new TLatex(0.16, 0.85, "red: circle of the TRUE radius");
      tx->SetNDC(); tx->SetTextSize(0.036); tx->SetTextColor(kRed + 1); tx->Draw();
      if (gf->GetN() > 1) {
         auto *t2 = new TLatex(0.16, 0.80, "green: the FITTED trajectory");
         t2->SetNDC(); t2->SetTextSize(0.036); t2->SetTextColor(kGreen + 2); t2->Draw();
      }
   }
   TString out = outDir + "dp_view7T_" + TString(stem(stem.Last('/') + 1, 99)) + ".png";
   cv->SaveAs(out);
   printf("\n  wrote %s\n\n", out.Data());
   fs->Close(); fr->Close(); if (ff) ff->Close();
}
