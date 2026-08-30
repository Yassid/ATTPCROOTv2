/// @file dp_viewproj_C14.C
/// @brief 7 T tracks in XY, XZ and YZ, chosen by PRA GeoTheta band. Truth-free.
///
///   root -b -q 'dp_viewproj_C14.C("/mnt/f/a1954_C14dp_hits7T/attpc","7 T AT-TPC")'
///
/// The XZ/YZ projections are the ones that matter for the 7 T failure. genfit is handed the
/// clusters sorted by z_lab, which only orders the track while the helix advances in z faster
/// than the z resolution. At 6 MHz sampling one time bucket is ~1.7 mm of z, and near
/// theta_lab = 90 deg the measured z-step per cluster falls to 0.22 mm -- so a whole arc of the
/// spiral collapses into one or two time buckets and the sort returns noise.
///
/// One row per event: XY (the spiral seen down the beam), then XZ and YZ (where the pitch lives).
/// Clusters in blue, genfit's smoothed trajectory in green. A diverging fit leaves the blue.
/// Each row is labelled with GeoTheta, the fitted angle and chi2/ndf, so it reads on its own.

#include <algorithm>
#include <vector>

void dp_viewproj_C14(TString stem = "/mnt/f/a1954_C14dp_hits7T/attpc", TString title = "7 T AT-TPC",
                     Double_t gtLo = 88, Double_t gtHi = 92, Int_t nShow = 3, TString tag = "band88_92",
                     TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);
   const double ZPAD = 1000.0, R2D = 57.29577951;

   TFile *fr = TFile::Open(stem + "_reco.root");
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TFile *ff = TFile::Open(stem + "_genfit.root");
   TTree *tf = (TTree *)ff->Get("cbmsim");
   TClonesArray *pe = nullptr, *te = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   auto *c = new TCanvas("cproj", "proj", 1650, 420 * nShow);
   c->Divide(3, nShow);
   Int_t shown = 0;
   Long64_t N = std::min(tr->GetEntries(), tf->GetEntries());

   for (Long64_t i = 0; i < N && shown < nShow; ++i) {
      tr->GetEntry(i);
      tf->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0 || !te || te->GetEntriesFast() == 0)
         continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      auto *fv = (AtTrackingEvent *)te->At(0);
      if (!ev || !fv)
         continue;
      AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) {
            nb = t.GetHitArray().size();
            best = const_cast<AtTrack *>(&t);
         }
      if (!best || nb < 200)
         continue;
      auto *cl = best->GetHitClusterArray();
      if (!cl || cl->size() < 20)
         continue;
      double gt = best->GetGeoTheta() * R2D;
      if (gt < gtLo || gt >= gtHi)
         continue;

      std::vector<double> X, Y, Z;
      for (auto &cc : *cl) {
         auto p = cc.GetPosition();
         X.push_back(p.X());
         Y.push_back(p.Y());
         Z.push_back(ZPAD - p.Z());
      }
      // the fitted trajectory, and chi2
      std::vector<double> FX, FY, FZ;
      double chi = -1, thF = -1;
      for (auto &ftk : fv->GetFittedTracks()) {
         if (!ftk)
            continue;
         const auto &md = ftk->GetTrackMetadata();
         if (md && md->GetNdf() > 0)
            chi = md->GetChi2() / md->GetNdf();
         thF = ftk->GetKinematics().theta * R2D;
         for (auto &sp : ftk->GetSmoothedPositions()) {
            FX.push_back(sp.X());
            FY.push_back(sp.Y());
            FZ.push_back(sp.Z());
         }
         break;
      }

      const char *an[3] = {"XY", "XZ", "YZ"};
      for (int a = 0; a < 3; ++a) {
         c->cd(shown * 3 + a + 1);
         gPad->SetGrid();
         // pick the two coordinates for this projection
         std::vector<double> *hu, *hv, *fu, *fv2;
         if (a == 0) { hu = &X; hv = &Y; fu = &FX; fv2 = &FY; }
         else if (a == 1) { hu = &Z; hv = &X; fu = &FZ; fv2 = &FX; }
         else { hu = &Z; hv = &Y; fu = &FZ; fv2 = &FY; }
         // frame on the DATA, padded, with the fit included so divergence is visible
         double u0 = *std::min_element(hu->begin(), hu->end()), u1 = *std::max_element(hu->begin(), hu->end());
         double v0 = *std::min_element(hv->begin(), hv->end()), v1 = *std::max_element(hv->begin(), hv->end());
         if (!fu->empty()) {
            u0 = std::min(u0, *std::min_element(fu->begin(), fu->end()));
            u1 = std::max(u1, *std::max_element(fu->begin(), fu->end()));
            v0 = std::min(v0, *std::min_element(fv2->begin(), fv2->end()));
            v1 = std::max(v1, *std::max_element(fv2->begin(), fv2->end()));
         }
         double du = std::max(5.0, 0.10 * (u1 - u0)), dv = std::max(5.0, 0.10 * (v1 - v0));
         TString ttl;
         if (a == 0)
            ttl = TString::Format("evt %lld  GeoTheta %.1f  fit %.1f  #chi^{2}/ndf %.0f;x [mm];y [mm]",
                                  i, gt, thF, chi);
         else
            ttl = TString::Format("%s;z_{lab} [mm];%s [mm]", an[a], a == 1 ? "x" : "y");
         auto *fm = gPad->DrawFrame(u0 - du, v0 - dv, u1 + du, v1 + dv, ttl);
         fm->GetXaxis()->SetTitleSize(0.045);
         fm->GetYaxis()->SetTitleSize(0.045);
         auto *g = new TGraph(hu->size(), hu->data(), hv->data());
         g->SetMarkerStyle(20);
         g->SetMarkerSize(0.45);
         g->SetMarkerColor(kAzure + 2);
         g->Draw("P same");
         if (!fu->empty()) {
            auto *gf = new TGraph(fu->size(), fu->data(), fv2->data());
            gf->SetLineColor(kGreen + 2);
            gf->SetLineWidth(2);
            gf->Draw("L same");
         }
         if (a == 0) {
            auto *lg = new TLegend(0.60, 0.78, 0.98, 0.96);
            lg->SetBorderSize(0);
            lg->SetFillStyle(0);
            lg->AddEntry(g, "clusters", "p");
            if (!fu->empty()) lg->AddEntry((TObject *)nullptr, "green = genfit", "");
            lg->Draw();
         }
      }
      ++shown;
   }
   TString out = outDir + "dp_proj_" + tag + ".png";
   c->SaveAs(out);
   printf("\n  %s : %d events drawn, GeoTheta %.0f-%.0f\n  wrote %s\n", title.Data(), shown, gtLo, gtHi, out.Data());
}
