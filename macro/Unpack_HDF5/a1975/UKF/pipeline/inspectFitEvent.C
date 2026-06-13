/// @file inspectFitEvent.C
/// @brief Event-by-event quality inspection of the AtGenfitter fit, with PID gating.
///        Mirrors the PUMA/inspect_fit.C 2D-projection style (no Eve/OpenGL): three
///        pads XY / ZX / ZY showing the per-track hit cloud, the fitted vertex and a
///        kinematic-direction stub, annotated with fit quality (KE, chi2/ndf, conv).
///
/// PID gate: each source track is matched (by trackID) to its Spyral PID result and
/// tested against the deuteron band (pid/deuteron_band.json). Tracks INSIDE the gate
/// are drawn solid + colored and their fit is fully annotated; tracks OUTSIDE are
/// dimmed grey. So we inspect specifically how the fit performs on the gated deuterons.
///
/// Hits are shown in the SAME lab frame the fitter used: z_lab = zPadPlane - z_digi.
///
///   root -b -q 'pipeline/inspectFitEvent.C("/tmp/run_0106_genfitter_val.root",0,20)'
///   root -b -q 'pipeline/inspectFitEvent.C(file, 7, 7)'   // just event 7
///   gatedOnly=true -> skip events with no in-gate track (faster scanning)

void inspectFitEvent(TString file = "/tmp/run_0106_genfitter_val.root", Int_t evtStart = 0, Int_t evtEnd = 20,
                     TString gateFile = "pid/deuteron_band.json", TString outDir = "/tmp/fitdisp/",
                     Bool_t gatedOnly = false, Double_t zPadPlane = 1000.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gSystem->mkdir(outDir, kTRUE);

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   printf("Loaded PID gate '%s' (Z=%d A=%d) from %s\n", pid.GetName().c_str(), pid.GetZ(), pid.GetA(),
          gateFile.Data());

   TFile f(file);
   TTree *t = (TTree *)f.Get("cbmsim");
   if (!t) { printf("no cbmsim tree\n"); return; }
   TClonesArray *trkEvt = nullptr, *pidEvt = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &trkEvt);
   t->SetBranchAddress("AtPIDEvent", &pidEvt);

   const Color_t palette[] = {kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1, kCyan + 2, kMagenta + 1, kAzure + 7};
   auto col = [&](int i) { return palette[i % 7]; };

   // proj 0: XY, 1: ZX (z horiz, x vert), 2: ZY (z horiz, y vert)
   auto pick = [&](int proj, double X, double Y, double Z, double &a, double &b) {
      if (proj == 0) { a = X; b = Y; }
      else if (proj == 1) { a = Z; b = X; }
      else { a = Z; b = Y; }
   };
   const char *axX[3] = {"X [mm]", "z_{lab} [mm]", "z_{lab} [mm]"};
   const char *axY[3] = {"Y [mm]", "X [mm]", "Y [mm]"};

   Int_t lo = std::max(0, evtStart), hi = std::min((Int_t)t->GetEntries() - 1, evtEnd);
   for (Int_t iev = lo; iev <= hi; ++iev) {
      t->GetEntry(iev);
      auto *te = (AtTrackingEvent *)trkEvt->At(0);
      auto *pe = (pidEvt && pidEvt->GetEntriesFast() > 0) ? (AtPIDEvent *)pidEvt->At(0) : nullptr;
      if (!te) continue;
      auto tracks = te->GetTrackArray();          // source AtTracks (with hits)
      auto &fits = te->GetFittedTracks();         // fitted results
      if (tracks.empty()) continue;

      // trackID -> spyral PID + in-gate flag
      std::map<int, AtTools::AtSpyralResult> pidById;
      std::map<int, bool> inGate;
      if (pe)
         for (auto &r : pe->GetSpyral()) {
            pidById[r.trackID] = r;
            inGate[r.trackID] = r.valid && pid.IsInside(r.sqrtdEdx, r.brho);
         }

      int nGated = 0;
      for (auto &kv : inGate) if (kv.second) nGated++;
      if (gatedOnly && nGated == 0) continue;

      // ---- per-track console table ----
      printf("\n===== event %d : %zu tracks, %zu fits, %d in deuteron gate =====\n", iev, tracks.size(), fits.size(),
             nGated);
      printf("  %-4s %-6s %-7s %-8s %-9s %-9s %-6s %-5s\n", "trk", "gate", "Rpra", "KE", "theta", "chi2/ndf", "conv",
             "nhit");
      for (auto &tr : tracks) {
         int id = tr.GetTrackID();
         bool g = inGate.count(id) ? inGate[id] : false;
         const AtFittedTrack *ft = nullptr;
         for (auto &fp : fits) if (fp && fp->GetTrackID() == id) { ft = fp.get(); break; }
         double ke = ft ? ft->GetKinematics(0).kineticEnergy : -1;
         double th = ft ? ft->GetKinematics(0).theta * TMath::RadToDeg() : -1;
         double c2 = -1; int conv = -1;
         if (ft) { auto &m = const_cast<AtFittedTrack *>(ft)->GetTrackMetadata();
                   if (m) { c2 = m->GetNdf() > 0 ? m->GetChi2() / m->GetNdf() : -1; conv = m->GetFitConverged(); } }
         printf("  %-4d %-6s %-7.1f %-8.2f %-9.2f %-9.2f %-6d %-5zu\n", id, g ? "DEUT" : "-", tr.GetGeoRadius(), ke, th,
                c2, conv, tr.GetHitArray().size());
      }

      // ---- 3-projection canvas ----
      TCanvas *c = new TCanvas(Form("c_ev%d", iev), Form("fit inspect ev %d", iev), 1800, 600);
      c->Divide(3, 1);
      for (int proj = 0; proj < 3; ++proj) {
         c->cd(proj + 1);
         gPad->SetGrid();
         double aLo, aHi, bLo, bHi;
         if (proj == 0) { aLo = bLo = -250; aHi = bHi = 250; }
         else { aLo = -50; aHi = zPadPlane + 50; bLo = -250; bHi = 250; }
         gPad->DrawFrame(aLo, bLo, aHi, bHi,
                         Form("ev %d  %s vs %s;%s;%s", iev, axY[proj], axX[proj], axX[proj], axY[proj]));

         // hits per track (lab frame), colored if gated else grey
         for (size_t i = 0; i < tracks.size(); ++i) {
            auto &tr = tracks[i];
            bool g = inGate.count(tr.GetTrackID()) ? inGate[tr.GetTrackID()] : false;
            auto &hits = tr.GetHitArray();
            if (hits.empty()) continue;
            auto *gr = new TGraph();
            for (auto &h : hits) {
               const auto &p = h->GetPosition();
               double a, b; pick(proj, p.X(), p.Y(), zPadPlane - p.Z(), a, b);
               gr->SetPoint(gr->GetN(), a, b);
            }
            gr->SetMarkerStyle(g ? 20 : 24);
            gr->SetMarkerSize(g ? 0.55 : 0.4);
            gr->SetMarkerColor(g ? col((int)i) : (kGray + 1));
            gr->Draw("P SAME");
         }

         // fitted vertex + kinematic stub for every fit (emphasize gated)
         for (auto &fp : fits) {
            if (!fp) continue;
            int id = fp->GetTrackID();
            bool g = inGate.count(id) ? inGate[id] : false;
            const auto &kin = fp->GetKinematics(0);
            const auto &v = fp->GetVertex(0);
            Color_t cc = g ? kRed + 1 : (kGray + 2);
            double a, b; pick(proj, v.X(), v.Y(), v.Z(), a, b);
            auto *m = new TMarker(a, b, 29);
            m->SetMarkerSize(g ? 2.0 : 1.2); m->SetMarkerColor(cc); m->Draw();
            if (kin.kineticEnergy > 0 && std::isfinite(kin.theta) && std::isfinite(kin.phi)) {
               const double L = 60.0;
               double ux = std::sin(kin.theta) * std::cos(kin.phi);
               double uy = std::sin(kin.theta) * std::sin(kin.phi);
               double uz = std::cos(kin.theta);
               double a1, b1; pick(proj, v.X() + L * ux, v.Y() + L * uy, v.Z() + L * uz, a1, b1);
               auto *l = new TLine(a, b, a1, b1);
               l->SetLineColor(cc); l->SetLineWidth(g ? 3 : 1); l->SetLineStyle(g ? 1 : 2); l->Draw();
            }
         }
      }
      // annotate quality on first pad
      c->cd(1);
      auto *tx = new TLatex();
      tx->SetNDC(); tx->SetTextSize(0.035); tx->SetTextColor(kRed + 1);
      tx->DrawLatex(0.13, 0.92, Form("ev %d  #deuteron-gated=%d/%zu", iev, nGated, tracks.size()));

      TString png = Form("%s/ev%05d.png", outDir.Data(), iev);
      c->SaveAs(png);
      delete c;
   }
   printf("\nPNGs written to %s\n", outDir.Data());
}
