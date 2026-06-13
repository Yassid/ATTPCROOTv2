/// @file contactSheet.C
/// @brief One-glance contact sheet of fitted events: a grid of XY-projection thumbnails,
///        one pad per (gated) event, each labeled ev# / theta / KE / chi2ndf. Lets you
///        scan the whole set visually, then open the full 3-projection PNG (inspectFitEvent.C)
///        for any event of interest. Static PNG (WSL has no working Eve/OpenGL).
///
///   root -b -q 'pipeline/contactSheet.C("/tmp/run_0106_genfitter_gated.root","pid/deuteron_band.json","/tmp/contact.png",true,6)'
///   gatedOnly=false -> every event with a fit; ncol = thumbnails per row.

void contactSheet(TString file = "/tmp/run_0106_genfitter_gated.root", TString gateFile = "pid/deuteron_band.json",
                  TString outPng = "/tmp/contact_sheet.png", Bool_t gatedOnly = true, Int_t ncol = 6,
                  Double_t zPadPlane = 1000.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TFile f(file);
   TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *trkEvt = nullptr, *pidEvt = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &trkEvt);
   t->SetBranchAddress("AtPIDEvent", &pidEvt);

   // --- pass 1: collect the events to show ---
   struct Item { Int_t iev; double th, ke, c2; int conv; };
   std::vector<Item> items;
   Long64_t N = t->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      auto *te = (AtTrackingEvent *)trkEvt->At(0);
      auto *pe = (pidEvt && pidEvt->GetEntriesFast() > 0) ? (AtPIDEvent *)pidEvt->At(0) : nullptr;
      if (!te) continue;
      std::map<int, bool> inGate;
      if (pe)
         for (auto &r : pe->GetSpyral())
            inGate[r.trackID] = r.valid && pid.IsInside(r.sqrtdEdx, r.brho);
      auto &fits = te->GetFittedTracks();
      // representative fitted track for the label: prefer an in-gate one
      const AtFittedTrack *rep = nullptr; bool anyGate = false;
      for (auto &fp : fits) {
         if (!fp) continue;
         bool g = inGate.count(fp->GetTrackID()) ? inGate[fp->GetTrackID()] : false;
         if (g) anyGate = true;
         if (!rep || g) rep = fp.get();
      }
      if (gatedOnly && !anyGate) continue;
      if (fits.empty()) continue;
      double th = rep ? rep->GetKinematics(0).theta * TMath::RadToDeg() : -1;
      double ke = rep ? rep->GetKinematics(0).kineticEnergy : -1;
      double c2 = -1; int cv = -1;
      if (rep) { auto &m = const_cast<AtFittedTrack *>(rep)->GetTrackMetadata();
                 if (m) { c2 = m->GetNdf() > 0 ? m->GetChi2() / m->GetNdf() : -1; cv = m->GetFitConverged(); } }
      items.push_back({(Int_t)i, th, ke, c2, cv});
   }
   printf("contact sheet: %zu events\n", items.size());
   if (items.empty()) return;

   Int_t nrow = (items.size() + ncol - 1) / ncol;
   TCanvas *c = new TCanvas("cs", "contact sheet", ncol * 300, nrow * 300);
   c->Divide(ncol, nrow, 0.001, 0.001);

   auto pickXY = [&](double X, double Y, double Z, double &a, double &b) { a = X; b = Y; };
   const Color_t pal[] = {kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1, kCyan + 2, kMagenta + 1, kAzure + 7};

   // --- pass 2: draw each thumbnail ---
   for (size_t k = 0; k < items.size(); ++k) {
      Int_t iev = items[k].iev;
      t->GetEntry(iev);
      auto *te = (AtTrackingEvent *)trkEvt->At(0);
      auto *pe = (pidEvt && pidEvt->GetEntriesFast() > 0) ? (AtPIDEvent *)pidEvt->At(0) : nullptr;
      std::map<int, bool> inGate;
      if (pe)
         for (auto &r : pe->GetSpyral())
            inGate[r.trackID] = r.valid && pid.IsInside(r.sqrtdEdx, r.brho);
      auto tracks = te->GetTrackArray();
      auto &fits = te->GetFittedTracks();

      c->cd(k + 1);
      gPad->SetGrid();
      gPad->DrawFrame(-150, -150, 150, 150, Form("ev %d;X;Y", iev));
      for (size_t i = 0; i < tracks.size(); ++i) {
         auto &tr = tracks[i];
         bool g = inGate.count(tr.GetTrackID()) ? inGate[tr.GetTrackID()] : false;
         auto &hits = tr.GetHitArray();
         if (hits.empty()) continue;
         auto *gr = new TGraph();
         for (auto &h : hits) { const auto &p = h->GetPosition(); double a, b; pickXY(p.X(), p.Y(), zPadPlane - p.Z(), a, b); gr->SetPoint(gr->GetN(), a, b); }
         gr->SetMarkerStyle(g ? 20 : 24);
         gr->SetMarkerSize(g ? 0.45 : 0.3);
         gr->SetMarkerColor(g ? pal[i % 7] : (kGray + 1));
         gr->Draw("P SAME");
      }
      for (auto &fp : fits) {
         if (!fp) continue;
         bool g = inGate.count(fp->GetTrackID()) ? inGate[fp->GetTrackID()] : false;
         if (!g) continue;
         const auto &kin = fp->GetKinematics(0); const auto &v = fp->GetVertex(0);
         double a, b; pickXY(v.X(), v.Y(), v.Z(), a, b);
         auto *m = new TMarker(a, b, 29); m->SetMarkerSize(1.6); m->SetMarkerColor(kRed + 1); m->Draw();
         if (kin.kineticEnergy > 0 && std::isfinite(kin.theta) && std::isfinite(kin.phi)) {
            double L = 55, ux = std::sin(kin.theta) * std::cos(kin.phi), uy = std::sin(kin.theta) * std::sin(kin.phi);
            double a1, b1; pickXY(v.X() + L * ux, v.Y() + L * uy, 0, a1, b1);
            auto *l = new TLine(a, b, a1, b1); l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->Draw();
         }
      }
      auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.07);
      Color_t qcol = (items[k].conv == 1 && items[k].ke > 0.5 && items[k].ke < 60 && items[k].c2 < 5) ? kGreen + 3 : kRed + 1;
      tx->SetTextColor(qcol);
      tx->DrawLatex(0.12, 0.82, Form("#theta=%.0f KE=%.1f", items[k].th, items[k].ke));
      tx->DrawLatex(0.12, 0.74, Form("#chi^{2}/n=%.2f c=%d", items[k].c2, items[k].conv));
   }
   c->SaveAs(outPng);
   printf("wrote %s  (%dx%d grid)\n", outPng.Data(), ncol, nrow);
}
