/// @file pid_overlay.C
/// @brief Mass-PID overlay: pions (branch 8) and kaons (branch 10) on ONE
///        Brho vs dE/dx panel. Each species has its own digi file (pure sample,
///        so no truth match needed). Bulk momenta are matched (both p~375 MeV/c)
///        so the dE/dx separation is purely the mass: a K+ at the same momentum
///        is slower (beta = p/E) -> more ionizing -> a higher dE/dx band.
///        dE/dx = Sum q / phi-sorted 3D arc length (fitter-independent);
///        signed Brho from the fit KE + the species mass.
/// Run: root -b -q pid_overlay.C
void pid_overlay(TString digiPi = "./data/output_digi_pi.root", Double_t mPi = 139.57039,
                 TString digiK = "./data/output_digi_K.root", Double_t mK = 493.677,
                 TString fitter = "AtTrackingEventUKF", TString outPng = "./data/pid_overlay.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   auto fill = [&](TString digiFile, double mass, TGraph *g) {
      TFile f(digiFile);
      TTree *t = (TTree *)f.Get("cbmsim");
      if (!t) { printf("missing %s\n", digiFile.Data()); return; }
      TClonesArray *te = new TClonesArray("AtTrackingEvent"); t->SetBranchAddress(fitter, &te);
      TClonesArray *pe = new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent", &pe);
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         if (!te->GetEntries() || !pe->GetEntries()) continue;
         auto &tracks = ((AtPatternEvent *)pe->At(0))->GetTrackCand();
         for (const auto &ft : ((AtTrackingEvent *)te->At(0))->GetFittedTracks()) {
            double KE = ft->GetKinematics(0).kineticEnergy;
            if (!(KE > 0)) continue;
            double p = std::sqrt(KE * KE + 2 * KE * mass);
            const auto &pinfo = ft->GetParticleInfo(0);
            int sgn = pinfo.charge != 0 ? (pinfo.charge > 0 ? 1 : -1) : 0;
            if (!sgn) { TString id = pinfo.idPDG; sgn = id.Contains("+") ? 1 : (id.Contains("-") ? -1 : 0); }
            if (!sgn) continue;
            double brho = sgn * (p / 1000.0) / 0.299792458;
            int tid = ft->GetTrackID();
            if (tid < 0 || tid >= (int)tracks.size()) continue;
            const auto &hits = tracks[tid].GetHitArray();
            if (hits.size() < 3) continue;
            auto cen = tracks[tid].GetGeoCenter();
            std::vector<std::tuple<double, double, double, double, double>> s;
            for (const auto &h : hits) { const auto &q = h->GetPosition();
               s.emplace_back(std::atan2(q.Y() - cen.second, q.X() - cen.first), q.X(), q.Y(), q.Z(), h->GetCharge()); }
            std::sort(s.begin(), s.end());
            double totQ = 0, len = 0;
            for (size_t i = 0; i < s.size(); ++i) { totQ += std::get<4>(s[i]); if (i) {
               double dx = std::get<1>(s[i]) - std::get<1>(s[i - 1]), dy = std::get<2>(s[i]) - std::get<2>(s[i - 1]), dz = std::get<3>(s[i]) - std::get<3>(s[i - 1]);
               len += std::sqrt(dx * dx + dy * dy + dz * dz); } }
            if (len > 0) g->SetPoint(g->GetN(), brho, totQ / len);
         }
      }
   };

   auto *gPi = new TGraph(); auto *gK = new TGraph();
   fill(digiPi, mPi, gPi);
   fill(digiK, mK, gK);

   // dE/dx axis capped at the 98th percentile over both samples
   std::vector<double> dv;
   for (auto *g : {gPi, gK}) for (int i = 0; i < g->GetN(); ++i) dv.push_back(g->GetY()[i]);
   std::sort(dv.begin(), dv.end());
   double dmax = dv.empty() ? 1 : dv[(size_t)(0.98 * dv.size())] * 1.15;
   double bmax = 3.0;

   auto *c = new TCanvas("pidov", "mass-PID overlay", 900, 700);
   c->SetGrid();
   c->DrawFrame(-bmax, 0, bmax, dmax,
                "PUMA mass PID: #pi vs K  (signed B#rho vs dE/dx);signed B#rho [T#upoint m];dE/dx [a.u./mm]");
   gK->SetMarkerColor(kRed + 1);   gK->SetMarkerStyle(20); gK->SetMarkerSize(0.5); gK->Draw("P");
   gPi->SetMarkerColor(kAzure + 2); gPi->SetMarkerStyle(20); gPi->SetMarkerSize(0.5); gPi->Draw("P");
   auto *leg = new TLegend(0.63, 0.78, 0.98, 0.93); leg->SetTextSize(0.032);
   leg->AddEntry(gPi, Form("#pi^{+}/#pi^{-}  (%d)", gPi->GetN()), "p");
   leg->AddEntry(gK, Form("K^{+}  (%d)", gK->GetN()), "p");
   leg->Draw();
   c->SaveAs(outPng);
   printf("mass-PID overlay (%s) -> %s   [pi %d, K %d]\n", fitter.Data(), outPng.Data(), gPi->GetN(), gK->GetN());
}
