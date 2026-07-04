/// @file pid_test8.C
/// @brief PID / charge-separation plot for the PUMA branch-8 sample. For each
///        FITTED track: signed Brho (from the fit momentum + fit charge) vs dE/dx
///        (Sum q / phi-sorted 3D arc length of the source AtTrack). Two panels
///        (UKF | genfit); points coloured by TRUTH charge (pi+ green, pi- magenta),
///        so a mis-charged fit shows as a wrong-colour point in the wrong Brho band.
///        Truth |p|=374.9 MeV/c => |Brho| = 1.25 T*m (dashed guide lines).
/// Run: root -b -q pid_test8.C
/// colorMode: "confusion" (default) colours each point by whether the FIT charge
/// matches truth (correct = grey, WRONG = red, so the mis-ID population pops).
/// "truth" colours by truth species (pi+ green, pi- magenta).
/// species: "pi" (branch 8) or "K"/"kaon" (branch 10) — sets the mass used for p
/// (from KE) and the truth |Brho| guide line. testEnergy overrides per-particle E.
void pid_test8(TString colorMode = "confusion", TString species = "pi", Double_t testEnergy = -1,
               TString digiFile = "./data/output_digi_both8.root",
               TString simFile = "./data/attpcsim.root", TString outPng = "./data/pid_test8.png")
{
   const bool confusion = (colorMode != "truth");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const bool isK = (species == "K" || species == "kaon");
   const double m_pi = isK ? 493.677 : 139.57039;
   const double E0 = (testEnergy > 0) ? testEnergy : (isK ? 0.777 : 0.4);
   const double p0 = std::sqrt(E0 * E0 - (m_pi / 1000) * (m_pi / 1000)) * 1000; // MeV/c
   const double brho0 = (p0 / 1000.0) / 0.299792458;                            // T*m
   const double kTol = 3.0;

   TFile fD(digiFile); TTree *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);  TTree *tS = (TTree *)fS.Get("cbmsim");
   TClonesArray *ukfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *gfArr = new TClonesArray("AtTrackingEvent");
   TClonesArray *patArr = new TClonesArray("AtPatternEvent");
   tD->SetBranchAddress("AtTrackingEventUKF", &ukfArr);
   tD->SetBranchAddress("AtTrackingEventGenfit", &gfArr);
   tD->SetBranchAddress("AtPatternEvent", &patArr);
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tS->SetBranchAddress("AtTpcPoint", &mcPts);
   tS->SetBranchAddress("MCTrack", &mcTrks);

   // one TGraph of (brho, dedx) per (fitter, category). confusion: 0=correct,1=WRONG.
   // truth: 0=pi+, 1=pi-.
   TGraph *g[2][2];
   for (int f = 0; f < 2; ++f) for (int s = 0; s < 2; ++s) g[f][s] = new TGraph();

   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());
   for (Long64_t e = 0; e < nE; ++e) {
      tD->GetEntry(e); tS->GetEntry(e);
      if (!patArr->GetEntries()) continue;
      auto *pat = (AtPatternEvent *)patArr->At(0);
      auto &tracks = pat->GetTrackCand();
      int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC); std::vector<int> mcPdg(nMC);
      for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcPts->At(k); mcX[k] = mp->GetX() * 10; mcY[k] = mp->GetY() * 10;
         int t = mp->GetTrackID(); auto *mt = (t >= 0 && t < mcTrks->GetEntries()) ? (AtMCTrack *)mcTrks->At(t) : nullptr; mcPdg[k] = mt ? mt->GetPdgCode() : 0; }

      TClonesArray *arr[2] = {ukfArr, gfArr};
      for (int f = 0; f < 2; ++f) {
         if (!arr[f]->GetEntries()) continue;
         auto *te = (AtTrackingEvent *)arr[f]->At(0);
         for (const auto &ft : te->GetFittedTracks()) {
            double KE = ft->GetKinematics(0).kineticEnergy;
            if (!(KE > 0)) continue;
            double p = std::sqrt(KE * KE + 2 * KE * m_pi);        // MeV/c
            // fit charge sign
            const auto &pinfo = ft->GetParticleInfo(0);
            int fitSign = pinfo.charge != 0 ? (pinfo.charge > 0 ? 1 : -1) : 0;
            if (!fitSign) { TString id = pinfo.idPDG; fitSign = id.Contains("+") ? 1 : (id.Contains("-") ? -1 : 0); }
            if (!fitSign) continue;
            double brho = fitSign * (p / 1000.0) / 0.299792458;   // signed T*m

            int tid = ft->GetTrackID();
            if (tid < 0 || tid >= (int)tracks.size()) continue;
            const auto &hits = tracks[tid].GetHitArray();
            if (hits.size() < 3) continue;
            // dE/dx = sum q / phi-sorted 3D arc length
            auto cen = tracks[tid].GetGeoCenter();
            std::vector<std::tuple<double, double, double, double, double>> s; // phi,x,y,z,q
            for (const auto &h : hits) { const auto &p3 = h->GetPosition();
               s.emplace_back(std::atan2(p3.Y() - cen.second, p3.X() - cen.first), p3.X(), p3.Y(), p3.Z(), h->GetCharge()); }
            std::sort(s.begin(), s.end());
            double totQ = 0, len = 0;
            for (size_t i = 0; i < s.size(); ++i) { totQ += std::get<4>(s[i]); if (i) {
               double dx = std::get<1>(s[i]) - std::get<1>(s[i - 1]), dy = std::get<2>(s[i]) - std::get<2>(s[i - 1]), dz = std::get<3>(s[i]) - std::get<3>(s[i - 1]);
               len += std::sqrt(dx * dx + dy * dy + dz * dz); } }
            if (!(len > 0)) continue;
            double dedx = totQ / len;

            // truth charge via (x,y) hit majority vote
            std::map<int, int> votes;
            for (const auto &h : hits) { const auto &p3 = h->GetPosition(); double best = kTol * kTol; int bp = 0;
               for (int k = 0; k < nMC; ++k) { double d2 = (p3.X() - mcX[k]) * (p3.X() - mcX[k]) + (p3.Y() - mcY[k]) * (p3.Y() - mcY[k]); if (d2 < best) { best = d2; bp = mcPdg[k]; } }
               if (bp) votes[bp]++; }
            int tp = 0, bv = 0; for (auto &kv : votes) if (kv.second > bv) { bv = kv.second; tp = kv.first; }
            if (!tp) continue;
            int truthSign = (tp > 0) ? 1 : -1;
            int idx = confusion ? ((fitSign == truthSign) ? 0 : 1)   // 0=correct, 1=WRONG
                                : ((tp > 0) ? 0 : 1);                 // 0=pi+, 1=pi-
            g[f][idx]->SetPoint(g[f][idx]->GetN(), brho, dedx);
         }
      }
   }

   auto *c = new TCanvas("pid", "PID branch-8", 1400, 640);
   c->Divide(2, 1);
   const char *fname[2] = {"UKF", "genfit"};
   for (int f = 0; f < 2; ++f) {
      c->cd(f + 1); gPad->SetGrid();
      // cap dE/dx axis at the 97th percentile so a few short-track outliers don't
      // compress the bulk (single-species pions -> dedx is a tight band).
      std::vector<double> dv;
      for (int s = 0; s < 2; ++s) for (int i = 0; i < g[f][s]->GetN(); ++i) dv.push_back(g[f][s]->GetY()[i]);
      std::sort(dv.begin(), dv.end());
      double dmax = dv.empty() ? 1 : dv[(size_t)(0.97 * dv.size())] * 1.15;
      double bAx = std::max(2.2, brho0 * 1.4);
      auto *fr = gPad->DrawFrame(-bAx, 0, bAx, dmax,
                                 Form("%s  PID: signed B#rho vs dE/dx (%s);signed B#rho [T#upoint m];dE/dx [a.u./mm]", fname[f],
                                      confusion ? "red = WRONG charge" : "colour = truth charge"));
      (void)fr;
      for (double bx : {-brho0, brho0}) { auto *l = new TLine(bx, 0, bx, dmax); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
      int n0 = g[f][0]->GetN(), n1 = g[f][1]->GetN();
      if (confusion) {
         // correct (grey, small) first, then WRONG (red, bigger) on top so they pop
         g[f][0]->SetMarkerColor(kGray + 2); g[f][0]->SetMarkerStyle(20); g[f][0]->SetMarkerSize(0.45); g[f][0]->Draw("P");
         g[f][1]->SetMarkerColor(kRed);      g[f][1]->SetMarkerStyle(20); g[f][1]->SetMarkerSize(0.7);  g[f][1]->Draw("P");
         auto *leg = new TLegend(0.60, 0.76, 0.98, 0.93); leg->SetTextSize(0.03);
         double acc = (n0 + n1) ? 100.0 * n0 / (n0 + n1) : 0;
         leg->AddEntry(g[f][0], Form("correct charge (%d)", n0), "p");
         leg->AddEntry(g[f][1], Form("WRONG charge (%d)", n1), "p");
         leg->AddEntry((TObject *)nullptr, Form("accuracy %.1f%%", acc), "");
         leg->Draw();
      } else {
         g[f][0]->SetMarkerColor(kGreen + 2);  g[f][0]->SetMarkerStyle(20); g[f][0]->SetMarkerSize(0.5); g[f][0]->Draw("P");
         g[f][1]->SetMarkerColor(kMagenta + 1); g[f][1]->SetMarkerStyle(20); g[f][1]->SetMarkerSize(0.5); g[f][1]->Draw("P");
         auto *leg = new TLegend(0.68, 0.78, 0.98, 0.93); leg->SetTextSize(0.03);
         leg->AddEntry(g[f][0], Form("truth #pi+ (%d)", n0), "p");
         leg->AddEntry(g[f][1], Form("truth #pi- (%d)", n1), "p");
         leg->Draw();
      }
   }
   if (outPng == "./data/pid_test8.png") // build a species/mode-specific default name
      outPng = Form("./data/pid_test%s%s.png", isK ? "10_K" : "8", confusion ? "_confusion" : "");
   c->SaveAs(outPng);
   printf("PID plot (%s) -> %s\n", colorMode.Data(), outPng.Data());
   for (int f = 0; f < 2; ++f) {
      int n0 = g[f][0]->GetN(), n1 = g[f][1]->GetN();
      if (confusion) printf("  %s: correct %d, WRONG %d  (accuracy %.1f%%)\n", fname[f], n0, n1, (n0 + n1) ? 100.0 * n0 / (n0 + n1) : 0);
      else printf("  %s: pi+ %d, pi- %d\n", fname[f], n0, n1);
   }
}
