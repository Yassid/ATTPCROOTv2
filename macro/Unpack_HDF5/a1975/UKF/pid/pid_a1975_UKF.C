/// @file pid_a1975_UKF.C
/// @brief Spyral-style particle ID for a1975 16C+p, ported to C++ (AtTools::AtPID).
///
/// Reads <run>_reco.root (AtPatternEvent from PRA), computes Spyral PID
/// observables per track with AtTools::AtPIDEstimator (brho, sqrt_dEdx, polar),
/// fills the PID plane (sqrt_dEdx vs brho) and the kinematics plane (polar vs
/// brho), and saves a PNG. If a gate JSON is supplied it is loaded as an
/// AtTools::AtParticleID and applied (overlaid + gated counts), exactly mirroring
/// Spyral's particle_id workflow — but entirely in C++.
///
/// Run:
///   root -b -q 'pid_a1975_UKF.C("run_0116")'                      // just the PID plot
///   root -b -q 'pid_a1975_UKF.C("run_0116", "proton_gate.json")'  // apply a gate

// AtTools::AtCut2D / AtParticleID / AtPIDEstimator are autoloaded from libAtTools
// (no #include needed — ROOT_INCLUDE_PATH is empty in this setup; the dictionary
// PCM provides the class definitions).

void pid_a1975_UKF(TString fileName = "run_0116", TString gateFile = "", Double_t bField = 2.85,
                   Double_t smallPadRadius = 152.0)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kRainBow);

   TString inputFile = fileName + "_reco.root";
   if (gSystem->AccessPathName(inputFile)) {
      printf("\033[1;31mERROR: %s not found (run unpackReco_a1975_UKF.C first)\033[0m\n", inputFile.Data());
      return;
   }

   TFile *f = TFile::Open(inputFile);
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *peArr = nullptr;
   // Only read AtPatternEvent — the file also holds the large AtRawEvent branch.
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &peArr);

   AtTools::AtPIDEstimator estimator(bField, smallPadRadius);

   // Optional gate
   AtTools::AtParticleID pid;
   bool haveGate = false;
   if (gateFile.Length() > 0) {
      pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
      haveGate = pid.IsValid();
      if (haveGate)
         printf("Applying gate '%s' (Z=%d A=%d, mass %.2f MeV) on axes %s vs %s\n", pid.GetName().c_str(), pid.GetZ(),
                pid.GetA(), pid.GetMassMeV(), pid.GetXAxis().c_str(), pid.GetYAxis().c_str());
   }

   // First pass: collect observables to set sensible axis ranges.
   std::vector<double> vbrho, vsqdedx, vpolar;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (peArr->GetEntries() == 0)
         continue;
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe)
         continue;
      for (auto &track : pe->GetTrackCand()) {
         auto r = estimator.Estimate(const_cast<AtTrack &>(track));
         if (!r.valid)
            continue;
         vbrho.push_back(r.brho);
         vsqdedx.push_back(r.sqrtdEdx);
         vpolar.push_back(r.polar * TMath::RadToDeg());
      }
   }
   if (vbrho.empty()) {
      printf("\033[1;31mNo valid PID observables.\033[0m\n");
      return;
   }
   auto pct = [](std::vector<double> v, double p) {
      std::sort(v.begin(), v.end());
      return v[std::min(v.size() - 1, (size_t)(p * v.size()))];
   };
   // Physical display caps: beam-like tracks (sin theta -> 0) give brho -> inf,
   // a long tail that hides the real bands. Cap like Spyral's kinematics plot.
   double sqMax = std::min(pct(vsqdedx, 0.98) * 1.3, 20.0);
   double brMax = 3.0;
   printf("N=%zu tracks | brho: med=%.3f max(98%%)=%.3f Tm | sqrt_dEdx: med=%.1f max(98%%)=%.1f\n", vbrho.size(),
          pct(vbrho, 0.5), pct(vbrho, 0.98), pct(vsqdedx, 0.5), pct(vsqdedx, 0.98));

   TH2F *hpid = new TH2F("hpid", "Particle ID (C++ Spyral port);sqrt(dEdx) [#sqrt{ADC/mm}];B#rho [T m]", 200, 0, sqMax,
                         200, 0, brMax);
   TH2F *hkin = new TH2F("hkin", "Kinematics;#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200, 0, brMax);
   TH2F *hgated = haveGate ? new TH2F("hgated", "PID gated;sqrt(dEdx);B#rho [T m]", 200, 0, sqMax, 200, 0, brMax)
                          : nullptr;
   int nGated = 0, nTot = (int)vbrho.size();

   for (size_t k = 0; k < vbrho.size(); ++k) {
      hpid->Fill(vsqdedx[k], vbrho[k]);
      hkin->Fill(vpolar[k], vbrho[k]);
      if (haveGate && pid.IsInside(vsqdedx[k], vbrho[k])) {
         hgated->Fill(vsqdedx[k], vbrho[k]);
         ++nGated;
      }
   }

   TCanvas *c = new TCanvas("c", "pid", haveGate ? 1500 : 1100, 500);
   c->Divide(haveGate ? 3 : 2, 1);
   c->cd(1);
   gPad->SetLogz();
   hpid->Draw("colz");
   // overlay the gate polygon if present
   if (haveGate) {
      const auto &v = pid.GetCut().GetVertices();
      TPolyLine *pl = new TPolyLine(v.size() + 1);
      for (size_t k = 0; k < v.size(); ++k)
         pl->SetPoint(k, v[k].first, v[k].second);
      pl->SetPoint(v.size(), v[0].first, v[0].second);
      pl->SetLineColor(kRed);
      pl->SetLineWidth(2);
      pl->Draw();
   }
   c->cd(2);
   gPad->SetLogz();
   hkin->Draw("colz");
   if (haveGate) {
      c->cd(3);
      gPad->SetLogz();
      hgated->Draw("colz");
      printf("\033[1;32mGated %d / %d tracks inside '%s'\033[0m\n", nGated, nTot, pid.GetName().c_str());
   }
   TString png = "pid/plots/" + fileName + "_pid.png";
   c->SaveAs(png);
   printf("Saved %s\n", png.Data());
}
