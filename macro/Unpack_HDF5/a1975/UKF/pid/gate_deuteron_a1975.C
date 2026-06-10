/// @file gate_deuteron_a1975.C
/// @brief Draw a DEUTERON gate on the FAITHFUL Spyral PID plane (the upper/orange
///        band) and show the gated-track kinematics — the 16C(p,d) candidate
///        selection on the trustworthy first-arc+spline estimator.
///
/// Uses Bρ axis to 2.0 so the deuteron band ABOVE the proton band is actually
/// visible (the earlier pid_spyral/pid_stats plots capped Bρ ~0.5–1.0 and sliced
/// it off-scale). Reads the precomputed AtPIDEvent (Spyral results) from
/// <run>_pid.root, so it is fast (no refit). Gate = pid/deuteron_band.json
/// (sqrtdEdx vs brho), parallel to pid/proton_band.json.
///
///   root -b -q 'pid/gate_deuteron_a1975.C("run_0106,...,run_0115","/mnt/f/a1975/reco/")'

void gate_deuteron_a1975(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,"
                                            "run_0114,run_0115",
                         TString inDir = "/mnt/f/a1975/reco/", TString gateFile = "pid/deuteron_band.json")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TH2F *hS = new TH2F("hS", "Faithful Spyral PID + deuteron gate;#sqrt{dEdx};B#rho [T m]", 300, 0, 40, 400, 0, 2.0);
   TH2F *hkin =
      new TH2F("hkin", "Deuteron-gated kinematics;#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200, 0, 2.0);
   TH1F *hbr = new TH1F("hbr", "Deuteron-gated B#rho;B#rho [T m];tracks", 200, 0, 2.0);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nAll = 0, nIn = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = inDir + run + "_pid.root";
      if (gSystem->AccessPathName(fn)) {
         printf("skip (missing) %s\n", fn.Data());
         continue;
      }
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPIDEvent", &pe);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &r : ev->GetSpyral())
            if (r.valid) {
               ++nAll;
               hS->Fill(r.sqrtdEdx, r.brho);
               if (pid.IsInside(r.sqrtdEdx, r.brho)) {
                  ++nIn;
                  hkin->Fill(r.polar * TMath::RadToDeg(), r.brho);
                  hbr->Fill(r.brho);
               }
            }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nspyral-valid=%ld   in deuteron gate=%ld (%.2f%%)\n", nAll, nIn, 100.0 * nIn / std::max(1L, nAll));

   // gate polygon overlay (closed)
   auto &vts = pid.GetCut().GetVertices();
   TPolyLine *pl = new TPolyLine(vts.size() + 1);
   for (size_t k = 0; k < vts.size(); ++k)
      pl->SetPoint(k, vts[k].first, vts[k].second);
   pl->SetPoint(vts.size(), vts[0].first, vts[0].second);
   pl->SetLineColor(kRed);
   pl->SetLineWidth(3);

   TCanvas *c = new TCanvas("c", "deuteron_gate", 1550, 520);
   c->Divide(3, 1);
   c->cd(1);
   gPad->SetLogz();
   hS->Draw("colz");
   pl->Draw("same");
   c->cd(2);
   gPad->SetLogz();
   hkin->Draw("colz");
   c->cd(3);
   hbr->Draw();
   c->SaveAs("pid/plots/deuteron_gate.png");
   printf("saved pid/plots/deuteron_gate.png\n");
}
