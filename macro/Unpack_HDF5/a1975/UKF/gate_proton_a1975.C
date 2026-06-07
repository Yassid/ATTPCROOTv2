/// @file gate_proton_a1975.C
/// @brief Build / apply a proton-band PID gate on the DEFAULT (Spyral) PID plane.
///
/// Reads the AtPIDEvent branch (GetDefault() = Spyral first-arc+spline PID) from the
/// <run>_pid.root files, draws the dEdx-vs-brho plane, overlays a proton-band polygon
/// (AtTools::AtCut2D), reports how many tracks fall inside, and (write=true) saves
/// the gate as a Spyral-format AtParticleID JSON (proton, Z=1 A=1) for reuse.
///
///   root -b -q 'gate_proton_a1975.C("run_0106,...","/mnt/f/a1975/reco/", false)'  // draw plane only
///   root -b -q 'gate_proton_a1975.C(runs, dir, true)'                              // + apply + write JSON

#include <vector>

void gate_proton_a1975(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,"
                                          "run_0114,run_0115",
                       TString inDir = "/mnt/f/a1975/reco/", Bool_t applyGate = false)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   // ---- proton-band polygon, traced from the Spyral plane (sqrt(dEdx), brho Tm) ----
   // sqrt(dEdx) on x (Spyral native PID axis). Upper edge along the band top from
   // high-brho/low-sqrtdEdx down to low-brho/high-sqrtdEdx, then lower edge back.
   std::vector<std::pair<double, double>> poly = {
      {4.5, 1.05}, {6.0, 0.70}, {8.0, 0.45}, {10.5, 0.32}, {14.0, 0.235}, {18.0, 0.185},
      {24.0, 0.150}, {31.0, 0.130}, {39.0, 0.118}, // upper edge
      {39.0, 0.090}, {31.0, 0.098}, {24.0, 0.110}, {18.0, 0.130}, {14.0, 0.160},
      {10.5, 0.205}, {8.0, 0.285}, {6.0, 0.42}, {4.5, 0.62}}; // lower edge back up
   AtTools::AtCut2D cut("proton_band", poly, "sqrtdEdx", "brho");

   TH2F *h = new TH2F("h", "DEFAULT PID (Spyral) + proton gate;#sqrt{dEdx};B#rho [T m]", 400, 0, 40, 400, 0, 1.4);
   TH2F *hkin = new TH2F("hkin", "gated proton kinematics;#theta [deg];B#rho [T m]", 180, 0, 180, 250, 0, 0.45);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nAll = 0, nIn = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = inDir + run + "_pid.root";
      if (gSystem->AccessPathName(fn))
         continue;
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
         for (auto &r : ev->GetDefault()) { // GetDefault() == Spyral (adopted default)
            if (!r.valid)
               continue;
            h->Fill(r.sqrtdEdx, r.brho);
            ++nAll;
            if (applyGate && cut.IsInside(r.sqrtdEdx, r.brho)) {
               ++nIn;
               hkin->Fill(r.polar * TMath::RadToDeg(), r.brho);
            }
         }
      }
      f->Close();
   }
   printf("\nDefault(Spyral) valid tracks: %ld\n", nAll);
   if (applyGate)
      printf("inside proton gate: %ld  (%.1f%%)\n", nIn, nAll ? 100.0 * nIn / nAll : 0);

   // ---- overlay the polygon ----
   TGraph *g = new TGraph;
   for (size_t k = 0; k < poly.size(); ++k)
      g->SetPoint(k, poly[k].first, poly[k].second);
   g->SetPoint(poly.size(), poly[0].first, poly[0].second); // close
   g->SetLineColor(kRed);
   g->SetLineWidth(3);

   TCanvas *c = new TCanvas("c", "protongate", applyGate ? 1500 : 850, 640);
   if (applyGate)
      c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz();
   h->Draw("colz");
   g->Draw("L same");
   if (applyGate) {
      c->cd(2);
      hkin->Draw("colz");
   }
   c->SaveAs("proton_gate.png");
   printf("saved proton_gate.png\n");

   if (applyGate) {
      AtTools::AtParticleID pid(cut, 1, 1); // proton: Z=1, A=1
      pid.WriteJSON("proton_band.json");
      printf("wrote proton_band.json (Spyral-format AtParticleID gate)\n");
   }
}
