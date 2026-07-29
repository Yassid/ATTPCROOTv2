/// @file brho_theta_C15.C
/// @brief Bρ vs θ_lab for a2091 15C from PRA track candidates (AtSpyralPID, no fit).
///        Quick-look kinematic locus straight off <run>_reco.root.
///
///   root -b -q 'pid/brho_theta_C15.C("run_0138,run_0058", "/home/yassid/a2091_C15_reco/")'
void brho_theta_C15(TString runsCSV = "run_0138", TString inDir = "/home/yassid/a2091_C15_reco/", Long64_t nEvents = -1,
                     double bField = 2.85, double arclenMin = 0, TString outTag = "", double brhoMax = 2.5)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a2091/UKF/pid/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH2F *hBt = new TH2F("hBt", "15C  B#rho vs #theta_{lab}  (PRA tracks, Spyral);#theta_{lab} [deg];B#rho [T m]", 180,
                        0, 180, 400, 0, brhoMax);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nTracks = 0, nValid = 0, nRuns = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString rf = inDir + run + "_reco.root";
      if (gSystem->AccessPathName(rf))
         continue;
      TFile *f = TFile::Open(rf);
      TTree *t = (TTree *)f->Get("cbmsim");
      if (!t) { f->Close(); continue; }
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("AtPatternEvent*", 1);
      TClonesArray *peArr = nullptr;
      t->SetBranchAddress("AtPatternEvent", &peArr);
      Long64_t nE = t->GetEntries();
      if (nEvents > 0 && nEvents < nE)
         nE = nEvents;
      for (Long64_t i = 0; i < nE; ++i) {
         t->GetEntry(i);
         if (peArr->GetEntries() == 0)
            continue;
         auto *pe = (AtPatternEvent *)peArr->At(0);
         if (!pe)
            continue;
         for (auto &trk : pe->GetTrackCand()) {
            ++nTracks;
            AtTrack &tr = const_cast<AtTrack &>(trk);
            auto r = spy.Estimate(tr);
            if (!r.valid)
               continue;
            if (r.arclength < arclenMin)
               continue;
            ++nValid;
            hBt->Fill(r.polar * TMath::RadToDeg(), r.brho);
         }
      }
      f->Close();
      ++nRuns;
   }
   printf("\n== 15C Brho-vs-theta ==  runs=%ld  tracks=%ld  valid=%ld (%.1f%%)  arclenMin=%.0f\n", nRuns, nTracks,
          nValid, nTracks ? 100.0 * nValid / nTracks : 0, arclenMin);

   TCanvas *c = new TCanvas("c", "brho_theta", 900, 700);
   c->SetRightMargin(0.13);
   hBt->Draw("colz");
   TString png = plotDir + "brho_theta_C15" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
