/// @file pid_C14.C
/// @brief PID plane (sqrt(dEdx) vs Brho) for a1954 14C(p,p') from the faithful Spyral
///        estimator (AtTools::AtSpyralPID), computed on the PRA track candidates in
///        <run>_reco.root. NO fit needed — this is the step to inspect BEFORE fitting.
///        Spyral PID is charge-sign agnostic, so it is unaffected by the Bz handedness.
///
///   root -b -q 'pid/pid_C14.C("run_0055")'
///   root -b -q 'pid/pid_C14.C("run_0055,run_0058,run_0060", "/home/yassid/a1954_C14_reco/")'

// Parse the "vertices":[[x,y],...] array of an AtParticleID JSON gate into a TGraph
// (closed), so we can overlay the reference band. Minimal reader (no JSON lib).
static TGraph *LoadGateGraph(const char *path)
{
   std::ifstream in(path);
   if (!in)
      return nullptr;
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto pos = s.find("vertices");
   if (pos == std::string::npos)
      return nullptr;
   pos = s.find('[', pos); // opening bracket of the array
   auto endArr = s.find("]\n", pos);
   std::vector<double> nums;
   const char *p = s.c_str() + pos;
   const char *end = s.c_str() + s.size();
   while (p < end) {
      char *np = nullptr;
      double v = strtod(p, &np);
      if (np == p) {
         ++p;
         continue;
      }
      nums.push_back(v);
      p = np;
      if (s.find("vertices") != std::string::npos && (size_t)(p - s.c_str()) > endArr + 2 && endArr != std::string::npos)
         break;
   }
   auto *g = new TGraph();
   for (size_t i = 0; i + 1 < nums.size(); i += 2)
      g->SetPoint(g->GetN(), nums[i], nums[i + 1]);
   if (g->GetN() > 0) {
      double x0, y0;
      g->GetPoint(0, x0, y0);
      g->SetPoint(g->GetN(), x0, y0); // close it
   }
   g->SetLineColor(kRed);
   g->SetLineWidth(2);
   return g;
}

void pid_C14(TString runsCSV = "run_0055", TString inDir = "/home/yassid/a1954_C14_reco/", Long64_t nEvents = -1,
              double bField = 2.85, TString outTag = "",
              TString refGate = "/home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_HDF5/a1975/UKF/pid/proton_band.json",
              double sqrtMax = 80, double brhoMax = 2.5, double dedxMax = 4000)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   // Axis ranges are parameters (sqrtMax/brhoMax/dedxMax); red overlay = 16C+p proton band.
   TH2F *hS = new TH2F("hS", "14C Spyral PID (red = 16C+p proton band);#sqrt{dEdx};B#rho [T m]", 400, 0, sqrtMax, 400,
                       0, brhoMax);
   TH2F *hSd = new TH2F("hSd", "14C Spyral PID;dEdx [counts];B#rho [T m]", 400, 0, dedxMax, 400, 0, brhoMax);
   TH2F *hKt = new TH2F("hKt", "dEdx vs #theta_{lab};#theta_{lab} [deg];dEdx [counts]", 180, 0, 180, 400, 0, dedxMax);

   TFile *fo = new TFile((plotDir + "pid_C14" + outTag + ".root").Data(), "RECREATE");
   TNtuple *nt = new TNtuple("spid", "spyral pid", "dedx:sqrtdedx:brho:polar:arclen:npts:vtxr:vtxz:radius");

   TObjArray *runs = runsCSV.Tokenize(",");
   long nTracks = 0, nValid = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString rf = inDir + run + "_reco.root";
      if (gSystem->AccessPathName(rf)) {
         printf("skip %s (no %s)\n", run.Data(), rf.Data());
         continue;
      }
      TFile *f = TFile::Open(rf);
      TTree *t = (TTree *)f->Get("cbmsim");
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
            ++nValid;
            double th = r.polar * TMath::RadToDeg();
            hS->Fill(r.sqrtdEdx, r.brho);
            hSd->Fill(r.dEdx, r.brho);
            hKt->Fill(th, r.dEdx);
            double vtxr = TMath::Hypot(r.vertex.X(), r.vertex.Y());
            nt->Fill(r.dEdx, r.sqrtdEdx, r.brho, th, r.arclength, r.nPoints, vtxr, r.vertex.Z(), r.radius);
         }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }
   fo->cd();
   nt->Write();
   printf("\n===== 14C Spyral PID  [%s]  =====\n", runsCSV.Data());
   printf("tracks=%ld   Spyral-valid=%ld (%.1f%%)\n", nTracks, nValid, nTracks ? 100.0 * nValid / nTracks : 0);

   TGraph *gGate = LoadGateGraph(refGate.Data());
   if (gGate)
      printf("overlaying reference gate %s (%d vertices)\n", refGate.Data(), gGate->GetN() - 1);

   TCanvas *c = new TCanvas("c", "pid", 1500, 500);
   c->Divide(3, 1);
   c->cd(1);
   hS->Draw("colz");
   if (gGate)
      gGate->Draw("L same");
   c->cd(2);
   hSd->Draw("colz");
   c->cd(3);
   hKt->Draw("colz");
   TString png = plotDir + "pid_C14" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
   fo->Close();
}
