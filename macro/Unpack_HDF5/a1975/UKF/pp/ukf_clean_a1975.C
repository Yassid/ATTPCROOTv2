/// @file ukf_clean_a1975.C
/// @brief Clean UKF proton kinematics: UKF-fitted KE vs theta for IC-gated 16C beam
///        events, proton-PID-gated, good fits.
///
/// Self-contained in <run>_ukf.root (AtTrackingEvent: keeps the original PRA tracks
/// via GetTrackArray() AND the UKF fits via GetFittedTracks()) + <run>_FRIB.root (IC).
/// Each fitted track is linked to its source PRA track by GetTrackID(); we re-run the
/// default (Spyral) PID estimator on that track to apply the proton gate
/// (proton_band.json), the event's ion chamber for the 16C beam gate, and a chi2/ndf
/// cut on the fit. Draws the clean UKF KE-vs-theta and the KE spectrum.
///
///   root -b -q 'ukf_clean_a1975.C("run_0106,...","/mnt/f/a1975/reco/")'

#include <map>
#include <vector>

void ukf_clean_a1975(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,"
                                       "run_0114,run_0115",
                    TString inDir = "/mnt/f/a1975/reco/", TString gateFile = "pid/proton_band.json", double icMin = 950,
                    double icMax = 1350, int icTbLo = 1000, int icTbHi = 1350, double chi2Cut = 5.0, double bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH2F *hkin = new TH2F("hkin", "UKF clean proton kinematics (IC+PID+fit);#theta_{lab} [deg];KE [MeV]", 180, 0, 180,
                         250, 0, 25);
   TH1F *hke = new TH1F("hke", "UKF proton KE (clean);KE [MeV];protons", 200, 0, 25);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nFit = 0, nIC = 0, nProton = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString uf = inDir + run + "_ukf.root", ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(uf) || gSystem->AccessPathName(ff)) {
         printf("skip %s (missing ukf or FRIB)\n", run.Data());
         continue;
      }
      TFile *fu = TFile::Open(uf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      TClonesArray *te = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      TFile *fc = TFile::Open(ff);
      TTree *tc = (TTree *)fc->Get("cbmsim");
      TClonesArray *re = nullptr;
      tc->SetBranchAddress("AtRawEvent", &re);

      Long64_t N = std::min(tu->GetEntries(), tc->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         tc->GetEntry(i);
         double ic = -1;
         if (re->GetEntries() > 0) {
            auto *raw = (AtRawEvent *)re->At(0);
            if (raw && !raw->GetGenTraces().empty()) {
               auto &adc = raw->GetGenTraces()[0]->GetADC();
               double mx = -1e9;
               for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
                  mx = std::max(mx, adc[b]);
               ic = mx;
            }
         }
         tu->GetEntry(i);
         if (te->GetEntries() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         bool icOK = (ic >= icMin && ic <= icMax);
         // map PRA trackID -> original AtTrack (for the PID gate)
         std::vector<AtTrack> orig = ev->GetTrackArray();
         std::map<int, AtTrack *> byID;
         for (auto &tr : orig)
            byID[tr.GetTrackID()] = &tr;

         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            ++nFit;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double c2n = ndf > 0 ? chi2 / ndf : 1e9;
            double ke = k.kineticEnergy, th = k.theta * TMath::RadToDeg();
            if (ke <= 0 || ke > 1000 || c2n > chi2Cut)
               continue;
            if (!icOK)
               continue;
            ++nIC;
            // proton PID gate on the source track
            auto it = byID.find(ft->GetTrackID());
            if (it == byID.end())
               continue;
            auto r = spy.Estimate(*it->second);
            if (!r.valid || !pid.IsInside(r.sqrtdEdx, r.brho))
               continue;
            ++nProton;
            hkin->Fill(th, ke);
            hke->Fill(ke);
         }
      }
      fu->Close();
      fc->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nUKF good fits=%ld  +IC-gate=%ld  +proton-gate=%ld\n", nFit, nIC, nProton);

   TCanvas *c = new TCanvas("c", "ukfclean", 1300, 560);
   c->Divide(2, 1);
   c->cd(1);
   hkin->Draw("colz");
   c->cd(2);
   hke->Draw();
   c->SaveAs("pp/plots/ukf_clean.png");
   printf("saved ukf_clean.png\n");
}
