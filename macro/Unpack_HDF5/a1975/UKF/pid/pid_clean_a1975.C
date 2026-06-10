/// @file pid_clean_a1975.C
/// @brief Clean, IC-gated PID + proton-gated kinematics on the DEFAULT (Spyral) PID.
///
/// For each run it index-matches <run>_pid.root (AtPIDEvent, GetDefault()=Spyral) to
/// <run>_FRIB.root (AtRawEvent; ion chamber = generic trace[0]). Keeps tracks with
/// the IC amplitude inside the 16C beam window and passing quality cuts, then applies
/// the proton gate (proton_band.json, sqrtdEdx vs brho). Draws: IC spectrum, the clean
/// sqrt(dEdx) PID plane with the proton polygon, and the proton-gated brho-vs-theta
/// kinematics (theta reflected 180-theta for display).
///
///   root -b -q 'pid_clean_a1975.C("run_0106,...","/mnt/f/a1975/reco/")'

#include <vector>

void pid_clean_a1975(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,"
                                       "run_0114,run_0115",
                    TString inDir = "/mnt/f/a1975/reco/", TString gateFile = "pid/proton_band.json", double icMin = 950,
                    double icMax = 1350, int icTbLo = 1000, int icTbHi = 1350, int minPts = 15, double maxVtxR = 40)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data()); // proton gate (sqrtdEdx vs brho)

   TH1F *hic = new TH1F("hic", "IC amplitude (quality-cut tracks);IC amp [ADC];tracks", 250, 0, 2500);
   TH2F *hplane = new TH2F("hplane", "Clean PID (quality+IC) + proton gate;#sqrt{dEdx};B#rho [T m]", 350, 0, 40, 350, 0,
                           1.0);
   TH2F *hkin = new TH2F("hkin", "Proton-gated kinematics (IC+quality);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 350,
                         0, 1.0);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nQual = 0, nIC = 0, nProton = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString pf = inDir + run + "_pid.root", ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(pf) || gSystem->AccessPathName(ff)) {
         printf("skip %s (missing pid or FRIB)\n", run.Data());
         continue;
      }
      TFile *fp = TFile::Open(pf);
      TTree *tp = (TTree *)fp->Get("cbmsim");
      TClonesArray *pe = nullptr;
      tp->SetBranchAddress("AtPIDEvent", &pe);
      TFile *fc = TFile::Open(ff);
      TTree *tc = (TTree *)fc->Get("cbmsim");
      TClonesArray *re = nullptr;
      tc->SetBranchAddress("AtRawEvent", &re);

      Long64_t N = std::min(tp->GetEntries(), tc->GetEntries());
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
         tp->GetEntry(i);
         if (pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &r : ev->GetDefault()) {
            if (!r.valid)
               continue;
            double vtxr = std::sqrt(r.vertex.X() * r.vertex.X() + r.vertex.Y() * r.vertex.Y());
            if (r.nPoints < minPts || vtxr > maxVtxR)
               continue; // quality
            ++nQual;
            hic->Fill(ic);
            if (ic < icMin || ic > icMax)
               continue; // IC gate (16C beam)
            ++nIC;
            hplane->Fill(r.sqrtdEdx, r.brho);
            if (pid.IsInside(r.sqrtdEdx, r.brho)) { // proton gate
               ++nProton;
               hkin->Fill(180.0 - r.polar * TMath::RadToDeg(), r.brho); // reflected theta
            }
         }
      }
      fp->Close();
      fc->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nquality=%ld  +IC-gate=%ld  +proton-gate=%ld\n", nQual, nIC, nProton);

   TCanvas *c = new TCanvas("c", "clean", 1550, 520);
   c->Divide(3, 1);
   c->cd(1);
   gPad->SetLogy();
   hic->Draw();
   TLine *l1 = new TLine(icMin, 0, icMin, hic->GetMaximum());
   TLine *l2 = new TLine(icMax, 0, icMax, hic->GetMaximum());
   l1->SetLineColor(kRed);
   l2->SetLineColor(kRed);
   l1->Draw();
   l2->Draw();
   c->cd(2);
   gPad->SetLogz();
   hplane->Draw("colz");
   c->cd(3);
   gPad->SetLogz();
   hkin->Draw("colz");
   c->SaveAs("pid/plots/pid_clean.png");
   printf("saved pid_clean.png\n");
}
