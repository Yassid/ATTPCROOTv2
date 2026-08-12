/// @file points_data_vs_sim.C
/// @brief Points per track in the DATA against the SIMULATION, at the reconstruction level.
///
/// The fMinPoints = 30 cut inside AtSpyralPID is identical for data and simulation -- neither
/// AtPIDTask nor AtGenfitter ever calls SetMinPoints. It costs the simulation 10.1 % of its tracks.
/// If it cost the data the same, the data's forward-theta_cm yield would be suppressed the same
/// way, and it is not: theta_cm 10-15 deg holds the largest peak in the measured distribution
/// while the simulation calls that region 6 % efficient. Either the data's point clouds are richer
/// than the simulation's, or that reading is wrong. This measures it.
///
/// HOW. AtSpyralPID::Estimate() is run in memory over AtPatternEvent, exactly as AtPIDTask drives
/// it, with fMinPoints forced to 1 so the points count is actually assigned. nPoints, like polar
/// and radius, is only written in the final block of Estimate() (AtSpyralPID.cxx:465), so any
/// track that exits early carries nPoints = 0 and would otherwise pile up at zero and be read as
/// "no points". Tracks that still fail at fMinPoints = 1 are counted separately, not plotted.
///
/// NOT A LIKE-FOR-LIKE PHYSICS SAMPLE. The simulation is one recoil proton per event; the data
/// carries beam, reaction products and whatever else the track finder returned, with no IC gate
/// applied here. This compares the point clouds the finder hands to the PID, which is the quantity
/// the cut acts on -- it is not a statement about relative proton yields.
///
///   root -b -q 'points_data_vs_sim.C()'

void points_data_vs_sim(TString simDir = "/mnt/f/a1975_C16_pp_pid",
                        TString simTags = "s2001,s2002,s2003",
                        TString datDir = "/mnt/f/a1975/reco",
                        TString datRuns = "run_0106,run_0107,run_0108",
                        Long64_t maxEvt = 4000, Double_t bField = 2.85,
                        TString png = "plots/points_data_vs_sim.png")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   spy.SetMinPoints(1); // so nPoints is assigned; the real cut is applied offline below

   auto *hPsim = new TH1D("hPsim", "", 120, 0, 240);
   auto *hPdat = new TH1D("hPdat", "", 120, 0, 240);
   auto *hCsim = new TH1D("hCsim", "", 120, 0, 240);
   auto *hCdat = new TH1D("hCdat", "", 120, 0, 240);

   auto sweep = [&](TString dir, TString list, TString suffix, TH1D *hP, TH1D *hC, const char *what) {
      long nTrk = 0, nEarly = 0;
      TObjArray *ta = list.Tokenize(",");
      for (int it = 0; it < ta->GetEntries(); ++it) {
         TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
         TString fn = dir + "/" + tg + suffix;
         if (gSystem->AccessPathName(fn)) { printf("  skip %s\n", fn.Data()); continue; }
         TFile *f = TFile::Open(fn);
         TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
         if (!t) { if (f) f->Close(); continue; }
         // the data reco file also carries the raw/event branches; only the pattern is needed
         t->SetBranchStatus("*", 0);
         t->SetBranchStatus("AtPatternEvent*", 1);
         TClonesArray *pa = nullptr;
         t->SetBranchAddress("AtPatternEvent", &pa);
         Long64_t N = std::min(t->GetEntries(), maxEvt);
         for (Long64_t i = 0; i < N; ++i) {
            t->GetEntry(i);
            if (!pa || !pa->GetEntriesFast()) continue;
            auto *pe = (AtPatternEvent *)pa->At(0);
            if (!pe) continue;
            for (auto &track : pe->GetTrackCand()) {
               AtTrack &tr = const_cast<AtTrack &>(track);
               auto res = spy.Estimate(tr);
               ++nTrk;
               hC->Fill(res.nClusters);
               if (!res.valid) { ++nEarly; continue; }
               hP->Fill(res.nPoints);
            }
         }
         f->Close();
         printf("  %-10s %-8s tracks so far %ld\n", what, tg.Data(), nTrk);
      }
      delete ta;
      printf("  %s: %ld tracks, %ld (%.1f %%) exit early even at fMinPoints=1 (no nPoints)\n\n", what, nTrk, nEarly,
             100.0 * nEarly / std::max(1L, nTrk));
      return nTrk;
   };

   printf("\n");
   sweep(simDir, simTags, "_reco.root", hPsim, hCsim, "SIM");
   sweep(datDir, datRuns, "_reco.root", hPdat, hCdat, "DATA");

   auto frac = [](TH1D *h, double thr) {
      double tot = h->Integral(1, h->GetNbinsX());
      if (tot <= 0) return 0.0;
      return 100.0 * h->Integral(1, h->FindBin(thr) - 1) / tot;
   };
   printf("  points per track (tracks with a usable estimate)\n");
   printf("    SIM  mean %6.1f   median-ish bin %5.0f\n", hPsim->GetMean(), hPsim->GetBinCenter(hPsim->GetMaximumBin()));
   printf("    DATA mean %6.1f   median-ish bin %5.0f\n\n", hPdat->GetMean(), hPdat->GetBinCenter(hPdat->GetMaximumBin()));
   printf("  clusters per track\n");
   printf("    SIM  mean %6.1f\n    DATA mean %6.1f\n\n", hCsim->GetMean(), hCdat->GetMean());
   printf("  fraction BELOW a given fMinPoints  (what the cut would throw away)\n");
   printf("    thr      SIM      DATA\n");
   for (double thr : {30., 25., 20., 15., 10.})
      printf("    %3.0f    %5.1f %%   %5.1f %%\n", thr, frac(hPsim, thr), frac(hPdat, thr));

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TCanvas *c = new TCanvas("cPD", "points per track", 1400, 600);
   c->Divide(2, 1);
   auto style = [](TH1D *h, int col, const char *ti) {
      h->SetLineColor(col); h->SetLineWidth(2); h->SetTitle(ti);
      if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
   };
   c->cd(1);
   gPad->SetGridy();
   style(hPsim, kAzure + 2, "points per track;n points used by AtSpyralPID;fraction of tracks");
   style(hPdat, kRed + 1, "");
   hPsim->Draw("hist"); hPdat->Draw("hist same");
   auto *l30 = new TLine(30, 0, 30, hPsim->GetMaximum() * 1.05);
   l30->SetLineStyle(2); l30->SetLineColor(kBlack); l30->Draw();
   auto *l15 = new TLine(15, 0, 15, hPsim->GetMaximum() * 1.05);
   l15->SetLineStyle(3); l15->SetLineColor(kGray + 2); l15->Draw();
   auto *lg = new TLegend(0.5, 0.7, 0.88, 0.88);
   lg->AddEntry(hPsim, "simulation", "l");
   lg->AddEntry(hPdat, "data", "l");
   lg->AddEntry(l30, "fMinPoints = 30", "l");
   lg->Draw();
   c->cd(2);
   gPad->SetGridy();
   style(hCsim, kAzure + 2, "clusters per track;n clusters;fraction of tracks");
   style(hCdat, kRed + 1, "");
   hCsim->Draw("hist"); hCdat->Draw("hist same");

   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(here + "/" + png);
   TString ro = png; ro.ReplaceAll(".png", ".root");
   TFile fo(here + "/" + ro, "RECREATE");
   hPsim->Write(); hPdat->Write(); hCsim->Write(); hCdat->Write();
   fo.Close();
   printf("\n  wrote %s\n\n", png.Data());
}
