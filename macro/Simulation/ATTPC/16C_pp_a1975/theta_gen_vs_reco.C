/// @file theta_gen_vs_reco.C
/// @brief Generated proton polar angle against the RECONSTRUCTED polar angle, no fitting.
///
/// Truth from <tag>_sim.root, reconstruction from <tag>_pid.root (pidPass_a1975.C -> AtPIDTask on
/// the pattern tracks). Every Spyral entry in the event is plotted against that event's true
/// proton angle, with no truth matching applied -- so the diagonal band is the proton and anything
/// off it is either a second track or a mis-measured one. Filtering first would hide exactly the
/// thing this plot is for.
///
/// The reference line is theta_reco = 180 - theta_true, not the identity: the reconstructed polar
/// is measured against the opposite z sense. Drawing y = x instead makes a correct reconstruction
/// look like a total failure.
///
///   root -b -q 'theta_gen_vs_reco.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void theta_gen_vs_reco(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                       TString png = "plots/theta_gen_vs_reco.png")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double m_p = 938.272;

   auto *h = new TH2D("hGR", "16C(p,p) protons, no fit;#theta_{lab} generated [deg];#theta_{lab} reconstructed [deg]",
                      180, 0, 180, 180, 0, 180);
   long nPair = 0, nZero = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fs = dir + "/" + tg + "_sim.root", fp = dir + "/" + tg + "_pid.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fp)) { printf("  skip %s\n", tg.Data()); continue; }
      TFile *Fs = TFile::Open(fs), *Fp = TFile::Open(fp);
      TTree *ts = Fs ? (TTree *)Fs->Get("cbmsim") : nullptr;
      TTree *tp = Fp ? (TTree *)Fp->Get("cbmsim") : nullptr;
      if (!ts || !tp) { if (Fs) Fs->Close(); if (Fp) Fp->Close(); continue; }
      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tp->SetBranchAddress("AtPIDEvent", &pe);

      Long64_t N = std::min(ts->GetEntries(), tp->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i); tp->GetEntry(i);
         double thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp > 0) thT = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (thT < 0 || !pe || !pe->GetEntriesFast()) continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev) continue;
         for (auto &sp : ev->GetSpyral()) {
            if (sp.polar <= 0) { ++nZero; continue; }
            h->Fill(thT, sp.polar * TMath::RadToDeg());
            ++nPair;
         }
      }
      Fs->Close(); Fp->Close();
      printf("  %-8s done\n", tg.Data());
   }
   delete ta;

   // residual about the 180 - x line, for the entries that sit on it
   double sum = 0, sum2 = 0; long n = 0;
   for (int bx = 1; bx <= h->GetNbinsX(); ++bx)
      for (int by = 1; by <= h->GetNbinsY(); ++by) {
         double c = h->GetBinContent(bx, by);
         if (c <= 0) continue;
         double d = (180.0 - h->GetYaxis()->GetBinCenter(by)) - h->GetXaxis()->GetBinCenter(bx);
         if (std::fabs(d) > 10) continue;
         sum += c * d; sum2 += c * d * d; n += (long)c;
      }
   double mean = n ? sum / n : 0, rms = n ? std::sqrt(sum2 / n - mean * mean) : 0;
   printf("\n  pairs plotted %ld   (polar == 0, dropped: %ld)\n", nPair, nZero);
   printf("  on the 180-x band (|d| < 10 deg): %ld entries, mean %.2f deg, rms %.2f deg\n\n", n, mean, rms);

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TCanvas *c = new TCanvas("cGR", "gen vs reco", 820, 760);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   h->Draw("colz");
   auto *l = new TLine(0, 180, 180, 0); // theta_reco = 180 - theta_gen
   l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->SetLineStyle(2);
   l->Draw();
   auto *lid = new TLine(0, 0, 180, 180); // identity, for contrast
   lid->SetLineColor(kGray + 2); lid->SetLineWidth(1); lid->SetLineStyle(3);
   lid->Draw();
   auto *lg = new TLegend(0.15, 0.78, 0.55, 0.88);
   lg->AddEntry(l, "#theta_{reco} = 180 - #theta_{gen}", "l");
   lg->AddEntry(lid, "identity", "l");
   lg->Draw();
   gSystem->mkdir(here + "/" + gSystem->DirName(png), kTRUE);
   c->SaveAs(here + "/" + png);
   printf("  wrote %s\n\n", png.Data());
}
