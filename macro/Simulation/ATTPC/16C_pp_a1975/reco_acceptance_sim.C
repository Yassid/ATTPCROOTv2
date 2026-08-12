/// @file reco_acceptance_sim.C
/// @brief 16C(p,p) acceptance at the RECONSTRUCTION level -- no fitting anywhere in the chain.
///
/// Inputs are <tag>_sim.root (truth) and <tag>_pid.root, the output of pidPass_a1975.C, which runs
/// AtPIDTask on the pattern tracks and never calls a fitter. Nothing here reads a genfitter file.
///
/// TWO CURVES, because "reconstructed" has two meanings and they should not be merged:
///   found : a Spyral entry exists whose polar angle matches the true proton. AtSpyralPID writes an
///           entry for every pattern track including the ones it rejects, so this measures whether
///           the track finder produced the proton's track at all.
///   valid : that entry additionally has valid == true, i.e. the Spyral analysis completed and the
///           track has a usable (dE/dx, Brho). This is what a PID gate can act on.
/// The difference between them is the cost of the PID analysis itself, separate from track finding.
///
/// HONEST LIMIT OF `found`. polar is filled during the Spyral analysis, so a track rejected at an
/// early exit can carry polar = 0 and will not match any truth angle. Such entries are counted and
/// printed; if that count is large, `found` is a lower limit on track finding rather than a
/// measurement of it. Reading it without reading that number would overstate the result.
///
/// theta_cm uses the same exact two-body chain as pp/angdist_gs_pp.C at the SAME Ebeam, so this
/// histogram shares an axis with the data yield.
///
///   root -b -q 'reco_acceptance_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

static double rc_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double rc_thetacm(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a4 = (std::cos(thetalab) * rc_omega2(s, m1 * m1, m2 * m2) * rc_omega2(u, m2 * m2, m3 * m3) -
                (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                  (2 * m2 * m2) +
               s + u - m2 * m2;
   if (a4 < 0)
      return std::nan("");
   double m4x = std::sqrt(a4);
   double t = m2 * m2 + m4x * m4x - 2 * m2 * Et4;
   double c = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
               (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
              (rc_omega2(s, m1 * m1, m2 * m2) * rc_omega2(s, m3 * m3, m4x * m4x));
   if (c < -1 || c > 1)
      return std::nan("");
   return (TMath::Pi() - std::acos(c)) * TMath::RadToDeg();
}

void reco_acceptance_sim(TString dir = "/mnt/f/a1975_C16_pp_pid",
                         TString tags = "s2001,s2002,s2003,s2004,s2005,s2006", Double_t Ebeam = 188.0,
                         Double_t dThetaMax = 10.0, Int_t nBins = 18, TString png = "plots/reco_acceptance_sim.png")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401, m_C16 = 16.0147 * u, m_p = 1.007825 * u;

   auto *hGen = new TH1D("hGen", "", nBins, 0, 180);
   auto *hFnd = new TH1D("hFnd", "", nBins, 0, 180);
   auto *hVal = new TH1D("hVal", "", nBins, 0, 180);
   hGen->Sumw2(); hFnd->Sumw2(); hVal->Sumw2();
   long nG = 0, nF = 0, nV = 0, nEntries = 0, nZeroPolar = 0;

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
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0) continue;
            keT = std::sqrt(pp * pp + m_p * m_p) - m_p;
            thT = std::acos(pz / pp);
            break;
         }
         if (keT <= 0) continue;
         double cm = rc_thetacm(m_C16, m_p, m_p, m_C16, Ebeam, thT, keT);
         if (std::isnan(cm)) continue;
         ++nG; hGen->Fill(cm);

         bool found = false, valid = false;
         if (pe && pe->GetEntriesFast()) {
            auto *ev = (AtPIDEvent *)pe->At(0);
            if (ev)
               for (auto &sp : ev->GetSpyral()) {
                  ++nEntries;
                  if (sp.polar <= 0) { ++nZeroPolar; continue; }
                  // HANDEDNESS. The reconstructed polar comes out as 180 - theta_true, measured
                  // against the opposite z sense (the simulation reverses drift z in digitization).
                  // Verified on s2001: truth 63.67 -> polar 118.05, 33.64 -> 147.17, 74.89 ->
                  // 107.20, residuals 0.1-2 deg. Matching sp.polar directly against theta_true
                  // gives an acceptance of exactly zero in every bin, which is how this was found.
                  double d = std::fabs((180.0 - sp.polar * TMath::RadToDeg()) - thT * TMath::RadToDeg());
                  if (d > dThetaMax) continue;
                  found = true;
                  if (sp.valid) valid = true;
               }
         }
         if (found) { ++nF; hFnd->Fill(cm); }
         if (valid) { ++nV; hVal->Fill(cm); }
      }
      Fs->Close(); Fp->Close();
      printf("  %-8s done\n", tg.Data());
   }
   delete ta;
   if (!nG) { printf("\033[1;31mno truth protons found\033[0m\n"); return; }

   printf("\n  generated %ld\n", nG);
   printf("  track found (polar match)   %ld  (%.1f %%)\n", nF, 100.0 * nF / nG);
   printf("  valid Spyral result         %ld  (%.1f %%)\n", nV, 100.0 * nV / nG);
   printf("  Spyral entries seen %ld, of which polar == 0 (unmatchable): %ld (%.1f %%)\n\n", nEntries, nZeroPolar,
          100.0 * nZeroPolar / std::max(1L, nEntries));

   auto mk = [&](TH1D *n, const char *nm, int col, int mrk) {
      auto *e = (TH1D *)n->Clone(nm);
      e->Divide(n, hGen, 1, 1, "B");
      e->SetDirectory(nullptr);
      e->SetTitle("^{16}C(p,p) reconstruction acceptance (no fit);#theta_{cm} [deg];accepted / generated");
      e->SetMinimum(0); e->SetMaximum(1.05);
      e->SetLineColor(col); e->SetMarkerColor(col); e->SetMarkerStyle(mrk); e->SetLineWidth(2);
      return e;
   };
   TH1D *aF = mk(hFnd, "hAcc_found", kGreen + 3, 24);
   TH1D *aV = mk(hVal, "hAcc_valid", kAzure + 2, 20);

   printf("  theta_cm      gen     track found        valid PID\n");
   for (int b = 1; b <= nBins; ++b) {
      double g = hGen->GetBinContent(b);
      if (g < 1) continue;
      printf("  %3.0f-%-3.0f %8.0f   %5.1f +- %3.1f %%    %5.1f +- %3.1f %%\n", hGen->GetBinLowEdge(b),
             hGen->GetBinLowEdge(b) + hGen->GetBinWidth(b), g, 100 * aF->GetBinContent(b), 100 * aF->GetBinError(b),
             100 * aV->GetBinContent(b), 100 * aV->GetBinError(b));
   }

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TCanvas *c = new TCanvas("cR", "reco acceptance", 900, 650);
   gPad->SetGridy();
   aF->Draw("E1"); aV->Draw("E1 same");
   auto *lg = new TLegend(0.15, 0.15, 0.58, 0.30);
   lg->AddEntry(aF, "track found (no fit)", "lp");
   lg->AddEntry(aV, "+ valid Spyral PID", "lp");
   lg->Draw();
   gSystem->mkdir(here + "/" + gSystem->DirName(png), kTRUE);
   c->SaveAs(here + "/" + png);
   TString ro = png; ro.ReplaceAll(".png", ".root");
   TFile fo(here + "/" + ro, "RECREATE");
   hGen->Write(); hFnd->Write(); hVal->Write(); aF->Write(); aV->Write();
   fo.Close();
   printf("\n  wrote %s\n         %s\n\n", png.Data(), ro.Data());
}
