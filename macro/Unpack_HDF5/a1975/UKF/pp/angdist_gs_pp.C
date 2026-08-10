/// @file angdist_gs_pp.C
/// @brief Ground-state angular distribution of 16C(p,p), selected by the hand-drawn Ex-theta_cm cut.
///
/// Uses exactly the plane the cut was drawn on -- same beam energy, same theta correction, same
/// chi2 -- so the polygon means what it meant when it was drawn. Change any of them and the cut
/// no longer corresponds to the data it is applied to.
///
///     Ebeam   188 MeV,  theta -> theta - (360/kcDenom)*(KE - kcPivot),  chi2/ndf < 5
///     cut     plots/ex_cut_gs.root   (TCutG in theta_cm vs Ex)
///
/// WHAT THIS DOES AND DOES NOT CORRECT. Counts are divided by the solid angle of the bin,
/// 2*pi*(cos(theta_lo) - cos(theta_hi)), which is not constant and whose omission would tilt the
/// whole distribution by sin(theta). NO acceptance correction is applied: no acceptance has been
/// simulated for a1975 (p,p). The result is therefore a SHAPE in arbitrary units, and it is a
/// shape that still carries the detector's angular efficiency.
///
/// The yield in each bin is simply the number of tracks inside the polygon. That is a valid
/// measure of the ground state only where the polygon contains the whole peak and little else;
/// the printed table gives the count per bin so thin bins are visible rather than implied.
///
///   root -b -q 'angdist_gs_pp.C("/path/pp_kin.root")'

void angdist_gs_pp(TString cache, TString cutFile = "plots/ex_cut_gs.root", TString cutName = "gs",
                   Double_t Ebeam = 188.0, Double_t kcDenom = 4000.0, Double_t kcPivot = 1.5,
                   Double_t chi2Max = 5.0, Double_t icLo = 1000, Double_t icHi = 1350, Double_t cmLo = 0,
                   Double_t cmHi = 150, Double_t dcm = 5.0, TString tag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double slope = 360.0 / kcDenom;

   TFile *fc = TFile::Open(here + "/" + cutFile);
   TCutG *cut = fc && !fc->IsZombie() ? (TCutG *)fc->Get(cutName) : nullptr;
   if (!cut) {
      printf("\033[1;31mno TCutG '%s' in %s\033[0m\n", cutName.Data(), cutFile.Data());
      return;
   }
   printf("\n  cut '%s': %d vertices\n", cutName.Data(), cut->GetN());

   TFile *f = TFile::Open(cache);
   TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("pk") : nullptr;
   if (!t) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   float ke, th, vz, c2, ic;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vz", &vz);
   t->SetBranchAddress("chi2ndf", &c2);
   t->SetBranchAddress("ic", &ic);

   const double u = 931.49401;
   const double mb = 16.0147013 * u, mt = 1.00782503 * u, m3 = 1.00782503 * u, m4 = 16.0147013 * u;
   auto om = [](double x, double y, double z) {
      return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
   };

   const int NB = (int)std::lround((cmHi - cmLo) / dcm);
   auto *hY = new TH1D("hY", "counts;#theta_{cm} [deg];counts", NB, cmLo, cmHi);
   hY->Sumw2();
   auto *hAll = new TH1D("hAll", "", NB, cmLo, cmHi);
   auto *h2 = new TH2D("h2", "selected events;#theta_{cm} [deg];E_{x} [MeV]", NB, cmLo, cmHi, 200, -3, 4);

   long nTot = 0, nIn = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(c2 < chi2Max && ic > icLo && ic < icHi))
         continue;
      double thc = (th - slope * (ke - kcPivot)) * TMath::DegToRad();
      double E1 = Ebeam + mb, E3 = ke + m3, E4 = E1 + mt - E3;
      double s = mb * mb + mt * mt + 2 * mt * E1, uu = mt * mt + m3 * m3 - 2 * mt * E3;
      double a = (std::cos(thc) * om(s, mb * mb, mt * mt) * om(uu, mt * mt, m3 * m3) -
                  (s - mb * mb - mt * mt) * (mt * mt + m3 * m3 - uu)) /
                    (2 * mt * mt) +
                 s + uu - mt * mt;
      if (a < 0)
         continue;
      double m4x = std::sqrt(a), ex = m4x - m4;
      double tt = mt * mt + m4x * m4x - 2 * mt * E4;
      double arg = (s * s + s * (2 * tt - mb * mb - mt * mt - m3 * m3 - m4x * m4x) +
                    (mb * mb - mt * mt) * (m3 * m3 - m4x * m4x)) /
                   (om(s, mb * mb, mt * mt) * om(s, m3 * m3, m4x * m4x));
      if (arg < -1 || arg > 1)
         continue;
      double cm = (TMath::Pi() - std::acos(arg)) * TMath::RadToDeg();
      ++nTot;
      hAll->Fill(cm);
      if (!cut->IsInside(cm, ex))
         continue;
      ++nIn;
      hY->Fill(cm);
      h2->Fill(cm, ex);
   }
   printf("  %ld protons on the plane, %ld inside the cut (%.1f%%)\n\n", nTot, nIn, 100.0 * nIn / std::max(1L, nTot));

   // counts -> counts per steradian. No acceptance: none exists for a1975 (p,p).
   auto *hD = new TH1D("hD", "d#sigma/d#Omega (shape);#theta_{cm} [deg];counts / sr  [arb]", NB, cmLo, cmHi);
   hD->Sumw2();
   printf("  theta_cm |  in cut |  all   | frac  |  dOmega [sr] | counts/sr\n");
   for (int b = 1; b <= NB; ++b) {
      double lo = hY->GetBinLowEdge(b), hi = lo + dcm;
      double n = hY->GetBinContent(b), e = hY->GetBinError(b), all = hAll->GetBinContent(b);
      double dO = 2 * TMath::Pi() * (std::cos(lo * TMath::DegToRad()) - std::cos(hi * TMath::DegToRad()));
      if (n <= 0 || dO <= 0)
         continue;
      hD->SetBinContent(b, n / dO);
      hD->SetBinError(b, e / dO);
      printf("  %3.0f-%3.0f  | %7.0f | %6.0f | %.3f | %12.4f | %10.4g\n", lo, hi, n, all, all > 0 ? n / all : 0, dO,
             n / dO);
   }

   TCanvas *c1 = new TCanvas("cad", "g.s. angular distribution", 1400, 600);
   c1->Divide(2, 1);
   c1->cd(1);
   gPad->SetLogz();
   h2->Draw("colz");
   cut->SetLineColor(kGreen + 2);
   cut->SetLineWidth(3);
   cut->Draw("L same");
   c1->cd(2);
   gPad->SetLogy();
   gPad->SetGridy();
   hD->SetMarkerStyle(20);
   hD->SetMarkerSize(1.2);
   hD->SetLineWidth(2);
   hD->SetTitle(TString::Format("^{16}C(p,p) g.s., E_{beam} = %.0f MeV (shape, no acceptance);"
                                "#theta_{cm} [deg];counts / sr  [arb]",
                                Ebeam));
   hD->Draw("E1");

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/angdist_gs_pp" + tag + ".png";
   c1->SaveAs(png);
   TFile fo(here + "/plots/angdist_gs_pp" + tag + ".root", "RECREATE");
   hY->Write("yield");
   hD->Write("dsdo_shape");
   h2->Write("ex_vs_thcm_selected");
   fo.Close();
   printf("\n  wrote %s\n         plots/angdist_gs_pp%s.root\n\n", png.Data(), tag.Data());
}
