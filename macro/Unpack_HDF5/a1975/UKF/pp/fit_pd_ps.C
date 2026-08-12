/// @file fit_pd_ps.C
/// @brief 16C(p,d)15C excitation spectrum with the 1n PHASE SPACE included, global then per angle.
///
/// Ports the 15C+p procedure (a2091/macros/C15_pp_ana_AngBins.C) to (p,d). The model is a functor
/// so the phase space can enter as a shape rather than an analytic term:
///     2 Gaussians (the BOUND states: g.s. and 0.740)
///   + 2 Breit-Wigners (the unbound structures near 3.4 and 4.9)
///   + one free scale on the 1n phase space, evaluated from a TGraph
/// Above Sn(15C) = 1.218 MeV the residual is unbound, so what sits under those structures is
/// continuum, not a straight line. Fitting it with pol1 -- as the first version of this analysis
/// did -- forces the Breit-Wigner areas to absorb continuum, and the bound peaks to absorb whatever
/// the background misses at the low-Ex end.
///
/// PHASE SPACE HANDLING, following the 15C+p macro exactly:
///   * PhaseSpace_16C_pd_1n.root, tree simulated_tree, branches Ex_cal / ThetaCM_cal / Weight_sim
///   * filled WEIGHTED, once inclusively and once per theta_cm bin, so each angle bin gets its own
///     phase-space shape -- the shape is angle dependent and a single inclusive one would be wrong
///   * Smooth()ed, then converted to a TGraph. The histogram spans the whole fit range, so it is
///     identically zero where the phase space has no events and the TGraph cannot extrapolate
///     nonsense below threshold.
///
/// Selection is the saved explorer configuration of 2026-08-12: Ebeam 185, chi2/ndf < 5,
/// KE 0-30, vertex z 0-500 mm, theta correction ON (360/2950.8 per MeV about KE = 27).
///
///   root -b -q 'pp/fit_pd_ps.C("/path/pd_kin.root")'

/// Four Gaussians plus the phase space. The unbound structures are Gaussians here rather than
/// Breit-Wigners: with a BW the ~3.4 term collapsed to Gamma = 0.02 (a delta spike leaving the real
/// 0.5 MeV-wide peak unfitted) and the ~4.9 term ran to whatever upper bound it was given, 2.0 then
/// 3.5, i.e. it was being used as extra background rather than as a resonance. A Gaussian cannot do
/// either. It is not the right lineshape for an unbound state, so these two widths are resolution
/// plus intrinsic width mixed together and should not be read as Gamma.
class PDSpectralModel {
public:
   TGraph *g_PS;
   PDSpectralModel(TGraph *g) : g_PS(g) {}
   double operator()(double *x, double *p)
   {
      double v = 0;
      v += p[0] * TMath::Gaus(x[0], p[1], p[2], false);          // g.s.
      v += p[3] * TMath::Gaus(x[0], p[4], p[5], false);          // 0.740
      v += p[6] * TMath::Gaus(x[0], p[7], p[8], false);          // ~3.4
      v += p[9] * TMath::Gaus(x[0], p[10], p[11], false);        // ~4.9
      v += p[12] * g_PS->Eval(x[0]);                             // 1n phase space
      return v;
   }
};

static TGraph *pdHistoToGraph(TH1F *h)
{
   auto *g = new TGraph();
   for (int i = 1; i <= h->GetNbinsX(); ++i) g->SetPoint(i - 1, h->GetBinCenter(i), h->GetBinContent(i));
   return g;
}

static double pd_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static std::pair<double, double> pd_kine(double Kp, double thl, double Ke)
{
   const double u = 931.49401;
   double m1 = 16.014701 * u, m2 = 1.007825 * u, m3 = 2.013553 * u, m4 = 15.010599 * u;
   double E1 = Kp + m1, E3 = Ke + m3, E4 = E1 + m2 - E3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, uu = m2 * m2 + m3 * m3 - 2 * m2 * E3;
   double a = (std::cos(thl) * pd_om2(s, m1 * m1, m2 * m2) * pd_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                 (2 * m2 * m2) +
              s + uu - m2 * m2;
   if (a < 0) return {std::nan(""), std::nan("")};
   double m4x = std::sqrt(a), ex = m4x - m4;
   double t = m2 * m2 + m4x * m4x - 2 * m2 * E4;
   double c = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
               (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
              (pd_om2(s, m1 * m1, m2 * m2) * pd_om2(s, m3 * m3, m4x * m4x));
   if (c < -1 || c > 1) return {std::nan(""), std::nan("")};
   return {ex, (TMath::Pi() - std::acos(c)) * TMath::RadToDeg()};
}

void fit_pd_ps(TString cache,
               TString psFile = "/mnt/c/Users/Yassid/Downloads/PhaseSpace_16C/PhaseSpace_16C_pd_1n.root",
               Double_t Ebeam = 185.0, Double_t kcDenom = 2950.8, Double_t kcPivot = 27.0,
               Double_t chi2max = 5.0, Double_t keLo = 0, Double_t keHi = 30, Double_t vzLo = 0,
               Double_t vzHi = 500, Double_t icMin = 950, Double_t icMax = 1350, Double_t cmLo = 15,
               Double_t cmHi = 80, Double_t dcm = 5.0, Double_t exLo = -1.0, Double_t exHi = 6.5,
               Int_t nEx = 150, Int_t slRebin = 2, Double_t lumiFull = 316.4, Double_t zFull = 940.0,
               TString tag = "ps20260812", Bool_t accWeightPS = kTRUE)
{
   gStyle->SetOptStat(0);
   ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double slope = 360.0 / kcDenom;
   const int NB = (int)std::lround((cmHi - cmLo) / dcm);

   // ---------------- data ----------------
   TFile *f = TFile::Open(cache);
   TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("pk") : nullptr;
   if (!t) { printf("\033[1;31mcannot open %s\033[0m\n", cache.Data()); return; }
   float ke, th, vz, c2, ic;
   t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vz", &vz);
   t->SetBranchAddress("chi2ndf", &c2); t->SetBranchAddress("ic", &ic);
   auto *hEx = new TH1F("hEx", "", nEx, exLo, exHi); hEx->SetDirectory(nullptr); hEx->Sumw2();
   // The slices carry ~1/13 of the statistics of the inclusive spectrum, so they are binned
   // coarser by slRebin. The phase-space slices are booked on the SAME axis further down -- if the
   // two ever differ, TGraph::Eval silently interpolates the continuum onto the wrong grid.
   const int nExSl = nEx / std::max(1, (int)slRebin);
   std::vector<TH1F *> hSl(NB);
   for (int k = 0; k < NB; ++k) {
      hSl[k] = new TH1F(Form("hSl%d", k), "", nExSl, exLo, exHi);
      hSl[k]->SetDirectory(nullptr); hSl[k]->Sumw2();
   }
   long n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(c2 < chi2max && ic > icMin && ic < icMax)) continue;
      if (!(ke > keLo && ke < keHi)) continue;
      if (!(vz > vzLo && vz < vzHi)) continue;
      double thc = (th - slope * (ke - kcPivot)) * TMath::DegToRad();
      auto [ex, tcm] = pd_kine(Ebeam, thc, ke);
      if (std::isnan(ex) || tcm < cmLo || tcm > cmHi) continue;
      ++n; hEx->Fill(ex);
      int k = (int)((tcm - cmLo) / dcm);
      if (k >= 0 && k < NB) hSl[k]->Fill(ex);
   }
   f->Close();
   printf("\n  %ld deuterons in theta_cm %.0f-%.0f after the saved selection\n", n, cmLo, cmHi);

   // ---------------- acceptance surface, for weighting the phase space ----------------
   // The phase space is GENERATOR-LEVEL: its theta_cm shape is a smooth hump peaking at 85-90 deg
   // with none of the acceptance's structure, and it is exactly zero above 155 deg where the
   // detector still accepts 0.43. So it has to be weighted by acceptance(Ex, theta_cm) before it
   // can be compared with data. Unweighted it overshoots the near-threshold gap by 1.5-1.7x while
   // being right at high Ex -- which is the acceptance falling off toward threshold.
   // The four simulated states sample the Ex axis at 0.000 / 0.740 / 3.100 / 4.660; interpolate
   // linearly between them in Ex, and hold the end values beyond.
   const double accEx[4] = {0.000, 0.740, 3.100, 4.660};
   const char *accT[4] = {"gs", "ex1", "ex2", "ex3"};
   TH1D *ACC[4] = {nullptr, nullptr, nullptr, nullptr};
   for (int g = 0; g < 4; ++g) {
      TFile *fa = TFile::Open(here + Form("/../../../../Simulation/ATTPC/16C_pd/diagnostics/acceptance_%s.root", accT[g]));
      if (fa && !fa->IsZombie()) {
         TH1D *A = (TH1D *)fa->Get(Form("hAcc_%s", accT[g]));
         if (A) { ACC[g] = (TH1D *)A->Clone(Form("ACCs%d", g)); ACC[g]->SetDirectory(nullptr); }
         fa->Close();
      }
      if (!ACC[g]) printf("\033[1;31m  missing acceptance for %s -- PS weighting disabled\033[0m\n", accT[g]);
   }
   bool haveAcc = ACC[0] && ACC[1] && ACC[2] && ACC[3];
   auto accAt = [&](double ex, double tc) {
      if (!haveAcc) return 1.0;
      int lo = 0;
      while (lo < 2 && ex > accEx[lo + 1]) ++lo;
      int hi = lo + 1;
      double w = (ex - accEx[lo]) / (accEx[hi] - accEx[lo]);
      w = std::max(0.0, std::min(1.0, w));
      double a0 = ACC[lo]->GetBinContent(ACC[lo]->FindBin(tc));
      double a1 = ACC[hi]->GetBinContent(ACC[hi]->FindBin(tc));
      return (1 - w) * a0 + w * a1;
   };

   // ---------------- phase space ----------------
   TFile *fps = TFile::Open(psFile);
   TTree *tps = fps && !fps->IsZombie() ? (TTree *)fps->Get("simulated_tree") : nullptr;
   if (!tps) { printf("\033[1;31mcannot open the phase space %s\033[0m\n", psFile.Data()); return; }
   double W, Eps, Tps;
   tps->SetBranchAddress("Weight_sim", &W);
   tps->SetBranchAddress("Ex_cal", &Eps);
   tps->SetBranchAddress("ThetaCM_cal", &Tps);
   auto *hPSfull = new TH1F("hPSfull", "", nEx, exLo, exHi); hPSfull->SetDirectory(nullptr);
   auto *hPSraw = new TH1F("hPSraw", "", nEx, exLo, exHi); hPSraw->SetDirectory(nullptr); // before weighting
   std::vector<TH1F *> hPS(NB);
   for (int k = 0; k < NB; ++k) {
      hPS[k] = new TH1F(Form("hPS%d", k), "", nExSl, exLo, exHi); hPS[k]->SetDirectory(nullptr);
   }
   const Long64_t nPS = tps->GetEntries(); // read BEFORE the file is closed -- the tree dies with it
   for (Long64_t i = 0; i < nPS; ++i) {
      tps->GetEntry(i);
      hPSraw->Fill(Eps, W);
      double wa = (accWeightPS && haveAcc) ? accAt(Eps, Tps) : 1.0;
      hPSfull->Fill(Eps, W * wa);
      int k = (int)((Tps - cmLo) / dcm);
      if (k >= 0 && k < NB) hPS[k]->Fill(Eps, W * wa);
   }
   fps->Close();
   printf("  phase space: %lld entries, integral raw %.4g -> acceptance-weighted %.4g  (x%.3f)\n", nPS,
          hPSraw->Integral(), hPSfull->Integral(),
          hPSraw->Integral() > 0 ? hPSfull->Integral() / hPSraw->Integral() : 0.0);
   if (accWeightPS && haveAcc) {
      printf("\n  phase-space shape, before vs after the acceptance weighting (both normalised)\n");
      printf("   Ex band    raw     weighted   ratio\n");
      double rN = hPSraw->Integral(), wN = hPSfull->Integral();
      for (auto r : std::vector<std::pair<double, double>>{{1.2, 2.0}, {2.0, 2.9}, {2.9, 3.9},
                                                           {3.9, 4.4}, {4.4, 5.6}, {5.6, 6.4}}) {
         double a = hPSraw->Integral(hPSraw->FindBin(r.first), hPSraw->FindBin(r.second)) / rN;
         double b = hPSfull->Integral(hPSfull->FindBin(r.first), hPSfull->FindBin(r.second)) / wN;
         printf("  %4.1f-%4.1f  %6.4f   %6.4f   %5.3f\n", r.first, r.second, a, b, a > 0 ? b / a : 0);
      }
   }
   hPSfull->Smooth();
   for (int k = 0; k < NB; ++k) {
      // a bin with no phase-space events borrows its neighbour's SHAPE (the scale is free anyway);
      // leaving it empty would silently drop the continuum from that angle bin
      if (hPS[k]->Integral() <= 0) {
         int src = (k + 1 < NB && hPS[k + 1]->Integral() > 0) ? k + 1 : (k > 0 ? k - 1 : -1);
         if (src >= 0) { hPS[k]->Add(hPS[src]); printf("  PS bin %d empty -> using bin %d shape\n", k, src); }
      }
      hPS[k]->Smooth();
   }
   TGraph *gPSfull = pdHistoToGraph(hPSfull);

   // ---------------- global fit ----------------
   TF1 *f2g = new TF1("f2g", "gaus(0)+gaus(3)", -0.8, 1.4);
   f2g->SetParameters(hEx->GetMaximum() * 0.7, 0.0, 0.22, hEx->GetMaximum(), 0.78, 0.20);
   hEx->Fit(f2g, "RQ0");
   TF1 *bw1 = new TF1("bw1", "gaus", 2.9, 3.9);
   bw1->SetParameters(100, 3.38, 0.25); hEx->Fit(bw1, "RQ0");
   TF1 *bw2 = new TF1("bw2", "gaus", 4.2, 5.8);
   bw2->SetParameters(30, 4.90, 0.55);  hEx->Fit(bw2, "RQ0+");
   printf("  pre-fit ~4.9 gaussian: A %.1f  E %.3f  sigma %.3f\n", bw2->GetParameter(0), bw2->GetParameter(1),
          bw2->GetParameter(2));

   auto *mFull = new PDSpectralModel(gPSfull);
   TF1 *F = new TF1("F", mFull, exLo, exHi, 13, "PDSpectralModel");
   double ini[13] = {f2g->GetParameter(0), f2g->GetParameter(1), f2g->GetParameter(2),
                     f2g->GetParameter(3), f2g->GetParameter(4), f2g->GetParameter(5),
                     bw1->GetParameter(0), bw1->GetParameter(1), bw1->GetParameter(2),
                     0.0, 4.90, 0.555,  // amplitude filled in just below, from the data excess
                     hEx->Integral() / std::max(1e-9, hPSfull->Integral()) * 0.2};
   {
      // scale the phase space on 5.6-6.4, where it describes the data, then take what the data has
      // ABOVE it on 4.3-5.6 as the seed for the broad Gaussian
      double sTail = hEx->Integral(hEx->FindBin(5.6), hEx->FindBin(6.4)) /
                     std::max(1e-9, hPSfull->Integral(hPSfull->FindBin(5.6), hPSfull->FindBin(6.4)));
      double dOver = hEx->Integral(hEx->FindBin(4.3), hEx->FindBin(5.6)) -
                     sTail * hPSfull->Integral(hPSfull->FindBin(4.3), hPSfull->FindBin(5.6));
      ini[12] = sTail;
      ini[9] = std::max(1.0, dOver) * hEx->GetBinWidth(1) / (ini[11] * std::sqrt(2 * TMath::Pi()));
      printf("  seed: PS scale %.4g from the 5.6-6.4 tail; excess on 4.3-5.6 = %.0f counts -> A_4.9 = %.1f\n",
             sTail, dOver, ini[9]);
   }
   for (int i = 0; i < 13; ++i) F->SetParameter(i, ini[i]);
   F->SetParLimits(0, 0, 1e7); F->SetParLimits(1, -0.35, 0.35); F->SetParLimits(2, 0.08, 0.35);
   F->SetParLimits(3, 0, 1e7); F->SetParLimits(4, 0.45, 1.15); F->SetParLimits(5, 0.08, 0.35);
   F->SetParLimits(6, 0, 1e7); F->SetParLimits(7, 3.0, 3.8);
   // WIDTHS ARE RESOLUTION AND CANNOT DEPEND ON THE BACKGROUND MODEL. Fitted WITHOUT the phase
   // space (4 gaussians + pol1) these came out 0.200 and 0.555. Left loose against a free phase
   // space the 3.4 peak collapsed to 0.109 while the broad one grew to 0.946 -- the two were
   // trading counts, and 0.109 is narrower than the bound-state resolution at the same angles,
   // which is impossible for a state at 3.4 MeV. Hold both near the no-PS values.
   F->SetParLimits(8, 0.17, 0.26);
   F->SetParLimits(9, 0, 1e7); F->SetParLimits(10, 4.3, 5.5);
   // The ~4.9 structure spans roughly 4.3-5.6, i.e. sigma ~0.35. Allowing sigma down to 0.10 let it
   // collapse to a 14-count spike: it and the FREE phase space are degenerate over that plateau, so
   // the minimizer parks the counts in the phase space and shrinks the Gaussian away. Holding the
   // width broad is what makes the two separable at all.
   F->SetParLimits(11, 0.40, 0.80);
   F->SetParLimits(12, 0, 1e7);
   const char *pn[13] = {"A_gs", "E_gs", "s_gs", "A_074", "E_074", "s_074",
                         "A_bw1", "E_bw1", "G_bw1", "A_bw2", "E_bw2", "G_bw2", "PS"};
   for (int i = 0; i < 13; ++i) F->SetParName(i, pn[i]);
   hEx->Fit(F, "RQ0"); hEx->Fit(F, "RQ0");

   const double bw = hEx->GetBinWidth(1);
   const double sq = std::sqrt(2 * TMath::Pi());
   printf("\n  GLOBAL FIT WITH PHASE SPACE   chi2/ndf = %.2f\n", F->GetChisquare() / std::max(1, F->GetNDF()));
   printf("  g.s.    E = %+.3f +- %.3f   sigma %.3f   area %6.0f +- %.0f\n", F->GetParameter(1), F->GetParError(1),
          F->GetParameter(2), F->GetParameter(0) * F->GetParameter(2) * sq / bw,
          F->GetParError(0) * F->GetParameter(2) * sq / bw);
   printf("  0.740   E = %+.3f +- %.3f   sigma %.3f   area %6.0f +- %.0f\n", F->GetParameter(4), F->GetParError(4),
          F->GetParameter(5), F->GetParameter(3) * F->GetParameter(5) * sq / bw,
          F->GetParError(3) * F->GetParameter(5) * sq / bw);
   printf("  ~3.4    E = %+.3f   sigma %.3f   area %6.0f\n", F->GetParameter(7), F->GetParameter(8),
          F->GetParameter(6) * F->GetParameter(8) * sq / bw);
   printf("  ~4.9    E = %+.3f   sigma %.3f   area %6.0f\n", F->GetParameter(10), F->GetParameter(11),
          F->GetParameter(9) * F->GetParameter(11) * sq / bw);
   printf("  PS scale = %.4g  (phase space integral in range = %.1f counts)\n", F->GetParameter(12),
          F->GetParameter(12) * hPSfull->Integral());
   const double sgG = F->GetParameter(2), sg7 = F->GetParameter(5);

   // ---------------- per-angle fits: means free, sigmas fixed (Fit A) ----------------
   auto *yGS = new TH1D("yGS", "", NB, cmLo, cmHi); yGS->SetDirectory(nullptr);
   auto *y74 = new TH1D("y74", "", NB, cmLo, cmHi); y74->SetDirectory(nullptr);
   auto *sGS = new TH1D("sGS", "", NB, cmLo, cmHi); sGS->SetDirectory(nullptr); // width-treatment syst
   auto *s74 = new TH1D("s74", "", NB, cmLo, cmHi); s74->SetDirectory(nullptr);
   auto *y34 = new TH1D("y34", "", NB, cmLo, cmHi); y34->SetDirectory(nullptr);
   auto *y49 = new TH1D("y49", "", NB, cmLo, cmHi); y49->SetDirectory(nullptr);
   auto *s34 = new TH1D("s34", "", NB, cmLo, cmHi); s34->SetDirectory(nullptr);
   auto *s49 = new TH1D("s49", "", NB, cmLo, cmHi); s49->SetDirectory(nullptr);
   printf("\n  PER-ANGLE (sigmas fixed to the global %.3f / %.3f, means free, PS scale free)\n", sgG, sg7);
   printf("  theta_cm | entries |  E_gs  |  N_gs   stat  syst |  E_074 |  N_074  stat  syst | c2/ndf\n");
   std::vector<TF1 *> FA(NB, nullptr);
   for (int k = 0; k < NB; ++k) {
      double lo = cmLo + k * dcm, hi = lo + dcm;
      if (hSl[k]->Integral() < 40) { printf("  %3.0f-%3.0f  | %7.0f |  too few\n", lo, hi, hSl[k]->Integral()); continue; }
      auto *mk = new PDSpectralModel(pdHistoToGraph(hPS[k]));
      FA[k] = new TF1(Form("FA%d", k), mk, exLo, exHi, 13, "PDSpectralModel");
      // par 5 is tied to par 2 after each fit pass (see below); TF1 has no native parameter tying,
      // so the two are equalised by iterating: fit, copy par2 -> par5 as fixed, refit.
      double sc = hSl[k]->Integral() / std::max(1e-9, hEx->Integral());
      for (int i = 0; i < 13; ++i) FA[k]->SetParameter(i, F->GetParameter(i) * ((i % 3 == 0 || i == 12) ? sc : 1.0));
      // SEED THE BROAD COMPONENT FROM THIS SLICE'S OWN TAIL, exactly as the global fit does. Scaling
      // the global parameters by the slice's statistics is not enough: the broad Gaussian and the
      // free phase space are degenerate over the 4-6 MeV plateau, and starting them anywhere near
      // each other lets the minimizer collapse the Gaussian to zero and leave the phase space
      // running below the data. That is what emptied the ~4.9 component at 60-65, 70-75 and 75-80.
      // THE GAP FIXES THE PHASE-SPACE SCALE. 1.4-2.9 MeV is pure continuum: the bound peaks sit
      // more than 3 sigma below it and the 3.4 peak about 10 sigma above, so nothing else
      // contributes there. Determining the scale from the gap and FIXING it removes the degeneracy
      // between the phase space and the broad Gaussian over 4-6 MeV -- the degeneracy that kept
      // collapsing the broad component at 55-65 and 75-80 and made the continuum over-fill the gap
      // to compensate. Everything above 3 MeV is then left to the peaks, which is where the data
      // actually distinguishes them.
      {
         double pg = hPS[k]->Integral(hPS[k]->FindBin(1.4), hPS[k]->FindBin(2.9));
         double dg = hSl[k]->Integral(hSl[k]->FindBin(1.4), hSl[k]->FindBin(2.9));
         double sGap = (pg > 0 && dg > 3) ? dg / pg : F->GetParameter(12) * sc; // too few -> global
         FA[k]->FixParameter(12, sGap);
         double over = hSl[k]->Integral(hSl[k]->FindBin(4.3), hSl[k]->FindBin(5.6)) -
                       sGap * hPS[k]->Integral(hPS[k]->FindBin(4.3), hPS[k]->FindBin(5.6));
         FA[k]->SetParameter(9, std::max(0.5, over) * hSl[k]->GetBinWidth(1) /
                                   (F->GetParameter(11) * std::sqrt(2 * TMath::Pi())));
         printf("      PS scale from the gap: %ld data / %.1f PS = %.4g%s\n", (long)dg, pg, sGap,
                (pg > 0 && dg > 3) ? "" : "  (too few counts -> global)");
      }
      // Per slice only the BOUND states are wanted, and a few hundred counts cannot constrain the
      // unbound structures as well: their means, widths and the phase-space shape are degenerate
      // over a region with almost no statistics, which is what produced chi2/ndf 11 and errors of
      // 1e7. Freeze the unbound POSITIONS and widths at the global values and let only their
      // amplitudes float, so they can still absorb the counts above 3 MeV without wandering.
      // WIDTHS FROZEN for the central value. Freeing one shared width fits better (chi2 0.70-1.74
      // against 0.79-2.64, and it recovers the real 0.12->0.30 resolution trend) but moves the
      // yields by up to a factor two in the backward bins, where the shared width is what separates
      // two barely-resolved peaks. The frozen fit is the central value; the difference between the
      // two is taken as a systematic below rather than choosing between them.
      FA[k]->FixParameter(2, sgG);  FA[k]->FixParameter(5, sg7);
      FA[k]->FixParameter(8, F->GetParameter(8));  FA[k]->FixParameter(11, F->GetParameter(11));
      // Unbound POSITIONS stay frozen. Giving them a +-0.20/0.35 window was tried and made things
      // worse everywhere (25-30 went 2.64 -> 3.79 with the g.s. sliding to +0.440 and vanishing):
      // the freedom gets spent on the broad-component/continuum degeneracy and the bound peaks pay.
      FA[k]->FixParameter(7, F->GetParameter(7));  FA[k]->FixParameter(10, F->GetParameter(10));
      FA[k]->SetParLimits(0, 0, 1e7); FA[k]->SetParLimits(3, 0, 1e7);
      FA[k]->SetParLimits(6, 0, 1e7); FA[k]->SetParLimits(9, 0, 1e7);
      // NO SetParLimits on par 12 -- it is FixParameter'd from the gap above, and setting limits
      // afterwards silently releases it. That is why the first attempt at this changed nothing.
      FA[k]->SetParLimits(1, -0.45, 0.45); FA[k]->SetParLimits(4, 0.45, 1.15);
      // POISSON LIKELIHOOD, not chi2. Even rebinned, a slice holds a few hundred counts spread
      // over 75 bins, so many bins are single-digit and the Gaussian-error assumption behind chi2
      // fails -- four slices came back with a zero or 1e7 error because Minuit could not build the
      // error matrix. "L" also uses the empty bins, which chi2 discards.
      // POISSON LIKELIHOOD. chi2 discards empty bins, so in a slice with nothing above 4 MeV the
      // broad component is unconstrained and free to invent counts -- at 15-20 it put 26 where the
      // data has 0. And Neyman chi2, weighting by sqrt(N_data), biases low in sparse regions, which
      // is the 1.3-1.7 shortfall the model showed on 4.3-5.6 in EVERY slice from 25 deg up. "L"
      // cures both: it uses the empty bins and is unbiased at low counts.
      hSl[k]->Fit(FA[k], "RQ0L"); hSl[k]->Fit(FA[k], "RQ0L");

      // ---- systematic: repeat with ONE SHARED FREE bound width and keep the difference ----
      auto *mkB = new PDSpectralModel(pdHistoToGraph(hPS[k]));
      TF1 *FB = new TF1(Form("FB%d", k), mkB, exLo, exHi, 13, "PDSpectralModel");
      for (int i = 0; i < 13; ++i) FB->SetParameter(i, FA[k]->GetParameter(i));
      FB->FixParameter(7, F->GetParameter(7));  FB->FixParameter(10, F->GetParameter(10));
      FB->FixParameter(8, F->GetParameter(8));  FB->FixParameter(11, F->GetParameter(11));
      FB->SetParLimits(0, 0, 1e7); FB->SetParLimits(3, 0, 1e7);
      FB->SetParLimits(6, 0, 1e7); FB->SetParLimits(9, 0, 1e7);
      FB->FixParameter(12, FA[k]->GetParameter(12));
      FB->SetParLimits(1, -0.45, 0.45); FB->SetParLimits(4, 0.45, 1.15);
      FB->SetParameter(2, 0.5 * (sgG + sg7)); FB->SetParLimits(2, 0.08, 0.35);
      FB->FixParameter(5, 0.5 * (sgG + sg7));
      hSl[k]->Fit(FB, "RQ0L");
      for (int it = 0; it < 3; ++it) { FB->FixParameter(5, FB->GetParameter(2)); hSl[k]->Fit(FB, "RQ0L"); }
      const double bwS = hSl[k]->GetBinWidth(1); // slice bin width, NOT the inclusive one
      double n0 = FA[k]->GetParameter(0) * sgG * sq / bwS, e0 = FA[k]->GetParError(0) * sgG * sq / bwS;
      double n1 = FA[k]->GetParameter(3) * sg7 * sq / bwS, e1 = FA[k]->GetParError(3) * sg7 * sq / bwS;
      const double swB = FB->GetParameter(2);
      double b0 = FB->GetParameter(0) * swB * sq / bwS, b1 = FB->GetParameter(3) * swB * sq / bwS;
      double s0 = std::fabs(n0 - b0), s1 = std::fabs(n1 - b1); // width-treatment systematic
      // Minuit sometimes returns 0 or ~1e7 for these -- the error matrix failed while the central
      // value converged fine. A zero error is far more dangerous than a large one: it would enter
      // the cross-section as an infinitely precise point and dominate any later fit. Fall back to
      // the Poisson error on the fitted area and say so.
      // POISSON FLOOR, not just a zero-check. Minuit returned errors like 0.4 on a 74-count area
      // (a 2.4 % error where counting statistics alone give 11.6 %) -- not exactly zero, so a
      // zero-test misses them, and they would enter the cross-section as near-infinitely precise
      // points and dominate any subsequent fit. A peak area of N counts cannot be known better
      // than sqrt(N).
      double f0 = std::sqrt(std::max(n0, 1.0)), f1 = std::sqrt(std::max(n1, 1.0));
      bool fb0 = !(e0 > f0 && e0 < 10 * n0), fb1 = !(e1 > f1 && e1 < 10 * n1);
      if (fb0) e0 = f0;
      if (fb1) e1 = f1;
      if (fb0 || fb1) printf("  %3.0f-%3.0f  error below the Poisson floor -> floored: %s%s\n", lo, hi,
                             fb0 ? "g.s. " : "", fb1 ? "0.740" : "");
      yGS->SetBinContent(k + 1, n0); yGS->SetBinError(k + 1, e0);
      y74->SetBinContent(k + 1, n1); y74->SetBinError(k + 1, e1);
      sGS->SetBinContent(k + 1, s0); s74->SetBinContent(k + 1, s1);
      // the two unbound structures. Their Gaussian areas depend on where the continuum is drawn,
      // which the gap-fixed phase space now pins -- but they are still resonances fitted with a
      // Gaussian, so treat them as indicative rather than as yields on the same footing as the
      // bound states.
      {
         double a2 = FA[k]->GetParameter(6) * F->GetParameter(8) * sq / bwS;
         double e2 = FA[k]->GetParError(6) * F->GetParameter(8) * sq / bwS;
         double a3 = FA[k]->GetParameter(9) * F->GetParameter(11) * sq / bwS;
         double e3 = FA[k]->GetParError(9) * F->GetParameter(11) * sq / bwS;
         double b2 = FB->GetParameter(6) * F->GetParameter(8) * sq / bwS;
         double b3 = FB->GetParameter(9) * F->GetParameter(11) * sq / bwS;
         if (!(e2 > std::sqrt(std::max(a2, 1.0)) && e2 < 10 * a2)) e2 = std::sqrt(std::max(a2, 1.0));
         if (!(e3 > std::sqrt(std::max(a3, 1.0)) && e3 < 10 * a3)) e3 = std::sqrt(std::max(a3, 1.0));
         y34->SetBinContent(k + 1, a2); y34->SetBinError(k + 1, e2);
         y49->SetBinContent(k + 1, a3); y49->SetBinError(k + 1, e3);
         s34->SetBinContent(k + 1, std::fabs(a2 - b2)); s49->SetBinContent(k + 1, std::fabs(a3 - b3));
      }
      // per-window closure: what the data has vs what the model puts there. A model/data far from
      // 1 in the continuum windows is a mis-shared broad component, not a statistical fluctuation.
      {
         auto win = [&](double a2, double b2) {
            double d = hSl[k]->Integral(hSl[k]->FindBin(a2), hSl[k]->FindBin(b2));
            double m = FA[k]->Integral(a2, b2) / bwS;
            return std::make_pair(d, m);
         };
         auto w1 = win(1.4, 2.9), w2 = win(4.3, 5.6), w3 = win(5.6, 6.4);
         printf("      windows  1.4-2.9 d/m %4.0f/%-4.0f=%4.2f | 4.3-5.6 %4.0f/%-4.0f=%4.2f | 5.6-6.4 %4.0f/%-4.0f=%4.2f\n",
                w1.first, w1.second, w1.second > 0 ? w1.first / w1.second : 0,
                w2.first, w2.second, w2.second > 0 ? w2.first / w2.second : 0,
                w3.first, w3.second, w3.second > 0 ? w3.first / w3.second : 0);
      }
      printf("  %3.0f-%3.0f  | %7.0f | %+6.3f | %5.0f %5.0f %5.0f | %+6.3f | %5.0f %5.0f %5.0f | %5.2f\n", lo, hi,
             hSl[k]->Integral(), FA[k]->GetParameter(1), n0, e0, s0, FA[k]->GetParameter(4), n1, e1, s1,
             FA[k]->GetChisquare() / std::max(1, FA[k]->GetNDF()));
   }

   // ---------------- acceptance + luminosity ----------------
   const double lumi = lumiFull * (vzHi - vzLo) / zFull;
   printf("\n  luminosity %.1f mb^-1 x %.0f/%.0f mm = %.1f mb^-1\n", lumiFull, vzHi - vzLo, zFull, lumi);
   const char *aTag[4] = {"gs", "ex1", "ex2", "ex3"};
   const char *aNm[4] = {"g.s. (1/2+)", "0.740 (5/2+)", "~3.4 (unbound)", "~4.9 (unbound)"};
   TH1D *Y[4] = {yGS, y74, y34, y49}, *S[4] = {sGS, s74, s34, s49};
   TH1D *X[4] = {nullptr, nullptr, nullptr, nullptr}, *XS[4] = {nullptr, nullptr, nullptr, nullptr};
   auto dOm = [](double a, double b) {
      return 2 * TMath::Pi() * (std::cos(a * TMath::DegToRad()) - std::cos(b * TMath::DegToRad()));
   };
   for (int g = 0; g < 4; ++g) {
      TFile *fa = TFile::Open(here + Form("/../../../../Simulation/ATTPC/16C_pd/diagnostics/acceptance_%s.root", aTag[g]));
      TH1D *A = fa && !fa->IsZombie() ? (TH1D *)fa->Get(Form("hAcc_%s", aTag[g])) : nullptr;
      if (!A) { printf("  no acceptance for %s\n", aTag[g]); if (fa) fa->Close(); continue; }
      A = (TH1D *)A->Clone(Form("A%d", g)); A->SetDirectory(nullptr); fa->Close();
      X[g] = (TH1D *)Y[g]->Clone(Form("X%d", g)); X[g]->Reset(); X[g]->SetDirectory(nullptr);
      XS[g] = (TH1D *)Y[g]->Clone(Form("XS%d", g)); XS[g]->Reset(); XS[g]->SetDirectory(nullptr);
      printf("\n  %s\n  theta_cm |  N fitted    |  acc  |  dsigma/dOmega [mb/sr]  stat   syst\n", aNm[g]);
      for (int b = 1; b <= Y[g]->GetNbinsX(); ++b) {
         double lo = Y[g]->GetBinLowEdge(b), hi = lo + Y[g]->GetBinWidth(b);
         double nn = Y[g]->GetBinContent(b), en = Y[g]->GetBinError(b);
         if (nn <= 0 || en > nn) { if (nn > 0) printf("  %3.0f-%3.0f  | %5.0f +-%4.0f | unusable fit\n", lo, hi, nn, en);
                                   continue; }
         double a = A->GetBinContent(A->FindBin(Y[g]->GetBinCenter(b)));
         double ea = A->GetBinError(A->FindBin(Y[g]->GetBinCenter(b)));
         if (a < 0.15) { printf("  %3.0f-%3.0f  | %5.0f        | %.3f | DROPPED\n", lo, hi, nn, a); continue; }
         double v = nn / dOm(lo, hi) / a / lumi;
         double e = v * std::sqrt(std::pow(en / nn, 2) + std::pow(ea / a, 2));
         double vs = S[g]->GetBinContent(b) / dOm(lo, hi) / a / lumi;
         X[g]->SetBinContent(b, v); X[g]->SetBinError(b, e);
         XS[g]->SetBinContent(b, v); XS[g]->SetBinError(b, vs);
         printf("  %3.0f-%3.0f  | %5.0f +-%4.0f | %.3f |  %8.4g  +-%.3g  +-%.3g\n", lo, hi, nn, en, a, v, e, vs);
      }
   }

   // ---------------- draw ----------------
   gSystem->mkdir(here + "/plots", kTRUE);
   TCanvas *cp = new TCanvas("cp", "ps before/after", 950, 700);
   gPad->SetLogy();
   auto *pr = (TH1F *)hPSraw->Clone("pr"); auto *pw = (TH1F *)hPSfull->Clone("pw");
   pr->SetDirectory(nullptr); pw->SetDirectory(nullptr);
   if (pr->Integral() > 0) pr->Scale(hEx->Integral(hEx->FindBin(4.4), hEx->FindBin(6.4)) /
                                     std::max(1e-9, pr->Integral(pr->FindBin(4.4), pr->FindBin(6.4))));
   if (pw->Integral() > 0) pw->Scale(hEx->Integral(hEx->FindBin(4.4), hEx->FindBin(6.4)) /
                                     std::max(1e-9, pw->Integral(pw->FindBin(4.4), pw->FindBin(6.4))));
   hEx->SetTitle("1n phase space before/after acceptance weighting;E_{x} [MeV];counts");
   hEx->SetLineColor(kBlack); hEx->SetLineWidth(2); hEx->SetMinimum(0.5); hEx->Draw("hist");
   pr->SetLineColor(kRed + 1); pr->SetLineWidth(2); pr->SetLineStyle(2); pr->Draw("hist same");
   pw->SetLineColor(kGreen + 3); pw->SetLineWidth(2); pw->Draw("hist same");
   auto *lp = new TLegend(0.52, 0.70, 0.89, 0.88);
   lp->AddEntry(hEx, "data", "l");
   lp->AddEntry(pr, "PS generator level", "l");
   lp->AddEntry(pw, "PS x acceptance(E_{x},#theta_{cm})", "l");
   lp->AddEntry((TObject *)nullptr, "both normalised on 4.4-6.4 MeV", "");
   lp->Draw();
   cp->SaveAs(here + "/plots/fit_pd_ps_beforeafter_" + tag + ".png");

   TCanvas *cg = new TCanvas("cg", "global", 1000, 700);
   gPad->SetLogy();
   hEx->SetTitle("^{16}C(p,d)^{15}C with 1n phase space;E_{x} [MeV];counts");
   hEx->SetLineWidth(2); hEx->SetMinimum(0.5); hEx->Draw("hist");
   F->SetLineColor(kRed + 1); F->SetLineWidth(2); F->SetNpx(600); F->Draw("same");
   auto *gps = new TGraph();
   for (int i = 1; i <= hPSfull->GetNbinsX(); ++i)
      gps->SetPoint(i - 1, hPSfull->GetBinCenter(i), F->GetParameter(12) * hPSfull->GetBinContent(i));
   gps->SetLineColor(kGreen + 3); gps->SetLineWidth(2); gps->SetLineStyle(2); gps->Draw("L same");
   // each peak on its own, so the decomposition can be judged rather than taken on trust
   const char *cmpNm[4] = {"g.s. (1/2^{+})", "0.740 (5/2^{+})", "~3.4", "~4.9"};
   int cmpCol[4] = {kAzure + 2, kMagenta + 2, kOrange + 7, kCyan + 2};
   TF1 *cmp[4];
   for (int g = 0; g < 4; ++g) {
      cmp[g] = new TF1(Form("cmp%d", g), "gaus", exLo, exHi);
      cmp[g]->SetParameters(F->GetParameter(3 * g), F->GetParameter(3 * g + 1), F->GetParameter(3 * g + 2));
      cmp[g]->SetLineColor(cmpCol[g]); cmp[g]->SetLineWidth(2); cmp[g]->SetLineStyle(7);
      cmp[g]->SetNpx(600); cmp[g]->Draw("same");
   }
   auto *lg = new TLegend(0.58, 0.58, 0.89, 0.88);
   lg->SetFillStyle(0);
   lg->AddEntry(hEx, "data", "l"); lg->AddEntry(F, "total fit", "l");
   for (int g = 0; g < 4; ++g) lg->AddEntry(cmp[g], cmpNm[g], "l");
   lg->AddEntry(gps, "1n phase space", "l"); lg->Draw();
   cg->SaveAs(here + "/plots/fit_pd_ps_global_" + tag + ".png");

   int nc = 4, nr = (NB + nc - 1) / nc;
   TCanvas *cs = new TCanvas("cs", "slices", 340 * nc, 260 * nr);
   cs->Divide(nc, nr);
   for (int k = 0; k < NB; ++k) {
      cs->cd(k + 1); gPad->SetLogy();
      hSl[k]->SetTitle(Form("#theta_{cm} %.0f-%.0f;E_{x} [MeV];counts", cmLo + k * dcm, cmLo + (k + 1) * dcm));
      hSl[k]->SetLineWidth(2); hSl[k]->SetMinimum(0.5); hSl[k]->Draw("hist");
      if (FA[k]) {
         FA[k]->SetLineColor(kRed + 1); FA[k]->SetNpx(600); FA[k]->Draw("same");
         int cc[4] = {kAzure + 2, kMagenta + 2, kOrange + 7, kCyan + 2};
         for (int g = 0; g < 4; ++g) {
            auto *cg2 = new TF1(Form("cs%d_%d", k, g), "gaus", exLo, exHi);
            cg2->SetParameters(FA[k]->GetParameter(3 * g), FA[k]->GetParameter(3 * g + 1),
                               FA[k]->GetParameter(3 * g + 2));
            cg2->SetLineColor(cc[g]); cg2->SetLineWidth(1); cg2->SetLineStyle(7); cg2->SetNpx(400);
            cg2->Draw("same");
         }
         auto *gp2 = new TGraph();
         for (int i = 1; i <= hPS[k]->GetNbinsX(); ++i)
            gp2->SetPoint(i - 1, hPS[k]->GetBinCenter(i), FA[k]->GetParameter(12) * hPS[k]->GetBinContent(i));
         gp2->SetLineColor(kGreen + 3); gp2->SetLineWidth(2); gp2->SetLineStyle(2); gp2->Draw("L same");
      }
   }
   cs->SaveAs(here + "/plots/fit_pd_ps_slices_" + tag + ".png");

   TCanvas *cx = new TCanvas("cx", "xsec", 950, 700);
   gPad->SetLogy(); gPad->SetGridy();
   int col[4] = {kBlack, kRed + 1, kGreen + 3, kOrange + 7}, mk2[4] = {20, 24, 21, 25};
   bool first = true;
   auto *lx = new TLegend(0.60, 0.70, 0.89, 0.88);
   for (int g = 0; g < 4; ++g) {
      if (!X[g]) continue;
      X[g]->SetMarkerStyle(mk2[g]); X[g]->SetMarkerColor(col[g]); X[g]->SetLineColor(col[g]); X[g]->SetLineWidth(2);
      X[g]->SetTitle("^{16}C(p,d)^{15}C;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
      XS[g]->SetFillColorAlpha(col[g], 0.25); XS[g]->SetLineColor(0); XS[g]->SetMarkerSize(0);
      XS[g]->SetTitle("^{16}C(p,d)^{15}C;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
      XS[g]->Draw(first ? "E2" : "E2 same");
      X[g]->Draw("E1 same"); first = false;
      lx->AddEntry(X[g], aNm[g], "lp");
   }
   lx->Draw();
   cx->SaveAs(here + "/plots/xsec_pd_ps_" + tag + ".png");

   TFile fo(here + "/plots/fit_pd_ps_" + tag + ".root", "RECREATE");
   hEx->Write("ex_all"); hPSfull->Write("ps_full");
   yGS->Write("yield_gs"); y74->Write("yield_074");
   for (int g = 0; g < 4; ++g) if (X[g]) { X[g]->Write(Form("dsdo_%s", aTag[g]));
                                            XS[g]->Write(Form("dsdo_%s_syst", aTag[g])); }
   fo.Close();
   printf("\n  wrote plots/fit_pd_ps_{global,slices}_%s.png, xsec_pd_ps_%s.png, fit_pd_ps_%s.root\n\n",
          tag.Data(), tag.Data(), tag.Data());
}
