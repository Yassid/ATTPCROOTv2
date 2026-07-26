/// @file explore_C14.C
/// @brief Interactive explorer for the a1954 14C(p,p') UKF protons. Loads the cached
/// kinematics ntuple written by pp/ex_C14.C (e.g. plots/proton_kin_clean161.root) into
/// memory and lets you change the BEAM ENERGY, the cuts and every histogram binning with
/// GUI controls; Ex and theta_cm are RECOMPUTED from (KE, theta_lab) on every redraw.
/// 2D TCanvas only (works over X11 under WSLg; no Eve/OpenGL).
///
/// Run interactively (NOT with -b). Plain ROOT is enough -- no config.sh, no ATTPCROOT
/// libraries, no VMCWORKDIR: the macro only reads floats out of the cache ntuple, and it
/// finds the default cache relative to its own location.
///   root -l /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954/UKF/pp/explore_C14.C
///
/// Another cache / another channel, e.g. 14C(p,d)13C:
///   root -l 'pp/explore_C14.C("plots/proton_kin_pd.root",2.014102,13.003355,"14C(p,d)13C")'
///
/// Controls
///   row 1 : Ebeam, chi2/ndf max, theta_lab window, KE window            (physics cuts)
///   row 2 : Ex bins/lo/hi, g.s. fit window, [Redraw] [Zero g.s.] [Save] (Ex axis + fit)
///   row 3 : 2D binnings, KE max, theta_cm slice for the bottom-right panel
/// [Zero g.s.] solves for the beam energy that puts the elastic peak at Ex = 0.

#include <vector>

/// directory holding this macro, however it was invoked (no VMCWORKDIR needed)
static TString macroDir()
{
   TString self = gInterpreter ? gInterpreter->GetCurrentMacroName() : "";
   if (self.IsNull())
      self = __FILE__;
   TString d = gSystem->DirName(self);
   return d.IsNull() ? TString(".") : d;
}

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics (same expression as ex_C14.C): returns {Ex, theta_cm[deg]}
static std::pair<double, double> kine2b(double m1, double m2, double m3, double m4, double Kp, double thl, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1, uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * omega2(s, m1 * m1, m2 * m2) * omega2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                 (2 * m2 * m2) +
              s + uu - m2 * m2;
   if (a < 0)
      return {NAN, NAN};
   double m4x = std::sqrt(a), ex = m4x - m4;
   double t = m2 * m2 + m4x * m4x - 2 * m2 * Et4;
   double tcm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
                                         (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
                                        (omega2(s, m1 * m1, m2 * m2) * omega2(s, m3 * m3, m4x * m4x)));
   return {ex, tcm * TMath::RadToDeg()};
}

/// ejectile KE vs lab angle locus for p(14C,14C*)p at excitation Ex (CM-angle scan)
static TGraph *kinLine(double Eb, double Ex, double m1, double m2, double m3, double m4_0, Color_t col, int style)
{
   double m4 = m4_0 + Ex;
   double E1 = Eb + m1, P = std::sqrt(E1 * E1 - m1 * m1), W = E1 + m2, s = W * W - P * P;
   auto *g = new TGraph();
   g->SetLineColor(col);
   g->SetLineWidth(2);
   g->SetLineStyle(style);
   if (s <= (m3 + m4) * (m3 + m4))
      return g; // below threshold
   double rs = std::sqrt(s), E3cm = (s + m3 * m3 - m4 * m4) / (2 * rs);
   double p3cm = std::sqrt(std::max(0., E3cm * E3cm - m3 * m3));
   double beta = P / W, gamma = W / rs;
   for (double tc = 0; tc <= 180; tc += 0.5) {
      double c = std::cos(tc * TMath::DegToRad()), sn = std::sin(tc * TMath::DegToRad());
      double Elab = gamma * (E3cm + beta * p3cm * c), pz = gamma * (p3cm * c + beta * E3cm), pperp = p3cm * sn;
      double th = std::atan2(pperp, pz) * TMath::RadToDeg(), ke = Elab - m3;
      if (ke > 0 && th >= 0)
         g->SetPoint(g->GetN(), th, ke);
   }
   return g;
}

class C14Explorer : public TObject {
public:
   C14Explorer(TString cache, double mEjectAmu, double mResidAmu, TString tag, double Ebeam0)
   {
      const double u = 931.49401;
      fMbeam = 14.003242 * u;
      fMtarg = 1.007825 * u;
      fMeject = mEjectAmu * u;
      fMresid = mResidAmu * u;
      fTag = tag;
      fCacheName = gSystem->BaseName(cache);

      TFile *f = TFile::Open(cache);
      if (!f || f->IsZombie()) {
         std::cerr << "\n\033[1;31m[explore_C14] cannot open cache file:\033[0m " << cache << "\n"
                   << "  Build it first with pp/ex_C14.C, or pass the path explicitly:\n"
                   << "    root -l 'explore_C14.C(\"/path/to/proton_kin_clean161.root\")'\n\n";
         return;
      }
      TTree *t = (TTree *)f->Get("pk");
      if (!t) {
         std::cerr << "\n\033[1;31m[explore_C14] no ntuple 'pk' in\033[0m " << cache << "\n\n";
         return;
      }
      float ke = 0, theta = 0, chi2ndf = 0;
      t->SetBranchAddress("ke", &ke);
      t->SetBranchAddress("theta", &theta);
      t->SetBranchAddress("chi2ndf", &chi2ndf);
      Long64_t N = t->GetEntries();
      fKe.reserve(N);
      fTh.reserve(N);
      fC2.reserve(N);
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         fKe.push_back(ke);
         fTh.push_back(theta);
         fC2.push_back(chi2ndf);
      }
      f->Close();
      std::cout << "C14Explorer: loaded " << fKe.size() << " ejectiles from " << cache << "\n";
      MakeGui(Ebeam0);
      Redraw();
   }

   void Redraw()
   {
      const double Eb = fEbeam->GetNumber();
      const int exb = (int)fExBins->GetNumber();
      const double exLo = fExLo->GetNumber(), exHi = fExHi->GetNumber();
      const int thLabB = (int)fThLabBins->GetNumber(), keB = (int)fKEBins->GetNumber();
      const int thCmB = (int)fThCMBins->GetNumber();
      const double keMax = fKEMax->GetNumber();
      const double tcLo = fThCmLo->GetNumber(), tcHi = fThCmHi->GetNumber();

      for (const char *n : {"hEx", "hKT", "hEt", "hExSel"})
         delete gROOT->FindObject(n);
      auto *hEx = new TH1F("hEx", Form("%s  E_{x}   (E_{beam} = %.2f MeV);E_{x} [MeV];counts", fTag.Data(), Eb), exb,
                           exLo, exHi);
      auto *hKT = new TH2F("hKT", "ejectile KE vs #theta_{lab};#theta_{lab} [deg];KE [MeV]", thLabB, 0, 95, keB, 0,
                           keMax);
      auto *hEt = new TH2F("hEt", "E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", thCmB, 0, 180, exb, exLo, exHi);
      auto *hExSel = new TH1F("hExSel", Form("E_{x} for %.0f#circ < #theta_{cm} < %.0f#circ;E_{x} [MeV];counts", tcLo,
                                             tcHi),
                              exb, exLo, exHi);
      hEx->SetLineColor(kBlue + 1);
      hExSel->SetLineColor(kAzure + 2);
      hExSel->SetFillColorAlpha(kAzure + 2, 0.25);

      long n = 0;
      for (size_t i = 0; i < fKe.size(); ++i) {
         if (!Pass(i))
            continue;
         hKT->Fill(fTh[i], fKe[i]);
         auto [ex, thcm] = kine2b(fMbeam, fMtarg, fMeject, fMresid, Eb, fTh[i] * TMath::DegToRad(), fKe[i]);
         if (std::isnan(ex))
            continue;
         ++n;
         hEx->Fill(ex);
         hEt->Fill(thcm, ex);
         if (thcm >= tcLo && thcm <= tcHi)
            hExSel->Fill(ex);
      }

      // g.s. peak fit in the user window
      double mean = 0, sig = 0;
      const double fLo = fFitLo->GetNumber(), fHi = fFitHi->GetNumber();
      if (hEx->GetEntries() > 50) {
         hEx->Fit("gaus", "Q0", "", fLo, fHi);
         if (auto *g = hEx->GetFunction("gaus")) {
            mean = g->GetParameter(1);
            sig = g->GetParameter(2);
         }
      }
      fLabel->SetText(Form("N = %ld    g.s. peak:  mean = %+.3f   sigma = %.3f   FWHM = %.3f MeV", n, mean, sig,
                           2.3548 * sig));

      fCanvas->cd(1);
      gPad->SetLogy(fLogY->IsDown() ? 1 : 0);
      hEx->Draw("hist");
      if (auto *g = hEx->GetFunction("gaus")) {
         g->SetNpx(500);
         g->SetLineColor(kRed + 1);
         g->Draw("same");
      }
      fCanvas->cd(2);
      gPad->SetLogz(1);
      hKT->Draw("colz");
      for (auto *g : fLines)
         delete g;
      fLines.clear();
      if (fKinLines->IsDown()) { // reference loci at the current Ebeam
         struct L {
            double ex;
            Color_t col;
            int st;
         };
         for (L l : {L{0, kRed + 1, 1}, L{6.09, kGreen + 2, 2}, L{6.59, kGray + 2, 2}, L{7.34, kGray + 2, 2}}) {
            TGraph *g = kinLine(Eb, l.ex, fMbeam, fMtarg, fMeject, fMresid, l.col, l.st);
            if (g->GetN() > 0) {
               g->Draw("L same");
               fLines.push_back(g);
            } else
               delete g;
         }
         auto *tx = new TLatex();
         tx->SetNDC();
         tx->SetTextSize(0.045);
         tx->SetTextColor(kRed + 1);
         tx->DrawLatex(0.45, 0.86, "E_{x}=0 (g.s.)");
         tx->SetTextColor(kGreen + 2);
         tx->DrawLatex(0.45, 0.81, "E_{x}=6.09 (1^{-})");
         tx->SetTextColor(kGray + 2);
         tx->DrawLatex(0.45, 0.76, "E_{x}=6.59, 7.34 MeV");
      }
      fCanvas->cd(3);
      gPad->SetLogz(1);
      hEt->Draw("colz");
      fCanvas->cd(4);
      gPad->SetLogy(fLogY->IsDown() ? 1 : 0);
      hExSel->Draw("hist");

      fCanvas->Modified();
      fCanvas->Update();
      gSystem->ProcessEvents();
   }

   /// solve for the beam energy that puts the elastic peak at Ex = 0 (secant iteration)
   void ZeroGs()
   {
      double e1 = fEbeam->GetNumber(), m1 = GsMean(e1);
      double e2 = e1 + 2.0, m2 = GsMean(e2);
      for (int it = 0; it < 12 && std::isfinite(m1) && std::isfinite(m2); ++it) {
         if (std::fabs(m2 - m1) < 1e-9)
            break;
         double e3 = e2 - m2 * (e2 - e1) / (m2 - m1);
         if (!std::isfinite(e3) || e3 < 10 || e3 > 2000)
            break;
         e1 = e2;
         m1 = m2;
         e2 = e3;
         m2 = GsMean(e2);
         printf("  iter %2d : Ebeam = %8.3f MeV -> g.s. mean = %+.4f MeV\n", it, e2, m2);
         if (std::fabs(m2) < 1e-4)
            break;
      }
      printf("[Zero g.s.] Ebeam = %.3f MeV  (%.3f MeV/u)\n", e2, e2 / 14.0);
      fEbeam->SetNumber(e2);
      Redraw();
   }

   void Save()
   {
      TString dir = macroDir() + "/plots";
      if (gSystem->AccessPathName(dir))
         dir = "/tmp";
      TString p = dir + "/explore_C14.png";
      fCanvas->SaveAs(p);
      std::cout << "saved " << p << "\n";
   }

private:
   bool Pass(size_t i) const
   {
      if (fC2[i] > fChi2->GetNumber())
         return false;
      if (fTh[i] < fThLo->GetNumber() || fTh[i] > fThHi->GetNumber())
         return false;
      if (fKe[i] < fKeLo->GetNumber() || fKe[i] > fKeHi->GetNumber())
         return false;
      return true;
   }

   /// g.s. peak centroid at a trial beam energy, with the current cuts (no drawing)
   double GsMean(double Eb)
   {
      TH1F h("hGs", "", (int)fExBins->GetNumber(), fExLo->GetNumber(), fExHi->GetNumber());
      for (size_t i = 0; i < fKe.size(); ++i) {
         if (!Pass(i))
            continue;
         auto [ex, thcm] = kine2b(fMbeam, fMtarg, fMeject, fMresid, Eb, fTh[i] * TMath::DegToRad(), fKe[i]);
         if (!std::isnan(ex))
            h.Fill(ex);
      }
      if (h.GetEntries() < 50)
         return NAN;
      h.Fit("gaus", "Q0", "", fFitLo->GetNumber(), fFitHi->GetNumber());
      auto *g = h.GetFunction("gaus");
      return g ? g->GetParameter(1) : NAN;
   }

   TGNumberEntry *mkNum(TGCompositeFrame *bar, const char *lab, double val, double lo, double hi, int dig = 0)
   {
      bar->AddFrame(new TGLabel(bar, lab), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 2, 3, 3));
      auto *ne = new TGNumberEntry(bar, val, 7, -1, dig ? TGNumberFormat::kNESRealTwo : TGNumberFormat::kNESInteger,
                                   TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax, lo, hi);
      ne->Connect("ValueSet(Long_t)", "C14Explorer", this, "Redraw()"); // arrows redraw immediately
      bar->AddFrame(ne, new TGLayoutHints(kLHintsLeft, 2, 4, 3, 3));
      return ne;
   }

   void MakeGui(double Ebeam0)
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1700, 900);
      main->SetWindowName(fTag + " explorer  [" + fCacheName + "]");

      auto *bar = new TGHorizontalFrame(main);
      fEbeam = mkNum(bar, "Ebeam [MeV]", Ebeam0, 10, 2000, 1);
      fChi2 = mkNum(bar, "chi2/ndf <", 5, 0, 1000, 1);
      fThLo = mkNum(bar, "#theta_{lab} lo", 0, 0, 180, 1);
      fThHi = mkNum(bar, "#theta_{lab} hi", 180, 0, 180, 1);
      fKeLo = mkNum(bar, "KE lo", 0, 0, 500, 1);
      fKeHi = mkNum(bar, "KE hi", 1000, 0, 2000, 1);
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      auto *bar2 = new TGHorizontalFrame(main);
      fExBins = mkNum(bar2, "Ex bins", 200, 10, 5000);
      fExLo = mkNum(bar2, "Ex lo", -5, -50, 0, 1);
      fExHi = mkNum(bar2, "Ex hi", 25, 1, 100, 1);
      fFitLo = mkNum(bar2, "fit lo", -1.0, -20, 20, 1);
      fFitHi = mkNum(bar2, "fit hi", 1.0, -20, 20, 1);
      auto *bR = new TGTextButton(bar2, "  Redraw  ");
      bR->Connect("Clicked()", "C14Explorer", this, "Redraw()");
      bar2->AddFrame(bR, new TGLayoutHints(kLHintsLeft, 12, 4, 3, 3));
      auto *bZ = new TGTextButton(bar2, "  Zero g.s. -> Ebeam  ");
      bZ->Connect("Clicked()", "C14Explorer", this, "ZeroGs()");
      bar2->AddFrame(bZ, new TGLayoutHints(kLHintsLeft, 4, 4, 3, 3));
      auto *bS = new TGTextButton(bar2, "  Save PNG  ");
      bS->Connect("Clicked()", "C14Explorer", this, "Save()");
      bar2->AddFrame(bS, new TGLayoutHints(kLHintsLeft, 4, 4, 3, 3));
      main->AddFrame(bar2, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      auto *bar3 = new TGHorizontalFrame(main);
      fThLabBins = mkNum(bar3, "#theta_{lab} bins", 100, 10, 2000);
      fKEBins = mkNum(bar3, "KE bins", 150, 10, 2000);
      fKEMax = mkNum(bar3, "KE max", 40, 5, 500, 1);
      fThCMBins = mkNum(bar3, "#theta_{cm} bins", 90, 10, 2000);
      fThCmLo = mkNum(bar3, "#theta_{cm} slice lo", 0, 0, 180, 1);
      fThCmHi = mkNum(bar3, "hi", 180, 0, 180, 1);
      fKinLines = new TGCheckButton(bar3, "kin lines");
      fKinLines->SetState(kButtonDown);
      fKinLines->Connect("Clicked()", "C14Explorer", this, "Redraw()");
      bar3->AddFrame(fKinLines, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 4, 3, 3));
      fLogY = new TGCheckButton(bar3, "log y");
      fLogY->Connect("Clicked()", "C14Explorer", this, "Redraw()");
      bar3->AddFrame(fLogY, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 6, 4, 3, 3));
      main->AddFrame(bar3, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      auto *bar4 = new TGHorizontalFrame(main);
      fLabel = new TGLabel(bar4, "                                                                              ");
      bar4->AddFrame(fLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 4, 2, 2));
      main->AddFrame(bar4, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      auto *ec = new TRootEmbeddedCanvas("ec_be12", main, 1680, 780);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();
      fCanvas->Divide(2, 2);
      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }

   std::vector<float> fKe, fTh, fC2;
   double fMbeam{0}, fMtarg{0}, fMeject{0}, fMresid{0};
   TString fTag, fCacheName;
   TCanvas *fCanvas{nullptr};
   TGNumberEntry *fEbeam{nullptr}, *fChi2{nullptr}, *fThLo{nullptr}, *fThHi{nullptr}, *fKeLo{nullptr}, *fKeHi{nullptr},
      *fExBins{nullptr}, *fExLo{nullptr}, *fExHi{nullptr}, *fFitLo{nullptr}, *fFitHi{nullptr}, *fThLabBins{nullptr},
      *fKEBins{nullptr}, *fKEMax{nullptr}, *fThCMBins{nullptr}, *fThCmLo{nullptr}, *fThCmHi{nullptr};
   TGCheckButton *fKinLines{nullptr}, *fLogY{nullptr};
   std::vector<TGraph *> fLines;
   TGLabel *fLabel{nullptr};
   ClassDef(C14Explorer, 0);
};

/// defaults = 14C(p,p') on the a1954 14C+p cache at 161 MeV (11.5 MeV/u)
void explore_C14(TString cache = "", double mEjectAmu = 1.007825, double mResidAmu = 14.003242,
                  TString tag = "14C(p,p')", double Ebeam0 = 161.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   if (gROOT->IsBatch()) {
      std::cerr << "\n\033[1;31m[explore_C14] ROOT is in BATCH mode -- no window can open.\033[0m\n"
                << "  Start it WITHOUT -b:   root -l " << macroDir() << "/explore_C14.C\n\n";
      return;
   }
   if (!gSystem->Getenv("DISPLAY")) {
      std::cerr << "\n\033[1;31m[explore_C14] DISPLAY is not set -- no X server to draw on.\033[0m\n\n";
      return;
   }

   if (cache.IsNull()) { // look next to this macro first, then in the usual workspace spot
      const char *rel[] = {"/plots/proton_kin_clean161.root", "/plots/proton_kin.root"};
      for (const char *r : rel) {
         TString p = macroDir() + r;
         if (!gSystem->AccessPathName(p)) {
            cache = p;
            break;
         }
      }
      if (cache.IsNull() && gSystem->Getenv("VMCWORKDIR"))
         cache = TString(gSystem->Getenv("VMCWORKDIR")) +
                 "/macro/Unpack_HDF5/a1954/UKF/pp/plots/proton_kin_clean161.root";
      if (cache.IsNull()) {
         std::cerr << "\n\033[1;31m[explore_C14] no default cache found next to the macro (" << macroDir()
                   << "/plots/).\033[0m\n  Pass one explicitly: root -l 'explore_C14.C(\"/path/to/cache.root\")'\n\n";
         return;
      }
   }
   new C14Explorer(cache, mEjectAmu, mResidAmu, tag, Ebeam0);
}
