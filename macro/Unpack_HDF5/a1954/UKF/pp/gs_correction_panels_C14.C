/// @file gs_correction_panels_C14.C
/// @brief The angle correction as three standalone figures: the top row of gs_selection_C14.C.
///
///   1  Ex vs theta_cm BEFORE the theta correction, with the fixed |Ex| < exWin window drawn
///   2  the same AFTER theta' = theta - slope*(KE - pivot), label read from the cache's theta_corr
///   3  the g.s. locus (mode of Ex per 5 deg theta_cm bin) before and after
///
/// Same GENFIT+CATIMA caches and vertex slab as the kinematics plane (kine_plane_genfit_C14.C), so
/// the two figures describe the same tracks. exWin <= 0 draws no window.
///
///   root -l -b -q 'gs_correction_panels_C14.C()'

void gs_correction_panels_C14(TString rawCache = "plots/proton_kin_cat5.root",
                              TString corrCache = "plots/proton_kin_cat5_s013.root",
                              TString cut = "vertexz>10&&vertexz<490", Double_t exWin = 0.6,
                              TString outStem = "plots/gs_correction_C14")
{
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fr_ = TFile::Open(here + "/" + rawCache);
   TFile *fc_ = TFile::Open(here + "/" + corrCache);
   if (!fr_ || fr_->IsZombie() || !fc_ || fc_->IsZombie()) {
      printf("\033[1;31mmissing an input\033[0m\n");
      return;
   }
   TTree *tr = (TTree *)fr_->Get("pk");
   TTree *tc = (TTree *)fc_->Get("pk");
   // the correction applied lives in the corrected cache, not in this macro
   TString corrLabel = "#theta-corrected";
   if (auto *n = (TNamed *)fc_->Get("theta_corr")) {
      corrLabel = n->GetTitle();
      corrLabel.ReplaceAll("(KE-0.00)", " #times K_{p}"); // zero pivot: write it as a product
   }
   printf("raw  %s : %lld of %lld tracks pass \"%s\"\n", rawCache.Data(), tr->GetEntries(cut), tr->GetEntries(),
          cut.Data());
   printf("corr %s : %lld tracks, correction \"%s\"\n", corrCache.Data(), tc->GetEntries(cut), corrLabel.Data());

   auto save = [&](TCanvas *c, const char *suffix) {
      TString png = here + "/" + outStem + suffix + ".png";
      c->SaveAs(png);
      c->SaveAs(TString(png).ReplaceAll(".png", ".pdf"));
   };
   auto label = [](const char *txt) {
      auto *l = new TPaveText(0.42, 0.82, 0.85, 0.89, "NDC");
      l->SetFillColor(kWhite);
      l->SetLineColor(kGray + 1);
      l->SetTextSize(0.040);
      l->AddText(txt);
      l->Draw();
   };

   // ---- 1 & 2: Ex vs theta_cm, raw and corrected
   auto map2 = [&](TTree *t, const char *nm, const char *txt, const char *suffix) {
      auto *c = new TCanvas(TString("c_") + nm, nm, 1200, 900);
      c->SetLogz();
      c->SetRightMargin(0.13);
      auto *h = new TH2D(nm, ";#theta_{cm} [deg];E_{x} [MeV]", 90, 0, 180, 140, -4, 4);
      t->Draw(TString::Format("ex:thcm>>%s", nm), cut, "goff");
      h->SetDirectory(nullptr);
      h->Draw("colz");
      if (exWin > 0)
         for (double s : {-exWin, exWin}) {
            auto *l = new TLine(0, s, 180, s);
            l->SetLineColor(kRed + 1);
            l->SetLineWidth(3);
            l->SetLineStyle(2);
            l->Draw();
         }
      label(txt);
      save(c, suffix);
   };
   map2(tr, "hRaw", "before the angle correction", "_before");
   map2(tc, "hCor", corrLabel, "_after");

   // ---- 3: locus before/after
   auto locus = [&](TTree *t, int col, int mk) {
      auto *g = new TGraph();
      for (double lo = 20; lo < 145; lo += 5) {
         auto *h = new TH1D("hl", "", 160, -4, 3);
         t->Draw("ex>>hl", TString::Format("(%s)&&thcm>=%g&&thcm<%g", cut.Data(), lo, lo + 5), "goff");
         h->SetDirectory(nullptr);
         if (h->Integral() > 60) {
            h->Smooth(1);
            g->SetPoint(g->GetN(), lo + 2.5, h->GetBinCenter(h->GetMaximumBin()));
         }
         delete h;
      }
      g->SetMarkerStyle(mk);
      g->SetMarkerSize(1.4);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->SetLineWidth(3);
      return g;
   };
   auto *c3 = new TCanvas("c_locus", "locus", 1200, 900);
   auto *gR = locus(tr, kGray + 2, 24), *gC = locus(tc, kRed + 1, 20);
   gR->SetTitle(";#theta_{cm} [deg];g.s. E_{x} peak [MeV]");
   gR->GetYaxis()->SetRangeUser(-2.6, 2.2);
   gR->GetXaxis()->SetLimits(15, 150);
   gR->Draw("ALP");
   gC->Draw("LP same");
   for (double s : {-exWin, 0.0, exWin}) {
      if (s != 0 && exWin <= 0)
         continue;
      auto *l = new TLine(15, s, 150, s);
      l->SetLineColor(s == 0 ? kGray + 2 : kRed + 1);
      l->SetLineWidth(2);
      l->SetLineStyle(2);
      l->Draw();
   }
   auto *lg3 = new TLegend(0.42, 0.74, 0.88, 0.88);
   lg3->AddEntry(gR, "before the angle correction", "lp");
   lg3->AddEntry(gC, corrLabel, "lp");
   lg3->SetTextSize(0.036);
   lg3->Draw();
   printf("\n  theta_cm   locus raw   locus corrected\n");
   for (int i = 0; i < gR->GetN(); ++i)
      printf("  %6.1f %11.3f %14.3f\n", gR->GetPointX(i), gR->GetPointY(i), gC->Eval(gR->GetPointX(i)));
   save(c3, "_locus");
   printf("wrote %s{_before,_after,_locus}.{png,pdf}\n", (here + "/" + outStem).Data());
}
