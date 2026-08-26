/// @file gs_selection_C14.C
/// @brief The whole g.s. extraction on one page: what was selected, how, and what came out.
///
/// Six panels, left to right, top to bottom:
///   1  Ex vs theta_cm BEFORE the theta correction, with the |Ex|<0.6 window drawn. The elastic
///      locus visibly leaves the window past ~55 deg -- this is why the fixed window failed.
///   2  the same AFTER theta' = theta - slope*(KE - pivot). The locus is straight, so the window
///      holds the peak at every angle.
///   3  the fitted g.s. locus before and after, which is the quantity the slope was tuned on.
///   4  Ex spectra in four theta_cm slices (corrected), window shaded, so the selection is visible
///      rather than asserted.
///   5  raw counts in the window vs theta_cm.
///   6  acceptance-corrected dsigma/dOmega against FRESCO, plus the ratio.
///
///   root -b -q 'gs_selection_C14.C()'

void gs_selection_C14(TString rawCache = "plots/proton_kin_300gfx_nc.root",
                      TString corrCache = "plots/proton_kin_300gfx_nc_tc.root",
                      TString accDir = "/mnt/f/a1954_C14_acc_gf_nochi2/", Double_t exWin = 0.6,
                      TString frFile = "", TString tag = "GENFIT")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (frFile.IsNull())
      frFile = here + "/../fresco/outputs/p14C_el_161_dsdo.dat";

   TFile *fr_ = TFile::Open(here + "/" + rawCache);
   TFile *fc_ = TFile::Open(here + "/" + corrCache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   if (!fr_ || fr_->IsZombie() || !fc_ || fc_->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mmissing an input\033[0m\n");
      return;
   }
   TTree *tr = (TTree *)fr_->Get("pk");
   TTree *tc = (TTree *)fc_->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_gs_sum");
   if (!tr || !tc || !acc)
      return;

   TCanvas *c1 = new TCanvas("c1", "gs selection", 1650, 1000);
   c1->Divide(3, 2);

   // ---- 1 & 2: Ex vs theta_cm, raw and corrected
   auto map2 = [&](TTree *t, const char *nm, const char *ttl, int pad) {
      c1->cd(pad);
      gPad->SetLogz();
      auto *h = new TH2D(nm, TString::Format("%s;#theta_{cm} [deg];E_{x} [MeV]", ttl), 90, 0, 180, 140, -4, 4);
      t->Draw(TString::Format("ex:thcm>>%s", nm), "", "goff");
      h->SetDirectory(nullptr);
      h->Draw("colz");
      for (double s : {-exWin, exWin}) {
         auto *l = new TLine(0, s, 180, s);
         l->SetLineColor(kRed + 1);
         l->SetLineWidth(2);
         l->SetLineStyle(2);
         l->Draw();
      }
      return h;
   };
   map2(tr, "hRaw", "BEFORE: E_{x} vs #theta_{cm} (raw)", 1);
   // The correction actually applied lives in the CACHE, not here. This label was hardcoded at
   // 0.220/MeV(KE-3.5) and so mis-described every cache carrying a different slope -- it read
   // 0.220 while the data had 0.056. Take it from the corrected cache's own title instead.
   TString corrLabel = "AFTER: #theta-corrected";
   {
      TFile *ftmp = TFile::Open(here + "/" + corrCache);
      if (ftmp && !ftmp->IsZombie()) {
         TNamed *n = (TNamed *)ftmp->Get("theta_corr");
         if (n) corrLabel = TString("AFTER: ") + n->GetTitle();
         ftmp->Close();
      }
   }
   map2(tc, "hCor", corrLabel, 2);

   // ---- 3: locus before/after
   c1->cd(3);
   auto locus = [&](TTree *t, int col, int mk) {
      auto *g = new TGraph();
      int n = 0;
      for (double lo = 20; lo < 145; lo += 5) {
         auto *h = new TH1D(TString::Format("hl%d_%d", col, (int)lo), "", 160, -4, 3);
         t->Draw(TString::Format("ex>>hl%d_%d", col, (int)lo), TString::Format("thcm>=%g&&thcm<%g", lo, lo + 5),
                 "goff");
         h->SetDirectory(nullptr);
         if (h->Integral() > 60) {
            h->Smooth(1);
            g->SetPoint(n++, lo + 2.5, h->GetBinCenter(h->GetMaximumBin()));
         }
         delete h;
      }
      g->SetMarkerStyle(mk);
      g->SetMarkerColor(col);
      g->SetLineColor(col);
      g->SetLineWidth(2);
      return g;
   };
   auto *gR = locus(tr, kGray + 2, 24), *gC = locus(tc, kRed + 1, 20);
   gR->SetTitle("g.s. locus (the quantity the slope was tuned on);#theta_{cm} [deg];E_{x} peak [MeV]");
   gR->GetYaxis()->SetRangeUser(-2.6, 2.2);
   gR->GetXaxis()->SetLimits(15, 150);
   gR->Draw("ALP");
   gC->Draw("LP same");
   for (double s : {-exWin, 0.0, exWin}) {
      auto *l = new TLine(15, s, 150, s);
      l->SetLineColor(s == 0 ? kGray + 2 : kRed + 1);
      l->SetLineStyle(2);
      l->Draw();
   }
   auto *lg3 = new TLegend(0.15, 0.72, 0.55, 0.88);
   lg3->AddEntry(gR, "raw", "lp");
   lg3->AddEntry(gC, "#theta-corrected", "lp");
   lg3->Draw();

   // ---- 4: Ex spectra in slices, window shaded
   c1->cd(4);
   const int NS = 4;
   double slo[NS] = {30, 60, 90, 120}, shi[NS] = {40, 70, 100, 130};
   int cols[NS] = {kAzure + 2, kGreen + 3, kOrange + 7, kMagenta + 1};
   double ymax = 0;
   TH1D *hs[NS];
   for (int i = 0; i < NS; ++i) {
      hs[i] = new TH1D(TString::Format("hsl%d", i), "", 120, -3, 3);
      tc->Draw(TString::Format("ex>>hsl%d", i), TString::Format("thcm>=%g&&thcm<%g", slo[i], shi[i]), "goff");
      hs[i]->SetDirectory(nullptr);
      if (hs[i]->Integral() > 0)
         hs[i]->Scale(1.0 / hs[i]->Integral());
      hs[i]->SetLineColor(cols[i]);
      hs[i]->SetLineWidth(2);
      ymax = std::max(ymax, hs[i]->GetMaximum());
   }
   hs[0]->SetTitle("corrected E_{x} in #theta_{cm} slices (area-normalised);E_{x} [MeV];fraction / bin");
   hs[0]->SetMaximum(ymax * 1.25);
   hs[0]->Draw("hist");
   auto *bx = new TBox(-exWin, 0, exWin, ymax * 1.25);
   bx->SetFillColorAlpha(kRed + 1, 0.10);
   bx->SetLineColor(kRed + 1);
   bx->SetLineStyle(2);
   bx->Draw();
   for (int i = 0; i < NS; ++i)
      hs[i]->Draw("hist same");
   auto *lg4 = new TLegend(0.60, 0.62, 0.89, 0.88);
   for (int i = 0; i < NS; ++i)
      lg4->AddEntry(hs[i], TString::Format("#theta_{cm} %.0f-%.0f", slo[i], shi[i]), "l");
   lg4->AddEntry(bx, "selection |E_{x}|<0.6", "f");
   lg4->SetTextSize(0.030);
   lg4->Draw();

   // ---- 5 & 6: yield and dsigma/dOmega
   const int nb = acc->GetNbinsX();
   auto *yld = new TH1D("yld", "", nb, acc->GetXaxis()->GetXmin(), acc->GetXaxis()->GetXmax());
   auto *dsd = (TH1D *)yld->Clone("dsd");
   yld->Sumw2();
   dsd->Sumw2();
   yld->SetDirectory(nullptr);
   dsd->SetDirectory(nullptr);
   printf("\n===== %s, theta-corrected, selection |Ex| < %.2f =====\n", tag.Data(), exWin);
   printf("  theta_cm   raw counts   acceptance   corrected   dsigma/dOmega\n");
   for (int b = 1; b <= nb; ++b) {
      double lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b), ctr = acc->GetBinCenter(b);
      if (ctr < 20 || ctr > 145)
         continue;
      auto *h = new TH1D(TString::Format("hy%d", b), "", 120, -3, 3);
      tc->Draw(TString::Format("ex>>hy%d", b), TString::Format("thcm>=%g&&thcm<%g", lo, lo + wid), "goff");
      h->SetDirectory(nullptr);
      double y = h->Integral(h->FindBin(-exWin), h->FindBin(exWin));
      delete h;
      double A = acc->GetBinContent(b), s = std::sin(ctr * TMath::DegToRad());
      if (y <= 0 || A <= 0.05 || s <= 1e-3)
         continue;
      yld->SetBinContent(b, y);
      yld->SetBinError(b, std::sqrt(y));
      dsd->SetBinContent(b, y / A / s);
      dsd->SetBinError(b, std::sqrt(y) / A / s);
      printf("  %3.0f-%3.0f %10.0f %12.3f %11.0f %14.4g\n", lo, lo + wid, y, A, y / A, y / A / s);
   }

   auto *fr = new TGraph();
   {
      std::ifstream in(frFile.Data());
      double th, xs;
      int n = 0;
      while (in >> th >> xs)
         fr->SetPoint(n++, th, xs);
   }
   double sn = 0, sd = 0;
   for (int b = 1; b <= nb; ++b) {
      double c = dsd->GetBinCenter(b), y = dsd->GetBinContent(b), e = dsd->GetBinError(b);
      if (c < 30 || c > 50 || y <= 0 || e <= 0)
         continue;
      double f = fr->Eval(c);
      sn += y * f / (e * e);
      sd += f * f / (e * e);
   }
   double k = sn > 0 ? sd / sn : 1.0;

   c1->cd(5);
   gPad->SetLogy();
   yld->SetTitle("raw counts inside the window;#theta_{cm} [deg];counts / 5#circ");
   yld->GetXaxis()->SetRangeUser(15, 150);
   yld->SetMarkerStyle(20);
   yld->SetMarkerColor(kAzure + 2);
   yld->SetLineColor(kAzure + 2);
   yld->SetLineWidth(2);
   yld->SetMinimum(0.5);
   yld->Draw("E1");

   c1->cd(6);
   gPad->SetLogy();
   auto *frm = new TH1D("frm", TString::Format("%s g.s.: d#sigma/d#Omega vs FRESCO;#theta_{cm} [deg];mb/sr (shape)",
                                               tag.Data()),
                        1, 15, 150);
   frm->SetMinimum(0.5);
   frm->SetMaximum(3e3);
   frm->Draw();
   fr->SetLineColor(kBlack);
   fr->SetLineWidth(3);
   fr->Draw("L same");
   auto *dn = (TH1D *)dsd->Clone("dn");
   dn->Scale(k);
   dn->SetMarkerStyle(20);
   dn->SetMarkerColor(kRed + 1);
   dn->SetLineColor(kRed + 1);
   dn->SetLineWidth(2);
   dn->Draw("E1 same");
   auto *lg6 = new TLegend(0.42, 0.72, 0.89, 0.88);
   lg6->AddEntry(fr, "FRESCO DWBA (KD03)", "l");
   lg6->AddEntry(dn, TString::Format("%s, #theta-corr, acc-corrected", tag.Data()), "lp");
   lg6->SetTextSize(0.031);
   lg6->Draw();

   double c2 = 0;
   int nc = 0;
   printf("\n  theta_cm   FRESCO    data     ratio\n");
   for (int b = 1; b <= nb; ++b) {
      double c = dsd->GetBinCenter(b);
      if (c < 20 || c > 145 || dsd->GetBinContent(b) <= 0)
         continue;
      double f = fr->Eval(c), d = dsd->GetBinContent(b) * k;
      printf("  %3.0f-%3.0f %9.4g %9.4g %8.2f\n", dsd->GetBinLowEdge(b), dsd->GetBinLowEdge(b) + dsd->GetBinWidth(b),
             f, d, f > 0 ? d / f : 0);
      if (f > 0 && d > 0) {
         c2 += std::pow(std::log(d / f), 2);
         ++nc;
      }
   }
   printf("\n  rms of ln(data/FRESCO) = %.3f over %d bins\n", nc ? std::sqrt(c2 / nc) : 0, nc);

   TString png = here + "/plots/gs_selection_C14.png";
   c1->SaveAs(png);
   TFile fo(here + "/plots/gs_angdist_final.root", "RECREATE");
   yld->Write("yield_raw");
   dsd->Write("dsigma_dOmega");
   fo.Close();
   printf("wrote %s\n         plots/gs_angdist_final.root\n\n", png.Data());
}
