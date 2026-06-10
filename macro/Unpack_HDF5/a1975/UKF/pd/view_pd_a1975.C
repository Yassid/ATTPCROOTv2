/// @file view_pd_a1975.C
/// @brief 16C(p,d)15C viewer — kinematics + (calibrated) excitation energy, in one
///        canvas, fully parameterized so you can zoom / rebin / recalibrate freely.
///
/// Reads the cached per-track ntuple deuteron_kin.root (written by ex_pd_a1975.C;
/// columns ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx). The Ex panel is
/// shifted by `exShift` (default -0.615 = the g.s.->0 calibration) and known 15C
/// levels are drawn for reference.
///
/// Batch (saves PNG):
///   root -b -q 'pd/view_pd_a1975.C'
/// Interactive (mouse-zoom the pads):
///   root -l 'pd/view_pd_a1975.C'
///
/// Args:
///   exShift            Ex calibration offset [MeV] (default -0.615; 0 = raw)
///   thLo,thHi,keLo,keHi,nth,nke   kinematics panel window + binning (KE vs theta_lab)
///   exLo,exHi,nex      Ex panel window + binning
///   cut                optional TTree selection, e.g. "chi2ndf<2" or "theta>15&&theta<35"
///   showLevels         draw known 15C levels on the Ex panel (default true)
///   cacheFile, outPng

void view_pd_a1975(double exShift = -0.615, double thLo = 10, double thHi = 40, double keLo = 0, double keHi = 50,
                   int nth = 300, int nke = 300, double exLo = -3, double exHi = 11, int nex = 140, TString cut = "",
                   bool showLevels = true, TString cacheFile = "deuteron_kin.root", TString outPng = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   if (gSystem->AccessPathName(cacheFile)) {
      printf("ERROR: %s not found (run pd/ex_pd_a1975.C first).\n", cacheFile.Data());
      return;
   }
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   if (!t) {
      printf("ERROR: ntuple 'dk' not in %s\n", cacheFile.Data());
      return;
   }

   // kinematics: KE vs theta_lab
   TH2F *hk = new TH2F("hk", "^{16}C(p,d)^{15}C kinematics;#theta_{lab} [deg];KE_{d} [MeV]", nth, thLo, thHi, nke, keLo,
                       keHi);
   t->Project("hk", "ke:theta", cut);

   // calibrated Ex spectrum: (ex + exShift)
   TH1F *hx = new TH1F("hx", "^{15}C excitation energy (p,d);E_{x}(^{15}C) [MeV];deuteron candidates", nex, exLo, exHi);
   TString exExpr = Form("ex+(%g)", exShift);
   t->Project("hx", exExpr, cut);
   hx->SetLineColor(kBlue + 1);
   hx->SetLineWidth(2);

   long nsel = (long)hk->GetEntries();
   printf("plotted %ld tracks (exShift=%.3f%s)\n", nsel, exShift, cut.Length() ? Form(", cut: %s", cut.Data()) : "");

   TCanvas *c = new TCanvas("cview", "pd view", 1500, 640);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   hk->Draw("colz");
   c->cd(2);
   hx->Draw("hist");
   if (showLevels) {
      double lev[] = {0.0, 0.740, 3.103, 4.220, 4.657, 5.866};
      const char *lbl[] = {"g.s.", "0.74", "3.10", "4.22", "4.66", "5.87"};
      double ymax = hx->GetMaximum();
      for (int i = 0; i < 6; ++i) {
         if (lev[i] < exLo || lev[i] > exHi)
            continue;
         TLine *l = new TLine(lev[i], 0, lev[i], ymax);
         l->SetLineColor(kRed);
         l->SetLineStyle(2);
         l->Draw();
         TLatex *tx = new TLatex(lev[i] + 0.05, ymax * (0.97 - 0.06 * (i % 2)), lbl[i]);
         tx->SetTextColor(kRed + 1);
         tx->SetTextSize(0.03);
         tx->Draw();
      }
   }

   if (outPng.Length() == 0)
      outPng = "pd/plots/view_pd.png";
   c->SaveAs(outPng);
   printf("saved %s\n", outPng.Data());
}
