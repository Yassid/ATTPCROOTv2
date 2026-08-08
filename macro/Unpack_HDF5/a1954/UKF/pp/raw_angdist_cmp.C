/// @file raw_angdist_cmp.C
/// @brief UKF vs GENFIT measured (UNCORRECTED) angular distributions for 14C(p,p').
///
/// No acceptance anywhere in here -- this is the raw yield each fitter delivers, which is what
/// has to be understood before any correction means anything. Four caches are overlaid:
///   UKF chi2/ndf<5, UKF no chi2 cut, GENFIT chi2/ndf<5, GENFIT no chi2 cut.
///
/// Each fitter gets its OWN Ex window, centred on its OWN g.s. peak, because the two disagree by
/// ~0.17 MeV in absolute Ex scale and a common window would be a different physics selection for
/// each. Windows (measured from a gaussian fit to the g.s. of each cache):
///   UKF chi2<5   mu -0.0579   GENFIT chi2<5   mu +0.1164
///   UKF no-chi2  mu -0.0354   GENFIT no-chi2  mu +0.1672
/// elastic  = mu +- 0.6;  inelastic = mu + [5.558, 7.058]  (the offsets the production used).
///
/// Rows: counts per 5 deg bin, then counts/sin(theta_cm) -- an equal-theta bin carries a
/// 2*pi*sin(theta) solid-angle factor, so the second row is the dsigma/dOmega shape (uncorrected).
///
///   root -b -q 'raw_angdist_cmp.C()'

void raw_angdist_cmp(Double_t cmMin = 15.0, Double_t cmMax = 155.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const int NS = 4;
   const char *file[NS] = {"plots/proton_kin_300_ukf.root", "plots/proton_kin_300_ukf_nc.root",
                           "plots/proton_kin_300gfx.root", "plots/proton_kin_300gfx_nc.root"};
   const char *lbl[NS] = {"UKF, #chi^{2}/ndf<5", "UKF, no #chi^{2} cut", "GENFIT, #chi^{2}/ndf<5",
                          "GENFIT, no #chi^{2} cut"};
   const double mu[NS] = {-0.0579, -0.0354, +0.1164, +0.1672};
   const int col[NS] = {kAzure + 2, kAzure + 2, kRed + 1, kRed + 1};
   const int mk[NS] = {20, 24, 21, 25};
   const int lst[NS] = {1, 2, 1, 2};

   const int nb = 28;
   const double lo = 15, hi = 155;

   TH1D *hEl[NS], *hIn[NS], *sEl[NS], *sIn[NS];
   for (int i = 0; i < NS; ++i) {
      TFile *f = TFile::Open(here + "/" + file[i]);
      if (!f || f->IsZombie()) {
         printf("\033[1;31mmissing %s\033[0m\n", file[i]);
         return;
      }
      TTree *t = (TTree *)f->Get("pk");
      hEl[i] = new TH1D(TString::Format("hEl%d", i), "", nb, lo, hi);
      hIn[i] = new TH1D(TString::Format("hIn%d", i), "", nb, lo, hi);
      hEl[i]->Sumw2();
      hIn[i]->Sumw2();
      // TTree::Draw resolves ">>name" through the CURRENT directory, so the histograms must stay
      // attached to the open file until after the Draw; detaching first silently fills nothing.
      t->Draw(TString::Format("thcm>>hEl%d", i), TString::Format("ex>%f&&ex<%f", mu[i] - 0.6, mu[i] + 0.6), "goff");
      t->Draw(TString::Format("thcm>>hIn%d", i), TString::Format("ex>%f&&ex<%f", mu[i] + 5.558, mu[i] + 7.058),
              "goff");
      hEl[i]->SetDirectory(nullptr);
      hIn[i]->SetDirectory(nullptr);
      // strip the solid-angle factor of an equal-theta bin
      sEl[i] = (TH1D *)hEl[i]->Clone(TString::Format("sEl%d", i));
      sIn[i] = (TH1D *)hIn[i]->Clone(TString::Format("sIn%d", i));
      // the clones are born in the open file's directory and would be deleted by Close()
      sEl[i]->SetDirectory(nullptr);
      sIn[i]->SetDirectory(nullptr);
      for (auto *h : {sEl[i], sIn[i]})
         for (int b = 1; b <= nb; ++b) {
            double s = std::sin(h->GetBinCenter(b) * TMath::DegToRad());
            if (s > 1e-3) {
               h->SetBinContent(b, h->GetBinContent(b) / s);
               h->SetBinError(b, h->GetBinError(b) / s);
            }
         }
      f->Close();
   }

   printf("\n===== measured yield per 5 deg bin, NO acceptance correction =====\n");
   printf("            |------------- elastic --------------|----------- Ex 5.5-7.0 -----------|\n");
   printf("  theta_cm  | UKF c2<5  UKF nc  GF c2<5   GF nc  | UKF c2<5  UKF nc  GF c2<5   GF nc |  GF/UKF el (nc)\n");
   double tot[NS] = {0}, totIn[NS] = {0};
   for (int b = 1; b <= nb; ++b) {
      double c = hEl[0]->GetBinCenter(b);
      if (c < cmMin || c > cmMax)
         continue;
      printf("  %3.0f-%3.0f   |", hEl[0]->GetBinLowEdge(b), hEl[0]->GetBinLowEdge(b) + hEl[0]->GetBinWidth(b));
      for (int i = 0; i < NS; ++i)
         printf(" %7.0f", hEl[i]->GetBinContent(b));
      printf("  |");
      for (int i = 0; i < NS; ++i)
         printf(" %7.0f", hIn[i]->GetBinContent(b));
      printf("  |  %6.3f\n", hEl[1]->GetBinContent(b) > 0 ? hEl[3]->GetBinContent(b) / hEl[1]->GetBinContent(b) : 0);
      for (int i = 0; i < NS; ++i) {
         tot[i] += hEl[i]->GetBinContent(b);
         totIn[i] += hIn[i]->GetBinContent(b);
      }
   }
   printf("  TOTAL     |");
   for (int i = 0; i < NS; ++i)
      printf(" %7.0f", tot[i]);
   printf("  |");
   for (int i = 0; i < NS; ++i)
      printf(" %7.0f", totIn[i]);
   printf("  |  %6.3f\n", tot[1] > 0 ? tot[3] / tot[1] : 0);

   TCanvas *c1 = new TCanvas("c1", "raw angular distributions", 1500, 950);
   c1->Divide(2, 2);
   auto style = [&](TH1D *h, int i) {
      h->SetMarkerStyle(mk[i]);
      h->SetMarkerColor(col[i]);
      h->SetLineColor(col[i]);
      h->SetLineStyle(lst[i]);
      h->SetLineWidth(2);
      h->SetMarkerSize(1.1);
      h->GetXaxis()->SetRangeUser(cmMin, cmMax);
   };
   auto panel = [&](int pad, TH1D **h, const char *title, const char *ytitle) {
      c1->cd(pad);
      gPad->SetLogy();
      double ymax = 0;
      for (int i = 0; i < NS; ++i) {
         style(h[i], i);
         ymax = std::max(ymax, h[i]->GetMaximum());
      }
      h[0]->SetTitle(TString::Format("%s;#theta_{cm} [deg];%s", title, ytitle));
      h[0]->SetMinimum(0.5);
      h[0]->SetMaximum(ymax * 5);
      h[0]->Draw("E1");
      for (int i = 1; i < NS; ++i)
         h[i]->Draw("E1 same");
      auto *lg = new TLegend(0.50, 0.68, 0.89, 0.88);
      for (int i = 0; i < NS; ++i)
         lg->AddEntry(h[i], lbl[i], "lp");
      lg->SetTextSize(0.032);
      lg->Draw();
   };
   panel(1, hEl, "elastic window, measured yield", "counts / 5 deg");
   panel(2, hIn, "E_{x} 5.5-7.0 window, measured yield", "counts / 5 deg");
   panel(3, sEl, "elastic, d#sigma/d#Omega shape (uncorrected)", "counts / sin#theta");
   panel(4, sIn, "E_{x} 5.5-7.0, d#sigma/d#Omega shape (uncorrected)", "counts / sin#theta");

   TString png = here + "/plots/raw_angdist_cmp.png";
   c1->SaveAs(png);

   TFile fo(here + "/plots/raw_angdist_cmp.root", "RECREATE");
   const char *nm[NS] = {"ukf_c2", "ukf_nc", "gf_c2", "gf_nc"};
   for (int i = 0; i < NS; ++i) {
      hEl[i]->Write(TString::Format("elastic_%s", nm[i]));
      hIn[i]->Write(TString::Format("inelastic_%s", nm[i]));
   }
   fo.Close();
   printf("\nwrote %s and .root\n\n", png.Data());
}
