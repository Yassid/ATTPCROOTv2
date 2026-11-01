/// @file plot_sbpid_mass_dp.C
/// @brief Mass distributions per Brho range (2026-11-01): the species-blind fitted mass (pid_sbpid_dp.C, absolute
///        MeV scale, branch `mass`) against the 3D ln-mass variable ln(sqrt(dE/dx)/median * Brho) (s1_g, b1),
///        rescaled so its deuteron peak sits at m_d. Lines at m_p, m_d, m_t. Log y.
///   root -l -b -q 'plot_sbpid_mass_dp.C("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root","emc","plots/sbpid_mass.png")'
void plot_sbpid_mass_dp(TString in, TString mass, TString png, TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&b1>0&&bA>0")
{
   gStyle->SetOptStat(0);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double med = x[x.size() / 2];
   TH1D hd("hd", "", 600, -4, 2);
   w->Draw(Form("log(s1_g/%g*b1)>>hd", med), sel, "goff");
   hd.Smooth(2);
   const double dmode = hd.GetBinCenter(hd.GetMaximumBin());
   const double bl[] = {0.1, 0.3, 0.45, 0.6, 0.8};
   auto *c = new TCanvas("c", "c", 1400, 1000);
   c->Divide(2, 2);
   for (int b = 0; b < 4; ++b) {
      c->cd(b + 1);
      gPad->SetLogy();
      gPad->SetGrid();
      TString cut = sel + Form("&&bA>%g&&bA<%g", bl[b], bl[b + 1]);
      auto *h3 = new TH1D(Form("h3_%d", b), Form("B#rho %.2f-%.2f T m;mass [MeV];tracks", bl[b], bl[b + 1]), 120, 400, 5000);
      auto *hs = new TH1D(Form("hs_%d", b), "", 120, 400, 5000);
      w->Draw(Form("1875.613*exp(log(s1_g/%g*b1)-(%g))>>h3_%d", med, dmode, b), cut, "goff");
      w->Draw(Form("%s>>hs_%d", mass.Data(), b), cut + Form("&&%s>0", mass.Data()), "goff");
      h3->SetLineColor(kGray + 2); h3->SetLineWidth(2);
      hs->SetLineColor(kRed + 1); hs->SetLineWidth(2);
      double ymax = std::max(h3->GetMaximum(), hs->GetMaximum());
      h3->SetMaximum(3 * ymax); h3->SetMinimum(0.8);
      h3->Draw("hist");
      hs->Draw("hist same");
      for (double m : {938.272, 1875.613, 2808.921}) {
         auto *l = new TLine(m, 0.8, m, 3 * ymax);
         l->SetLineStyle(2); l->SetLineColor(kBlue + 1); l->Draw();
      }
      auto *lg = new TLegend(0.55, 0.75, 0.89, 0.89);
      lg->AddEntry(h3, "3D ln-mass (d peak at m_{d})", "l");
      lg->AddEntry(hs, "species-blind fitted mass", "l");
      lg->Draw();
   }
   c->SaveAs(png);
}
