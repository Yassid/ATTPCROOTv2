/// @file pid_width_plot.C
/// @brief Four PID variants from pid_width_dp.C on the SAME tracks: Spyral/3D x peak-amplitude/width-
///        corrected charge. Top: PID planes. Bottom: ln(sqrt(dE/dx)*Brho), proportional to ln(mass)
///        for Z=1, so p, d, t sit ln2 / ln1.5 apart whatever the momentum.
///   root -l -b -q 'pid_width_plot.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016.root","plots/pid_width_0016.png")'
void pid_width_plot(TString in, TString png, TString sel = "va&&vw&&v3&&sp3<6&&nh>=30")
{
   gStyle->SetOptStat(0);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   const char *nm[4] = {"Spyral, peak amplitude (production)", "Spyral, width-corrected", "3D, peak amplitude",
                        "3D, width-corrected"};
   const char *xv[4] = {"sqa", "sqw", "sq3a", "sq3w"}, *yv[4] = {"br", "bw", "b3", "b3"};
   auto *c = new TCanvas("c", "c", 2000, 1000);
   c->Divide(4, 2);
   for (int k = 0; k < 4; ++k) {
      // normalise each sqrt(dE/dx) by its own median so all four share axes
      w->Draw(xv[k], sel, "goff");
      std::vector<double> v(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(v.begin(), v.end());
      double med = v[v.size() / 2];
      c->cd(k + 1);
      gPad->SetLogz();
      w->Draw(Form("%s:%s/%g>>h2_%d(200,0,3,200,0,1.4)", yv[k], xv[k], med, k), sel, "colz");
      ((TH2 *)gDirectory->Get(Form("h2_%d", k)))->SetTitle(Form("%s;#sqrt{dE/dx} / median;B#rho [T m]", nm[k]));
      c->cd(k + 5);
      w->Draw(Form("log(%s/%g*%s)>>h1_%d(200,-3.5,1.5)", xv[k], med, yv[k], k), sel + Form("&&%s>0", yv[k]), "hist");
      ((TH1 *)gDirectory->Get(Form("h1_%d", k)))
         ->SetTitle(Form("%s;ln(#sqrt{dE/dx}#upoint B#rho)  #propto ln m;tracks", nm[k]));
   }
   c->SaveAs(png);
}
