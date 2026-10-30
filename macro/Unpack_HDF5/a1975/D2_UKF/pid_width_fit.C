/// @file pid_width_fit.C
/// @brief Resolution numbers for pid_width_dp.C's four variants: fit ln(sqrt(dE/dx)*Brho) (proportional
///        to ln m for Z=1) with three gaussians (p, d, t) + a constant, on the SAME tracks.
///        Reports each peak's sigma and the p-d and d-t separation in units of the combined sigma.
///   root -l -b -q 'pid_width_fit.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016.root")'
void pid_width_fit(TString in, TString sel = "va&&vw&&v3&&sp3<6&&nh>=30&&sq3c>0", int only = -1)
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   const char *nm[5] = {"Spyral amp (production)", "Spyral width-corr", "3D amp", "3D width-corr", "3D centreline"};
   const char *xv[5] = {"sqa", "sqw", "sq3a", "sq3w", "sq3c"}, *yv[5] = {"br", "bw", "b3", "b3", "b3"};
   printf("%-24s %6s | %7s %7s %7s | %7s %7s %7s | %6s %6s | %s\n", "variant", "N", "mu_p", "mu_d", "mu_t", "sig_p",
          "sig_d", "sig_t", "S(p,d)", "S(d,t)", "chi2/ndf");
   for (int k = 0; k < 5; ++k) {
      if (only >= 0 && k != only && k != 0) continue;
      w->Draw(xv[k], sel, "goff");
      std::vector<double> v(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(v.begin(), v.end());
      double med = v[v.size() / 2];
      TH1D h(Form("h%d", k), "", 100, -3.0, 0.5);
      w->Draw(Form("log(%s/%g*%s)>>h%d", xv[k], med, yv[k], k), sel + Form("&&%s>0", yv[k]), "goff");
      // seed from the deuteron mode; p and t at -ln2 and +ln1.5 from it
      int ib = h.GetMaximumBin();
      double md = h.GetBinCenter(ib), A = h.GetMaximum();
      TF1 g("g", "gaus(0)+gaus(3)+gaus(6)+[9]", -3.0, 0.5);
      g.SetParameters(0.25 * A, md - 0.69, 0.1, A, md, 0.08, 0.15 * A, md + 0.41, 0.12, 1);
      g.SetParLimits(1, md - 0.95, md - 0.45); g.SetParLimits(4, md - 0.15, md + 0.15); g.SetParLimits(7, md + 0.25, md + 0.65);
      for (int s : {2, 5, 8}) g.SetParLimits(s, 0.02, 0.4);
      for (int a : {0, 3, 6}) g.SetParLimits(a, 0, 3 * A);
      g.SetParLimits(9, 0, A);
      h.Fit(&g, "QNRL");
      double mp = g.GetParameter(1), mdd = g.GetParameter(4), mt = g.GetParameter(7);
      double sp = std::abs(g.GetParameter(2)), sd = std::abs(g.GetParameter(5)), st = std::abs(g.GetParameter(8));
      // TAIL: fraction of the tracks in the fit range carried by the flat term (not in any peak)
      double bw = h.GetBinWidth(1), inR = h.Integral();
      double flat = g.GetParameter(9) * h.GetNbinsX();
      printf("%-24s %6.0f | %7.3f %7.3f %7.3f | %7.3f %7.3f %7.3f | %6.2f %6.2f | %.2f | tail %.1f%%\n", nm[k], h.GetEntries(),
             mp, mdd, mt, sp, sd, st, (mdd - mp) / std::hypot(sp, sd), (mt - mdd) / std::hypot(sd, st),
             g.GetChisquare() / std::max(1, g.GetNDF()), 100 * flat / std::max(1.0, inR));
      (void)bw;
      for (int p = 0; p < 10; ++p)
         if (g.GetParError(p) > 0 && (std::abs(g.GetParameter(p) - 0.02) < 1e-4 || std::abs(g.GetParameter(p) - 0.4) < 1e-4))
            printf("    WARNING: parameter %d railed at %.4f\n", p, g.GetParameter(p));
   }
}
