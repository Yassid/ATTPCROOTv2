/// @file pid_sim_data.C
/// @brief Does a PID estimator behave the SAME in simulation and data? Protons only, species known
///        without any PID cut on both sides:
///          sim  = 16C(d,p) g.s. simulation: every light track (Brho < 1.2 T m) is a proton
///          data = BACKWARD tracks (theta_lab > 100 deg by both geometries): (d,p) protons by kinematics
///        Per estimator and theta bin: ln m relative to that sample's own median (gain-free), half-width
///        and the fraction > 0.3 off. Matching sim/data = an efficiency that transfers through the
///        acceptance; a mismatch is a systematic on the cross section.
///   root -l -b -q 'pid_sim_data.C("/mnt/f/a1975/pid3d_sim/pv_sim_gs.root","/mnt/f/a1975/pid3d_final/pv_all.root","Spyral:sqa:br;3D:s0:b0")'
void pid_sim_data(TString simF, TString dataF, TString vars)
{
   TFile fs(simF), fd(dataF);
   auto *ws = (TTree *)fs.Get("w"), *wd = (TTree *)fd.Get("w");
   TString common = "va&&v3&&sp3<6&&nh>=30&&sq3h>0";
   TString selSim = common + "&&br<1.2&&b3<1.2&&th3>100&&(180-pol)>100";
   TString selDat = common + "&&th3>100&&(180-pol)>100";
   printf("protons: sim %lld, data (backward) %lld\n", ws->GetEntries(selSim), wd->GetEntries(selDat));
   const int bins[4] = {100, 120, 140, 160};
   for (auto *o : *vars.Tokenize(";")) {
      auto *fl = ((TObjString *)o)->GetString().Tokenize(":");
      TString nm = ((TObjString *)fl->At(0))->GetString(), sq = ((TObjString *)fl->At(1))->GetString(),
              br = ((TObjString *)fl->At(2))->GetString();
      TString lm = Form("log((%s)*(%s))", sq.Data(), br.Data());
      auto stats = [&](TTree *w, TString sel, double &med, double &hw, double &far) {
         w->Draw(lm, sel, "goff");
         int n = w->GetSelectedRows();
         std::vector<double> x(w->GetV1(), w->GetV1() + n);
         std::sort(x.begin(), x.end());
         if (n < 10) { med = hw = far = -1; return n; }
         med = x[n / 2];
         hw = 0.5 * (x[n * 84 / 100] - x[n * 16 / 100]);
         int k = 0;
         for (double v : x) k += std::abs(v - med) > 0.3;
         far = 100.0 * k / n;
         return n;
      };
      double mS, hS, fS, mD, hD, fD;
      stats(ws, selSim, mS, hS, fS);
      stats(wd, selDat, mD, hD, fD);
      printf("\n%s   ALL theta>100: sim half-width %.3f far %.1f %% | data half-width %.3f far %.1f %% | ratio data/sim %.2f\n",
             nm.Data(), hS, fS, hD, fD, hD / hS);
      printf("  theta   | sim: n  med-shift  half-width | data: n  med-shift  half-width\n");
      for (int b = 0; b < 3; ++b) {
         TString c = Form("&&th3>=%d&&th3<%d", bins[b], bins[b + 1]);
         double ms, hs, fs2, md, hd, fd2;
         int ns = stats(ws, selSim + c, ms, hs, fs2), nd = stats(wd, selDat + c, md, hd, fd2);
         printf("  %3d-%3d | %6d  %+.3f   %.3f      | %6d  %+.3f   %.3f\n", bins[b], bins[b + 1], ns, ms - mS, hs, nd, md - mD,
                hd);
      }
   }
}
