/// @file eom_selfconsistent_dp.C
/// @brief Self-consistent EOM Brho for the PID (2026-10-31). pid_eom_dp.C fits each track under p, d and t
///        (e<m>b<k>); a wrong hypothesis biases Brho away from the truth (sim: proton under t +7.7 %, triton
///        under p -15 %). Here, per mode m (0 planar, 1 3D), each hypothesis k gives L_k = ln(hs_g/med * e<m>b<k>);
///        the band position mu_k of species k is the mode of L_k near its expected place, and the track takes the
///        hypothesis minimising |L_k - mu_k|. Writes bsc<m> (chosen Brho; the hit-refit hb where no EOM fit),
///        ksc<m> (chosen k, -1 none) next to all input branches.
///   root -l -b -q 'eom_selfconsistent_dp.C("/mnt/f/a1975/pid3d_eom/pv_gainnorm.root","/mnt/f/a1975/pid3d_eom/pv_sc.root")'
void eom_selfconsistent_dp(TString in, TString out, TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&hs_g>0&&hb>0")
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   w->Draw("hs_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double med = x[x.size() / 2];
   // band positions per mode and hypothesis
   double mu[2][3];
   for (int m = 0; m < 2; ++m) {
      TH1D hd("hd", "", 600, -4, 2);
      w->Draw(Form("log(hs_g/%g*e%db1)>>hd", med, m), sel + Form("&&e%db1>0", m), "goff");
      hd.Smooth(2);
      mu[m][1] = hd.GetBinCenter(hd.GetMaximumBin());
      for (int k : {0, 2}) {
         TH1D h("h", "", 600, -4, 2);
         w->Draw(Form("log(hs_g/%g*e%db%d)>>h", med, m, k), sel + Form("&&e%db%d>0", m, k), "goff");
         h.Smooth(2);
         if (k == 0) h.GetXaxis()->SetRangeUser(mu[m][1] - 0.8, mu[m][1] - 0.35);
         else h.GetXaxis()->SetRangeUser(mu[m][1] + 0.22, mu[m][1] + 0.55);
         mu[m][k] = h.GetBinCenter(h.GetMaximumBin());
      }
      printf("mode %d: band positions p %.3f d %.3f t %.3f (t-d %.3f, p-d %.3f)\n", m, mu[m][0], mu[m][1], mu[m][2],
             mu[m][2] - mu[m][1], mu[m][0] - mu[m][1]);
   }
   TFile fo(out, "RECREATE");
   auto *o = w->CloneTree(0);
   double hs_g, hb, eb[2][3], bsc[2];
   int ksc[2];
   w->SetBranchAddress("hs_g", &hs_g);
   w->SetBranchAddress("hb", &hb);
   for (int m = 0; m < 2; ++m)
      for (int k = 0; k < 3; ++k) w->SetBranchAddress(Form("e%db%d", m, k), &eb[m][k]);
   for (int m = 0; m < 2; ++m) { o->Branch(Form("bsc%d", m), &bsc[m]); o->Branch(Form("ksc%d", m), &ksc[m]); }
   long nk[2][4] = {{0}};
   for (Long64_t i = 0; i < w->GetEntries(); ++i) {
      w->GetEntry(i);
      for (int m = 0; m < 2; ++m) {
         bsc[m] = hb; ksc[m] = -1;
         if (hs_g > 0) {
            double best = 1e9;
            for (int k = 0; k < 3; ++k) {
               if (eb[m][k] <= 0) continue;
               double d = std::abs(std::log(hs_g / med * eb[m][k]) - mu[m][k]);
               if (d < best) { best = d; ksc[m] = k; bsc[m] = eb[m][k]; }
            }
         }
         nk[m][ksc[m] + 1]++;
      }
      o->Fill();
   }
   for (int m = 0; m < 2; ++m)
      printf("mode %d: chosen none %ld p %ld d %ld t %ld\n", m, nk[m][0], nk[m][1], nk[m][2], nk[m][3]);
   o->Write();
   fo.Close();
}
