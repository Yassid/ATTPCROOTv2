/// @file plot_pid_range_data_dp.C
/// @brief RANGE vs circle Brho on DATA (runs 0016/0026/0058/0088), stopping tracks only (2026-11-01). Left: density plane
///        with the median range curves of the simulated deuterons and tritons overlaid (range_dp.C on the sims).
///        Right: the same plane coloured by the dE/dx class of the circle PID (ln(s1_g * bA) relative to the deuteron
///        mode of each Brho slice): p-like < -0.45 green, deuteron band |.| < 0.15 blue, triton-like 0.25-0.55 red,
///        heavier > 0.55 black -- does the dE/dx triton region also have triton RANGES?
///        Inputs: pid_sbpid_dp.C merged tree + range_dp.C data tree (same rows, friend "r").
///   root -l -b -q 'plot_pid_range_data_dp.C("plots/pid_range_data.png")'
void plot_pid_range_data_dp(TString png, double zFarMax = 955, TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root",
                            TString rr = "/mnt/f/a1975/pid_range/R_all.root", TString R = "/mnt/f/a1975/pid_range/")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   w->AddFriend("r=r", rr);
   TString sel = Form("va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0&&r.ok&&r.cont&&r.zfar>10&&r.zfar<%g&&r.len>20", zFarMax);
   printf("rows pv %lld, range friend %lld; stopping tracks %lld\n", w->GetEntries(), ((TTree *)w->GetListOfFriends()->At(0) ? w->GetFriend("r")->GetEntries() : -1), w->GetEntries(sel));
   // sim median curves
   TGraph *gs[2];
   {
      TChain d("w"), t("w"), dr("r"), trr("r");
      for (int s = 5001; s <= 5004; ++s) { d.Add(Form("/mnt/f/a1975_C16_dd_sim/pv_el_s%d.root", s)); dr.Add(Form("%sr_el_s%d.root", R.Data(), s)); }
      for (const char *s : {"gs", "ex2"}) { t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root", s)); trr.Add(Form("%sr_%s.root", R.Data(), s)); }
      d.AddFriend(&dr, "r"); t.AddFriend(&trr, "r");
      TTree *T[2] = {&d, &t};
      TString ss = Form("va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&bA>0&&r.ok&&r.cont&&r.zfar>10&&r.zfar<%g&&r.len>20", zFarMax);
      for (int s = 0; s < 2; ++s) {
         gs[s] = new TGraph();
         for (double b = 0.10; b < 0.50; b += 0.02) {
            long n = T[s]->Draw("r.len", ss + Form("&&bA>%g&&bA<%g", b, b + 0.02), "goff");
            if (n < 8) continue;
            std::vector<double> v(T[s]->GetV1(), T[s]->GetV1() + n);
            std::sort(v.begin(), v.end());
            gs[s]->SetPoint(gs[s]->GetN(), b + 0.01, v[n / 2]);
         }
         gs[s]->SetLineWidth(3); gs[s]->SetLineColor(s ? kRed + 1 : kBlue + 1); gs[s]->SetLineStyle(s ? 1 : 7);
      }
   }
   auto *c = new TCanvas("c", "c", 1900, 850);
   c->Divide(2, 1);
   c->cd(1); gPad->SetLogz(); gPad->SetLogy(); gPad->SetRightMargin(0.12); gPad->SetGrid();
   const int nbx = 120, nby = 120;
   std::vector<double> ye(nby + 1);
   for (int i = 0; i <= nby; ++i) ye[i] = 20 * std::pow(3000. / 20, double(i) / nby);
   auto *h = new TH2D("h", "data, stopping tracks: range vs circle B#rho;B#rho_{circle} [T m];range [mm]", nbx, 0.05, 0.6, nby, ye.data());
   w->Draw("r.len:bA>>h", sel, "goff");
   h->Draw("colz");
   gs[0]->Draw("L"); gs[1]->Draw("L");
   auto *lg = new TLegend(0.14, 0.76, 0.55, 0.88);
   lg->AddEntry(gs[0], "sim deuterons (median)", "l"); lg->AddEntry(gs[1], "sim tritons (median)", "l"); lg->Draw();
   // classes from the circle dE/dx PID, deuteron mode per slice
   c->cd(2); gPad->SetLogy(); gPad->SetGrid();
   auto *fr = gPad->DrawFrame(0.05, 20, 0.6, 3000);
   fr->SetTitle("coloured by the dE/dx class (circle PID);B#rho_{circle} [T m];range [mm]");
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2];
   TString dm = "(";
   for (double b = 0.05; b < 0.6; b += 0.05) {
      TH1D hh("hh", "", 600, -4, 2);
      w->Draw(Form("log(s1_g/%g*bA)>>hh", md), sel + Form("&&bA>%g&&bA<%g", b, b + 0.05), "goff");
      hh.Smooth(2);
      dm += Form("(bA>%g&&bA<=%g)*%g+", b, b + 0.05, hh.GetBinCenter(hh.GetMaximumBin()));
   }
   dm += "0)";
   TString lnm = Form("(log(s1_g/%g*bA)-%s)", md, dm.Data());
   const char *cls[4] = {"<-0.45", ">-0.15&&%s<0.15", ">0.25&&%s<0.55", ">0.55"};
   const char *cn[4] = {"p-like", "deuteron band", "triton-like", "heavier (^{4}He?)"};
   const int cc[4] = {kGreen + 2, kBlue + 1, kRed + 1, kBlack};
   auto *lg2 = new TLegend(0.14, 0.72, 0.5, 0.88);
   for (int k = 0; k < 4; ++k) {
      TString cut = TString(lnm) + Form(cls[k], lnm.Data());
      long n = w->Draw("r.len:bA", sel + "&&" + cut, "goff");
      auto *g = new TGraph(n, w->GetV2(), w->GetV1());
      g->SetMarkerStyle(k == 2 || k == 3 ? 7 : 1); g->SetMarkerColorAlpha(cc[k], k == 1 ? 0.35 : 0.8); g->Draw("P");
      auto *gl = new TGraph(); gl->SetMarkerStyle(20); gl->SetMarkerColor(cc[k]);
      lg2->AddEntry(gl, Form("%s (%ld)", cn[k], n), "p");
      printf("class %-18s %ld\n", cn[k], n);
   }
   gs[0]->Draw("L"); gs[1]->Draw("L");
   lg2->Draw();
   c->SaveAs(png);
}
