/// @file vertex_cmp_ana.C -- reads vertex_cmp_dp.root. Per 100 mm vz slice, data vs sim: medians and spreads (p16-p84
/// half-width) of the Spyral vertex x, y, r and the genfit vertex x, y; counts. Then the transverse beam profile
/// for the bound sample only, and a PNG of x(vz), y(vz), r distributions and the vz distribution.
double q(std::vector<double> v, double f){ if(v.empty()) return -999; std::sort(v.begin(),v.end()); return v[std::min(v.size()-1,size_t(f*v.size()))]; }
void vertex_cmp_ana(TString extra = "1")
{
   gStyle->SetOptStat(0);
   TFile f("vertex_cmp_dp.root");
   auto t = (TNtuple *)f.Get("v");
   printf("cut: %s\n", extra.Data());
   for (int src = 0; src < 2; ++src) {
      printf("\n=== %s\n  vz       |    n   | Spyral x med (w)   y med (w)   r med  r p90 | genfit x med (w)   y med (w)\n", src ? "SIM g.s." : "DATA gated p");
      for (int z = 0; z < 10; ++z) {
         TString c = Form("src==%d&&svz>=%d&&svz<%d&&(%s)", src, z * 100, z * 100 + 100, extra.Data());
         t->Draw("svx:svy:sqrt(svx*svx+svy*svy):gvx:gvy", c, "goff para");
         int n = t->GetSelectedRows();
         if (n < 5) { printf("  %3d-%3d  | %5d\n", z * 100, z * 100 + 100, n); continue; }
         std::vector<double> v[5];
         for (int j = 0; j < 5; ++j) for (int i = 0; i < n; ++i) if (t->GetVal(j)[i] > -900) v[j].push_back(t->GetVal(j)[i]);
         auto w = [&](int j) { return 0.5 * (q(v[j], 0.84) - q(v[j], 0.16)); };
         printf("  %3d-%3d  | %5d  | %7.2f (%5.2f) %7.2f (%5.2f)  %5.2f  %5.2f | %7.2f (%5.2f) %7.2f (%5.2f)\n", z * 100, z * 100 + 100, n,
                q(v[0], .5), w(0), q(v[1], .5), w(1), q(v[2], .5), q(v[2], .9), q(v[3], .5), w(3), q(v[4], .5), w(4));
      }
   }
   auto c = new TCanvas("c", "c", 1600, 900);
   c->Divide(3, 2);
   const char *vars[4] = {"svx:svz", "svy:svz", "gvx:svz", "gvy:svz"};
   for (int k = 0; k < 4; ++k) {
      c->cd(k + 1);
      for (int s = 0; s < 2; ++s) {
         auto p = new TProfile(Form("p%d%d", k, s), Form("%s (median-ish: profile mean, |x|<40);vz (mm)", vars[k]), 20, 0, 1000, -40, 40);
         t->Project(p->GetName(), vars[k], Form("src==%d&&(%s)", s, extra.Data()));
         p->SetLineColor(s ? kRed : kBlack); p->SetMarkerColor(s ? kRed : kBlack); p->SetMinimum(-15); p->SetMaximum(15);
         p->Draw(s ? "same" : "");
      }
   }
   c->cd(5);
   for (int s = 0; s < 2; ++s) {
      auto h = new TH1D(Form("hr%d", s), "Spyral vertex r;r (mm);normalised", 60, 0, 30);
      t->Project(h->GetName(), "sqrt(svx*svx+svy*svy)", Form("src==%d&&(%s)", s, extra.Data()));
      if (h->Integral() > 0) h->Scale(1 / h->Integral());
      h->SetLineColor(s ? kRed : kBlack); h->Draw(s ? "hist same" : "hist");
   }
   c->cd(6);
   for (int s = 0; s < 2; ++s) {
      auto h = new TH1D(Form("hz%d", s), "Spyral vz;vz (mm);normalised to 200-600", 50, 0, 1000);
      t->Project(h->GetName(), "svz", Form("src==%d&&(%s)", s, extra.Data()));
      double n = h->Integral(h->FindBin(201), h->FindBin(599));
      if (n > 0) h->Scale(1 / n);
      h->SetLineColor(s ? kRed : kBlack); h->Draw(s ? "hist same" : "hist");
   }
   c->SaveAs("vertex_cmp_dp.png");
}
