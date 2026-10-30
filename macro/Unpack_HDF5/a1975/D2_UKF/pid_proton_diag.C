/// @file pid_proton_diag.C
/// @brief Why is the 3D proton mass peak wider than Spyral's? Protons = the production (Spyral, amp)
///        ln-mass window. Splits ln m_3D - ln m_Spyral into its Brho part and its dE/dx part and
///        plots each against the 3D choices (vertex end, angle, small-pad path).
///   root -l -b -q 'pid_proton_diag.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016.root","plots/pid_proton_diag.png")'
void pid_proton_diag(TString in, TString png)
{
   gStyle->SetOptStat(0);
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   TString S = "va&&vw&&v3&&sp3<6&&nh>=30";
   auto med = [&](const char *v) {
      w->Draw(v, S, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      return x[x.size() / 2];
   };
   double ma = med("sqa"), m3 = med("sq3a");
   w->SetAlias("lmS", Form("log(sqa/%g*br)", ma));
   w->SetAlias("lm3", Form("log(sq3a/%g*b3)", m3));
   w->SetAlias("dB", "log(b3/br)");
   w->SetAlias("dE", Form("log(sq3a/%g/(sqa/%g))", m3, ma));
   TString P = S + "&&lmS>-2.3&&lmS<-1.75";
   auto rob = [&](const char *v, TString c) {
      w->Draw(v, c, "goff");
      int n = w->GetSelectedRows();
      std::vector<double> x(w->GetV1(), w->GetV1() + n);
      std::sort(x.begin(), x.end());
      return n > 10 ? std::array<double, 3>{x[n * 16 / 100], x[n / 2], x[n * 84 / 100]} : std::array<double, 3>{0, 0, 0};
   };
   printf("protons (production ln-m window): %lld\n", w->GetEntries(P));
   for (auto v : {"lmS", "lm3", "dB", "dE"}) {
      auto r = rob(v, P);
      printf("  %-4s  p16 %7.3f  med %7.3f  p84 %7.3f   half-width %.3f\n", v, r[0], r[1], r[2], 0.5 * (r[2] - r[0]));
   }
   printf("split by 3D end vs Spyral end, and by how the 3D end was chosen (0 sense,1 dca,2 rho):\n");
   for (TString c : {"agreeEnd==1", "agreeEnd==0", "how3==0", "how3==1", "how3==2"}) {
      auto a = rob("dE", P + "&&" + c), b = rob("dB", P + "&&" + c);
      printf("  %-12s n=%5lld   dE half-width %.3f (med %+.3f)   dB half-width %.3f (med %+.3f)\n", c.Data(),
             w->GetEntries(P + "&&" + c), 0.5 * (a[2] - a[0]), a[1], 0.5 * (b[2] - b[0]), b[1]);
   }
   auto *c = new TCanvas("c", "c", 1800, 1000);
   c->Divide(3, 2);
   c->cd(1);
   w->Draw("lm3:lmS>>s1(80,-2.4,-1.6,80,-2.6,-1.4)", P, "colz");
   ((TH2 *)gDirectory->Get("s1"))->SetTitle("protons: ln m, 3D vs Spyral;ln m (Spyral, production);ln m (3D)");
   c->cd(2);
   w->Draw("dE>>e1(80,-0.6,0.6)", P, "hist");
   w->Draw("dB>>b1(80,-0.6,0.6)", P, "hist same");
   auto *e1 = (TH1 *)gDirectory->Get("e1"), *b1 = (TH1 *)gDirectory->Get("b1");
   e1->SetTitle("3D / Spyral, protons;ln(3D / Spyral);tracks");
   e1->SetLineColor(kRed + 1); b1->SetLineColor(kBlue + 1);
   e1->SetMaximum(1.15 * std::max(e1->GetMaximum(), b1->GetMaximum()));
   auto *lg = new TLegend(0.12, 0.75, 0.5, 0.88);
   lg->AddEntry(e1, "#sqrt{dE/dx} part", "l"); lg->AddEntry(b1, "B#rho part", "l"); lg->Draw();
   c->cd(3);
   w->Draw("dE:th3>>t1(60,0,180,60,-0.6,0.6)", P, "colz");
   ((TH2 *)gDirectory->Get("t1"))->SetTitle("dE/dx part vs 3D #theta_{lab};#theta_{lab} 3D [deg];ln(3D/Spyral) #sqrt{dE/dx}");
   c->cd(4);
   w->Draw("dE>>a1(60,-0.6,0.6)", P + "&&agreeEnd==1", "hist");
   w->Draw("dE>>a0(60,-0.6,0.6)", P + "&&agreeEnd==0", "hist same");
   auto *a1 = (TH1 *)gDirectory->Get("a1"), *a0 = (TH1 *)gDirectory->Get("a0");
   a1->SetTitle("dE/dx part: same vertex end as Spyral (black) / opposite (red);ln(3D/Spyral) #sqrt{dE/dx};tracks");
   a0->SetLineColor(kRed + 1);
   c->cd(5);
   w->Draw("dE:pathIn3>>p1(60,0,300,60,-0.6,0.6)", P, "colz");
   ((TH2 *)gDirectory->Get("p1"))->SetTitle("dE/dx part vs 3D small-pad path;3D path on small pads [mm];ln(3D/Spyral) #sqrt{dE/dx}");
   c->cd(6);
   w->Draw("dB:th3>>q1(60,0,180,60,-0.6,0.6)", P, "colz");
   ((TH2 *)gDirectory->Get("q1"))->SetTitle("B#rho part vs 3D #theta_{lab};#theta_{lab} 3D [deg];ln(3D/Spyral) B#rho");
   c->SaveAs(png);
}
