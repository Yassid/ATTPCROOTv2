/// @file view_pid_peak_vs_integral_h2.C
/// @brief Interactive side-by-side of the two PID variables on the H2 (proton-target) float-tau re-reco (2026-11-11).
///   canvas 1: PEAK amplitude (sqa)   canvas 2: TRACE INTEGRAL (sqw)
///   each: sqrt(dE/dx) vs Brho | ln(sqrt(dE/dx) Brho) vs Brho (bands horizontal) | projection at Brho 0.85-1.0,
///   where the (p,d) deuterons sit (median Brho 0.91)
/// Input: /mnt/f/a1975/reco_h2_tau/pid_width_0*.root (tree w; va = Spyral valid, br = Brho, sqa/sqw).
///   setsid bash -c 'tail -f /dev/null | root -l view_pid_peak_vs_integral_h2.C' &
void view_pid_peak_vs_integral_h2(TString pat = "/mnt/f/a1975/reco_h2_tau/pid_width_0*.root", double bLo = 0.85,
                                  double bHi = 1.0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   auto *ch = new TChain("w");
   int nf = ch->Add(pat);
   printf("%d files, %lld tracks\n", nf, ch->GetEntries());
   const char *var[2] = {"sqa", "sqw"};
   const char *nm[2] = {"PEAK amplitude", "TRACE INTEGRAL"};
   for (int k = 0; k < 2; ++k) {
      auto *c = new TCanvas(Form("c%d", k), Form("H2 PID: %s (%d runs)", nm[k], nf), 50 + 60 * k, 50 + 40 * k, 1500, 520);
      c->Divide(3, 1);
      TString cut = Form("va>0&&br>0.05&&%s>0", var[k]);
      c->cd(1);
      gPad->SetLogz();
      ch->Draw(Form("%s:br>>h2_%d(300,0,2.5,300,0,60)", var[k], k), cut, "colz");
      auto *h2 = (TH2 *)gDirectory->Get(Form("h2_%d", k));
      h2->SetTitle(Form("%s;B#rho [T m];#sqrt{dE/dx} (%s)", nm[k], var[k]));
      c->cd(2);
      gPad->SetLogz();
      ch->Draw(Form("log(%s*br):br>>hl_%d(300,0,2.5,300,0.6,3.4)", var[k], k), cut, "colz");
      auto *hl = (TH2 *)gDirectory->Get(Form("hl_%d", k));
      hl->SetTitle(Form("%s;B#rho [T m];ln(#sqrt{dE/dx} B#rho)", nm[k]));
      auto *l1 = new TLine(bLo, 0.6, bLo, 3.4), *l2 = new TLine(bHi, 0.6, bHi, 3.4);
      for (auto *l : {l1, l2}) {
         l->SetLineColor(kRed);
         l->SetLineStyle(2);
         l->Draw();
      }
      c->cd(3);
      ch->Draw(Form("log(%s*br)>>hp_%d(240,0.6,3.0)", var[k], k), cut + Form("&&br>=%g&&br<%g", bLo, bHi), "hist");
      auto *hp = (TH1 *)gDirectory->Get(Form("hp_%d", k));
      hp->SetTitle(Form("%s, B#rho %.2f-%.2f (p,d deuterons);ln(#sqrt{dE/dx} B#rho);tracks", nm[k], bLo, bHi));
      hp->SetLineColor(k ? kBlue + 1 : kBlack);
      hp->SetLineWidth(2);
      c->Update();
   }
}
