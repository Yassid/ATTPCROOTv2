/// @file dv_compare.C
/// @brief Compare two 16C(d,t)15C caches reconstructed at different drift velocities.
///
/// Both sides must come from the identical PID-gated triton chain, so the only difference is
/// the drift velocity. Reports the two level positions and, as the primary number, how far the
/// whole spectrum moved by cross-correlation -- with a few thousand tritons the tallest-bin
/// position hops between structures and a windowed mean is pinned by its own window once the
/// shift approaches the window width, so neither survives a large step in dv. The correlation
/// uses every bin and degrades gracefully.
///
/// A sigma that lands on its fit limit is flagged rather than reported as a resolution: at
/// dv = 1.25 the 3.103 fit collapsed onto a single-bin fluctuation under every binning tried,
/// which means the level was no longer findable, not that it had become sharp.
///
///   root -b -q 'dv_compare.C("a.root","b.root","1.15","1.50")'

void dv_compare(TString cA, TString cB, TString labA = "A", TString labB = "B", double sigLimit = 0.08)
{
   gStyle->SetOptStat(0);
   TString sel = "chi2ndf<5&&ic>900&&ic<1300&&ke>5&&vertexz>50&&vertexz<700";
   TH1F *hh[2] = {nullptr, nullptr};
   TString fn[2] = {cA, cB}, nm[2] = {labA, labB};
   double gs[2] = {0, 0};

   printf("\n=== dv_compare: %s vs %s (Ebeam 180, identical gated chain) ===\n", labA.Data(), labB.Data());
   for (int i = 0; i < 2; ++i) {
      TFile *f = TFile::Open(fn[i]);
      if (!f || f->IsZombie()) { printf("cannot open %s\n", fn[i].Data()); return; }
      TTree *t = (TTree *)f->Get("pk");
      if (!t) { printf("no tree pk in %s\n", fn[i].Data()); return; }
      t->Draw(Form("ex>>h%d(90,-3,6)", i), sel, "goff");
      hh[i] = (TH1F *)gDirectory->Get(Form("h%d", i));
      hh[i]->SetDirectory(nullptr);

      TF1 g("g", "gaus(0)+pol1(3)", 1.9, 4.3);
      g.SetParameters(25, 3.1, 0.25, 10, 0);
      g.SetParLimits(1, 2.3, 4.2);
      g.SetParLimits(2, sigLimit, 1.2);
      hh[i]->Fit(&g, "QRN");
      TF1 d("d", "[0]*exp(-0.5*((x-[1])/[3])^2)+[2]*exp(-0.5*((x-[1]-0.740)/[3])^2)+[4]+[5]*x", -1.5, 2.4);
      d.SetParameters(25, 0, 30, 0.25, 10, 0);
      d.SetParLimits(1, -3.0, 5.0);
      d.SetParLimits(3, 0.10, 0.9);
      d.SetParLimits(0, 0, 1e5);
      d.SetParLimits(2, 0, 1e5);
      hh[i]->Fit(&d, "QRN");
      gs[i] = d.GetParameter(1);
      bool pin = std::fabs(std::fabs(g.GetParameter(2)) - sigLimit) < 0.005;
      printf("%-8s N=%5.0f | g.s. %+.3f (sig %.3f) | 3.103 %+.3f (sig %.3f)%s\n", nm[i].Data(), hh[i]->GetEntries(),
             gs[i], std::fabs(d.GetParameter(3)), g.GetParameter(1), std::fabs(g.GetParameter(2)),
             pin ? "   <- sigma pinned: the level is NOT findable here" : "");
   }

   double bs = 0, bc = -1, bw = hh[0]->GetBinWidth(1);
   for (int sh = -60; sh <= 60; ++sh) {
      double c = 0, na = 0, nb = 0;
      for (int b = 1; b <= hh[0]->GetNbinsX(); ++b) {
         int j = b + sh;
         if (j < 1 || j > hh[0]->GetNbinsX()) continue;
         c += hh[0]->GetBinContent(b) * hh[1]->GetBinContent(j);
         na += hh[0]->GetBinContent(b) * hh[0]->GetBinContent(b);
         nb += hh[1]->GetBinContent(j) * hh[1]->GetBinContent(j);
      }
      double cc = (na > 0 && nb > 0) ? c / std::sqrt(na * nb) : 0;
      if (cc > bc) { bc = cc; bs = sh * bw; }
   }
   printf("\nspectrum shift %s -> %s : %+.2f MeV   (g.s. moved %+.3f)\n", labA.Data(), labB.Data(), bs, gs[1] - gs[0]);
   double dvA = labA.Atof(), dvB = labB.Atof();
   if (dvB != dvA) {
      double lever = bs / (dvB - dvA);
      printf("lever arm %.1f MeV per cm/us -> g.s. reaches zero at dv = %.3f\n", lever,
             dvA - gs[0] / lever);
   }
}
