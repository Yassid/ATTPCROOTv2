/// @file elastic_dip_C14.C
/// @brief Localise the measured elastic diffraction minimum and test it against each OMP.
///
/// The global shape metric in elastic_omp_C14.C is dominated by the dip, where a 2 deg shift in the
/// minimum swings the ratio by an order of magnitude -- so it conflates two questions. This macro
/// separates them: (a) WHERE is the measured minimum, on a fine 2.5 deg grid, and (b) how well does
/// each potential describe the shape AWAY from it, which is what the luminosity actually keys on.
///
///   root -b -q 'elastic_dip_C14.C()'

void elastic_dip_C14(TString cache = "plots/proton_kin_cat5_s013.root",
                     TString accDir = "/mnt/f/a1954_C14_acc_catima_z10_490/", Double_t kSig = 2.5,
                     Double_t dcm = 2.5, Double_t zMin = 10, Double_t zMax = 490, Double_t chi2Cut = 5.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = (TNtuple *)fd->Get("pk");
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   TH1D *acc = (TH1D *)fa->Get("hAcc_gs_sum");
   const double cmMin = 15, cmMax = 150;
   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   std::vector<TH1D *> hs(NB);
   std::vector<double> ctr(NB), mu(NB, 0), w(NB, 0);
   std::vector<bool> ok(NB, false);
   for (int b = 0; b < NB; ++b) { hs[b] = new TH1D(Form("hd%d", b), "", 200, -6, 4); ctr[b] = cmMin + (b + 0.5) * dcm; }
   float *v;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i); v = t->GetArgs();
      if (v[5] > chi2Cut || v[0] <= 0 || v[2] < zMin || v[2] > zMax) continue;
      int b = (int)((v[3] - cmMin) / dcm);
      if (b >= 0 && b < NB) hs[b]->Fill(v[4]);
   }
   for (int b = 0; b < NB; ++b) {
      if (hs[b]->Integral() < 40) continue;
      TH1D *h = (TH1D *)hs[b]->Clone(Form("sd%d", b)); h->Smooth(2);
      double ymax = h->GetMaximum();
      int bm = h->GetMaximumBin(), lo = bm, hi = bm;
      while (lo > 1 && h->GetBinContent(lo) > 0.5 * ymax) --lo;
      while (hi < h->GetNbinsX() && h->GetBinContent(hi) > 0.5 * ymax) ++hi;
      double fwhm = h->GetBinCenter(hi) - h->GetBinCenter(lo);
      if (fwhm <= 0) continue;
      double ww = fwhm / 2.355, m = h->GetBinCenter(bm);
      for (int it = 0; it < 5; ++it) {
         double s = 0, n = 0;
         for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double x = h->GetBinCenter(i);
            if (std::fabs(x - m) > 1.2 * ww) continue;
            s += x * h->GetBinContent(i); n += h->GetBinContent(i);
         }
         if (n > 0) m = s / n;
      }
      mu[b] = m; w[b] = ww; ok[b] = true;
   }
   for (int b = 0; b < NB; ++b) if (!ok[b]) {
      int l = b, r = b;
      while (l >= 0 && !ok[l]) --l;
      while (r < NB && !ok[r]) ++r;
      if (l < 0 && r >= NB) continue;
      if (l < 0) { mu[b] = mu[r]; w[b] = w[r]; }
      else if (r >= NB) { mu[b] = mu[l]; w[b] = w[l]; }
      else { double f = (double)(b - l) / (r - l); mu[b] = mu[l] + f*(mu[r]-mu[l]); w[b] = w[l] + f*(w[r]-w[l]); }
   }
   std::vector<double> ms = mu, ws = w;
   for (int b = 1; b < NB - 1; ++b) { ms[b] = (mu[b-1]+mu[b]+mu[b+1])/3; ws[b] = (w[b-1]+w[b]+w[b+1])/3; }

   std::vector<double> X, Y, E, N;
   for (int b = 0; b < NB; ++b) {
      double lo = ms[b] - kSig*ws[b], hi = ms[b] + kSig*ws[b];
      double y = hs[b]->Integral(hs[b]->FindBin(lo), hs[b]->FindBin(hi));
      double c = ctr[b];
      double dOm = 2*TMath::Pi()*(std::cos((c-dcm/2)*TMath::DegToRad()) - std::cos((c+dcm/2)*TMath::DegToRad()));
      double A = acc->GetBinContent(acc->FindBin(c));
      if (A <= 0.05 || c < 18 || c > 148) continue;
      X.push_back(c); Y.push_back(y / A / dOm); E.push_back(std::sqrt(std::max(y,1.0))/A/dOm); N.push_back(y);
   }
   printf("\n  ==== fine scan through the minimum (%.1f deg bins) ====\n", dcm);
   printf("  theta_cm   counts    dsdo [arb]\n");
   int im = -1; double ymin = 1e99;
   for (size_t j = 0; j < X.size(); ++j) {
      if (X[j] < 45 || X[j] > 75) continue;
      printf("  %7.1f %8.0f %13.4g\n", X[j], N[j], Y[j]);
      if (Y[j] < ymin) { ymin = Y[j]; im = j; }
   }
   // parabolic interpolation in log(sigma) on the three bins about the minimum
   double thDip = X[im];
   if (im > 0 && im + 1 < (int)X.size() && Y[im-1] > 0 && Y[im+1] > 0) {
      double y0 = std::log(Y[im-1]), y1 = std::log(Y[im]), y2 = std::log(Y[im+1]);
      double den = y0 - 2*y1 + y2;
      if (den > 0) thDip = X[im] + 0.5 * dcm * (y0 - y2) / den;
   }
   printf("\n  measured minimum: bin centre %.1f deg, parabolic interpolation %.1f deg\n", X[im], thDip);

   // ---- each potential: its own dip, and its shape away from the dip -------------------------
   const int NP = 5;
   const char *pk[NP] = {"K","V","G","P","M"};
   const char *pn[NP] = {"KD03","CH89","Becchetti-Greenlees","Perey","Menet"};
   TString pdir = here + "/../ptolemy/dat/";
   printf("\n  ==== dip position, and shape away from the dip ====\n");
   printf("  potential                dip[deg]  offset   L(70-148)   rms 20-50 & 70-148\n");
   std::vector<TGraph *> GG; std::vector<double> TD, LL, RR;
   for (int i = 0; i < NP; ++i) {
      auto *g = new TGraph(); std::ifstream in((pdir + TString::Format("el_omp_%s.dat", pk[i])).Data());
      double a, b; while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
      double dmin = 1e99, td = 0;
      for (int j = 0; j < g->GetN(); ++j) { double x = g->GetX()[j];
         if (x < 40 || x > 80) continue; if (g->GetY()[j] < dmin) { dmin = g->GetY()[j]; td = x; } }
      double s = 0; int n = 0;
      for (size_t j = 0; j < X.size(); ++j) {
         if (X[j] < 70 || X[j] > 148 || Y[j] <= 0 || N[j] < 5) continue;
         double f = g->Eval(X[j]); if (f > 0) { s += std::log(Y[j]/f); ++n; }
      }
      double L = std::exp(s/n);
      double r = 0; int m = 0;
      for (size_t j = 0; j < X.size(); ++j) {
         bool in1 = (X[j] >= 20 && X[j] <= 50), in2 = (X[j] >= 70 && X[j] <= 148);
         if ((!in1 && !in2) || Y[j] <= 0 || N[j] < 5) continue;
         double f = L * g->Eval(X[j]); if (f <= 0) continue;
         r += std::pow(std::log(Y[j]/f), 2); ++m;
      }
      printf("  %-24s %7.1f  %+6.1f %10.2f %14.3f\n", pn[i], td, td - thDip, L, std::sqrt(r/m));
      GG.push_back(g); TD.push_back(td); LL.push_back(L); RR.push_back(std::sqrt(r/m));
   }
   // Publish the per-potential luminosity. The consistent optical-model test needs each
   // potential's OWN L, and re-deriving it there would mean a second copy of this extraction --
   // exactly the kind of duplication that lets two numbers drift apart.
   {
      std::ofstream lo((here + "/plots/omp_luminosity.txt").Data());
      lo << "# potential key, name, dip[deg], L[counts/mb], rms away from dip\n";
      lo << "# measured by elastic_dip_C14.C; data dip at " << thDip << " deg\n";
      for (int i = 0; i < NP; ++i)
         lo << pk[i] << " " << pn[i] << " " << TD[i] << " " << LL[i] << " " << RR[i] << "\n";
      printf("  wrote plots/omp_luminosity.txt\n");
   }

   // ---- figures -----------------------------------------------------------------------------
   int col[NP] = {kBlack, kRed+1, kBlue+1, kGreen+2, kMagenta+1};
   std::vector<double> Xg, Yg, Eg;
   for (size_t j = 0; j < X.size(); ++j) if (N[j] >= 3) { Xg.push_back(X[j]); Yg.push_back(Y[j]); Eg.push_back(E[j]); }
   auto *gd = new TGraphErrors(Xg.size(), &Xg[0], &Yg[0], nullptr, &Eg[0]);
   gd->SetMarkerStyle(20); gd->SetMarkerSize(1.1); gd->SetLineWidth(2);

   auto *c1 = new TCanvas("cdip", "", 1400, 1000); c1->Divide(2, 2);
   // (a) zoom on the minimum, every potential on its own L
   c1->cd(1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
   auto *f1 = gPad->DrawFrame(40, 5, 82, 8000);
   f1->SetTitle("the diffraction minimum, potentials on their own L(70-148#circ);#theta_{cm} [deg];d#sigma/d#Omega [arb.]");
   for (int i = 0; i < NP; ++i) { auto *q = new TGraph();
      for (int j = 0; j < GG[i]->GetN(); ++j) q->SetPoint(j, GG[i]->GetX()[j], LL[i]*GG[i]->GetY()[j]);
      q->SetLineColor(col[i]); q->SetLineWidth(3); q->Draw("L same"); }
   gd->Draw("P same");
   auto *lm = new TLine(thDip, 5, thDip, 8000); lm->SetLineColor(kOrange+7);
   lm->SetLineWidth(3); lm->SetLineStyle(2); lm->Draw();
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
   tx.SetTextColor(kOrange+7); tx.DrawLatex(0.15, 0.30, Form("measured minimum %.1f#circ", thDip));
   tx.SetTextColor(kBlack);
   // (b) full range
   c1->cd(2); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
   auto *f2 = gPad->DrawFrame(15, 5, 152, 2e5);
   f2->SetTitle("full range;#theta_{cm} [deg];d#sigma/d#Omega [arb.]");
   auto *lg = new TLegend(0.30, 0.60, 0.93, 0.89);
   lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
   for (int i = 0; i < NP; ++i) { auto *q = new TGraph();
      for (int j = 0; j < GG[i]->GetN(); ++j) q->SetPoint(j, GG[i]->GetX()[j], LL[i]*GG[i]->GetY()[j]);
      q->SetLineColor(col[i]); q->SetLineWidth(2); q->Draw("L same");
      lg->AddEntry(q, Form("%s: dip %+.0f#circ, rms %.3f", pn[i], TD[i]-thDip, RR[i]), "l"); }
   gd->Draw("P same"); lg->Draw();
   // (c) dip offset
   c1->cd(3); gPad->SetGridy();
   auto *hb = new TH1D("hb", "dip position relative to the data;;#theta_{dip}^{OMP} - #theta_{dip}^{data} [deg]", NP, 0, NP);
   for (int i = 0; i < NP; ++i) { hb->SetBinContent(i+1, TD[i]-thDip); hb->GetXaxis()->SetBinLabel(i+1, pn[i]); }
   hb->SetFillColor(kAzure-4); hb->SetBarWidth(0.6); hb->SetBarOffset(0.2);
   hb->GetXaxis()->SetLabelSize(0.055); hb->SetMinimum(-1); hb->SetMaximum(6);
   hb->Draw("bar"); auto *z = new TLine(0,0,NP,0); z->SetLineColor(kRed+1); z->SetLineWidth(2); z->Draw();
   // (d) luminosity
   c1->cd(4); gPad->SetGridy();
   auto *hL = new TH1D("hL", "luminosity from each potential;;L [counts/mb]", NP, 0, NP);
   for (int i = 0; i < NP; ++i) { hL->SetBinContent(i+1, LL[i]); hL->GetXaxis()->SetBinLabel(i+1, pn[i]); }
   hL->SetFillColor(kOrange-3); hL->SetBarWidth(0.6); hL->SetBarOffset(0.2);
   hL->GetXaxis()->SetLabelSize(0.055); hL->SetMinimum(0); hL->SetMaximum(100);
   hL->Draw("bar");
   TLatex t4; t4.SetNDC(); t4.SetTextSize(0.048);
   t4.DrawLatex(0.15, 0.84, Form("spread %.0f%% -- every cross section scales with 1/L",
                                 100*(*std::max_element(LL.begin(),LL.end())/ *std::min_element(LL.begin(),LL.end())-1)));
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "05_elastic_dip_vs_omp.png");
   printf("\n  wrote %s05_elastic_dip_vs_omp.png\n", out.Data());
}
