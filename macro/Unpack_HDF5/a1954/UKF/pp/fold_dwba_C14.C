/// @file fold_dwba_C14.C
/// @brief Fold the measured theta_cm response into the DWBA before comparing it to the data.
///
/// FRESCO is a zero-resolution calculation. The measurement is that curve convolved with the
/// theta_cm response, and theta_cm is not a measured angle -- the analysis derives it from
/// (KE, theta_lab), so it carries the KE resolution too. The response (built by
/// 14C_pp/response_matrix_C14.C from truth-matched sim) has a sharp core but 0.4-4 % of events
/// beyond +-2 bins. With the elastic cross section spanning ~350x between the forward peak and
/// the 62 deg minimum, that tail is not a detail: it is the dominant contribution to the yield in
/// the minimum, and it flattens the secondary maximum. Comparing raw data to a raw DWBA therefore
/// overstates the disagreement precisely where the sideband extraction overshoots.
///
///     N_true(t)      = dsigma/dOmega(t) * sin(t) * dt
///     N_reco(r)      = SUM_t P(r|t) N_true(t)
///     folded(r)      = N_reco(r) / (sin(r) dt)
///
/// CORRECTED 2026-08-08: an earlier version of this file claimed the sim only populated
/// theta_true above ~15 deg and truncated the fold there. That was WRONG -- acc_batch.sh generates
/// 2-178 deg and the response file has generated AND reconstructed events from 0-5 deg up. The
/// truncation was self-imposed and it suppressed exactly the forward leakage it was meant to
/// bound.
///
/// The migration matrix is normalised per true-column, so P(r|t) is CONDITIONAL on the event being
/// reconstructed. The number that migrates therefore carries the acceptance:
///     N_reco(r) = SUM_t P(r|t) * A(t) * dsigma/dOmega(t) * sin(t) dt
/// and since the data is acceptance-CORRECTED, the prediction to compare against is that divided
/// by A(r) sin(r). Omitting A(t) (as the first version did) over-weights the forward region, where
/// the acceptance falls to nearly zero below 20 deg -- so it would have manufactured feed-down.
/// Both effects are now handled explicitly and the fold runs from foldMin = 0.
///
///   root -b -q 'fold_dwba_C14.C()'

/// normLo defaults to 30, NOT 20: the fold is truncated at foldMin, so nothing migrates INTO the
/// 20-30 deg bins from below and the folded curve is artificially depleted there (fold/raw = 0.79
/// against 0.96-1.02 everywhere else). Anchoring the normalisation on those bins biases every
/// ratio by 7-14 %. Do not lower normLo without also extending the simulation below 15 deg.
void fold_dwba_C14(TString frFile = "",
                   Double_t normLo = 30.0, Double_t normHi = 50.0, Double_t foldMin = 0.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   // FRESCO products live beside this analysis, in ../fresco; resolve them relative to
   // this macro so the comparison is reproducible outside one machine.
   if (frFile.IsNull())
      frFile = here + "/../fresco/outputs/p14C_el_161_dsdo.dat";
   TString simDir = here + "/../../../../Simulation/ATTPC/14C_pp/diagnostics/";

   auto *fr = new TGraph();
   {
      std::ifstream in(frFile.Data());
      double th, xs;
      int n = 0;
      while (in >> th >> xs)
         fr->SetPoint(n++, th, xs);
   }
   if (fr->GetN() == 0) {
      printf("\033[1;31mcannot read %s\033[0m\n", frFile.Data());
      return;
   }

   const int NS = 2;
   const char *rf[NS] = {"response_ukf.root", "response_genfit.root"};
   const char *df[NS] = {"elastic_sideband_gs_ukf.root", "elastic_sideband_gs_gf.root"};
   const char *lbl[NS] = {"UKF", "GENFIT"};
   const int col[NS] = {kAzure + 2, kRed + 1};

   TH1D *dat[NS] = {nullptr, nullptr}, *fold[NS] = {nullptr, nullptr};
   for (int i = 0; i < NS; ++i) {
      TFile *fR = TFile::Open(simDir + rf[i]);
      if (!fR || fR->IsZombie()) {
         printf("\033[1;31mmissing %s -- run response_matrix_C14.C first\033[0m\n", rf[i]);
         return;
      }
      auto *R = (TH2D *)fR->Get("response");
      if (!R)
         return;
      R->SetDirectory(nullptr);
      // the acceptance that goes with this response, needed on BOTH sides of the migration
      TFile *fA = TFile::Open(i == 0 ? "/mnt/f/a1954_C14_acc_nochi2/acceptance_merged_gs.root"
                                     : "/mnt/f/a1954_C14_acc_gf_nochi2/acceptance_merged_gs.root");
      TH1D *Acc = nullptr;
      if (fA && !fA->IsZombie()) {
         Acc = (TH1D *)fA->Get("hAcc_gs_sum");
         if (Acc) Acc = (TH1D *)Acc->Clone(TString::Format("accf%d", i));
         if (Acc) Acc->SetDirectory(nullptr);
         fA->Close();
      }
      if (!Acc) {
         printf("\033[1;31mno acceptance for the fold -- refusing to guess\033[0m\n");
         return;
      }
      fR->Close();

      TFile *fD = TFile::Open(here + "/plots/" + df[i]);
      if (!fD || fD->IsZombie()) {
         printf("\033[1;31mmissing %s -- run elastic_sideband_C14.C first\033[0m\n", df[i]);
         return;
      }
      auto *d = (TH1D *)fD->Get("dsigma_dOmega");
      dat[i] = (TH1D *)d->Clone(TString::Format("dat%d", i));
      dat[i]->SetDirectory(nullptr);
      fD->Close();

      const int nb = R->GetNbinsX();
      auto *f = new TH1D(TString::Format("fold%d", i), "", nb, R->GetXaxis()->GetXmin(), R->GetXaxis()->GetXmax());
      f->SetDirectory(nullptr);
      // N_true then migrate
      std::vector<double> Nt(nb + 1, 0.0);
      for (int t = 1; t <= nb; ++t) {
         double th = R->GetXaxis()->GetBinCenter(t);
         if (th < foldMin)
            continue;
         double s = std::sin(th * TMath::DegToRad());
         double a = Acc->GetBinContent(Acc->FindBin(th)); // events only migrate if reconstructed
         Nt[t] = fr->Eval(th) * s * a;
      }
      for (int r = 1; r <= nb; ++r) {
         double acc = 0;
         for (int t = 1; t <= nb; ++t)
            acc += R->GetBinContent(t, r) * Nt[t];
         double sr = std::sin(f->GetBinCenter(r) * TMath::DegToRad());
         double ar = Acc->GetBinContent(Acc->FindBin(f->GetBinCenter(r)));
         if (sr > 1e-3 && ar > 0.05)
            f->SetBinContent(r, acc / sr / ar); // data is acceptance-corrected, so undo A(r)
      }
      fold[i] = f;
   }

   // normalise data to the FOLDED curve over the forward window
   printf("\n  theta_cm |    raw DWBA   folded  fold/raw |");
   for (int i = 0; i < NS; ++i)
      printf("   %-6s  ratio |", lbl[i]);
   printf("\n");
   double k[NS];
   for (int i = 0; i < NS; ++i) {
      double sn = 0, sd = 0;
      for (int b = 1; b <= dat[i]->GetNbinsX(); ++b) {
         double c = dat[i]->GetBinCenter(b), y = dat[i]->GetBinContent(b), e = dat[i]->GetBinError(b);
         if (c < normLo || c > normHi || y <= 0 || e <= 0)
            continue;
         double ff = fold[i]->GetBinContent(fold[i]->FindBin(c));
         if (ff <= 0)
            continue;
         sn += y * ff / (e * e);
         sd += ff * ff / (e * e);
      }
      k[i] = sn > 0 ? sd / sn : 1.0;
   }
   double chi[NS] = {0, 0};
   int nch[NS] = {0, 0};
   for (int b = 1; b <= fold[0]->GetNbinsX(); ++b) {
      double c = fold[0]->GetBinCenter(b);
      if (c < 18 || c > 148)
         continue;
      double raw = fr->Eval(c), fo = fold[0]->GetBinContent(b);
      if (fo <= 0)
         continue;
      printf("  %3.0f-%3.0f | %10.4g %8.4g %9.2f |", fold[0]->GetBinLowEdge(b), fold[0]->GetBinLowEdge(b) + fold[0]->GetBinWidth(b), raw, fo,
             raw > 0 ? fo / raw : 0);
      for (int i = 0; i < NS; ++i) {
         int bd = dat[i]->FindBin(c);
         double y = dat[i]->GetBinContent(bd) * k[i];
         double fi = fold[i]->GetBinContent(b);
         printf(" %8.4g %6.2f |", y, fi > 0 ? y / fi : 0);
         if (y > 0 && fi > 0) {
            chi[i] += std::pow(std::log(y / fi), 2);
            ++nch[i];
         }
      }
      printf("\n");
   }
   printf("\nrms of ln(data / FOLDED DWBA):");
   for (int i = 0; i < NS; ++i)
      printf("  %s %.3f (%d bins)", lbl[i], nch[i] ? std::sqrt(chi[i] / nch[i]) : 0, nch[i]);
   printf("\n");

   TCanvas *c1 = new TCanvas("c1", "folded", 1300, 560);
   c1->Divide(2, 1);
   c1->cd(1);
   gPad->SetLogy();
   auto *frm = new TH1D("frm", "^{14}C(p,p) elastic vs resolution-folded DWBA;#theta_{cm} [deg];"
                               "d#sigma/d#Omega [mb/sr, shape]",
                        1, 15, 150);
   frm->SetMinimum(0.5);
   frm->SetMaximum(3e3);
   frm->Draw();
   fr->SetLineWidth(2);
   fr->SetLineColor(kGray + 2);
   fr->SetLineStyle(2);
   fr->Draw("L same");
   auto style = [](TH1D *h, int c, int m, int ls = 1) {
      h->SetMarkerStyle(m);
      h->SetMarkerColor(c);
      h->SetLineColor(c);
      h->SetLineWidth(2);
      h->SetLineStyle(ls);
   };
   auto *lg = new TLegend(0.40, 0.64, 0.89, 0.88);
   lg->AddEntry(fr, "DWBA, no resolution", "l");
   for (int i = 0; i < NS; ++i) {
      auto *fs = (TH1D *)fold[i]->Clone(TString::Format("fs%d", i));
      fs->SetLineColor(col[i]);
      fs->SetLineWidth(3);
      fs->SetLineStyle(2);
      fs->Draw("hist same");
      auto *dd = (TH1D *)dat[i]->Clone(TString::Format("dd%d", i));
      dd->Scale(k[i]);
      style(dd, col[i], i == 0 ? 20 : 21);
      dd->Draw("E1 same");
      lg->AddEntry(fs, TString::Format("DWBA folded with %s response", lbl[i]), "l");
      lg->AddEntry(dd, TString::Format("%s data, sideband", lbl[i]), "lp");
   }
   lg->SetTextSize(0.030);
   lg->Draw();

   c1->cd(2);
   gPad->SetLogy();
   auto *r0 = (TH1D *)dat[0]->Clone("r0");
   r0->Reset();
   r0->SetTitle("data / folded DWBA;#theta_{cm} [deg];ratio");
   r0->GetXaxis()->SetRangeUser(15, 150);
   r0->SetMinimum(0.05);
   r0->SetMaximum(20);
   r0->Draw();
   for (int i = 0; i < NS; ++i) {
      auto *r = (TH1D *)dat[i]->Clone(TString::Format("rr%d", i));
      r->Reset();
      for (int b = 1; b <= dat[i]->GetNbinsX(); ++b) {
         double c = dat[i]->GetBinCenter(b), y = dat[i]->GetBinContent(b) * k[i];
         double fi = fold[i]->GetBinContent(fold[i]->FindBin(c));
         if (y <= 0 || fi <= 0)
            continue;
         r->SetBinContent(b, y / fi);
         r->SetBinError(b, dat[i]->GetBinError(b) * k[i] / fi);
      }
      style(r, col[i], i == 0 ? 20 : 21);
      r->Draw("E1 same");
   }
   auto *one = new TLine(15, 1, 150, 1);
   one->SetLineStyle(2);
   one->SetLineColor(kGray + 2);
   one->Draw();
   lg->Draw();

   TString png = here + "/plots/fold_dwba_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
