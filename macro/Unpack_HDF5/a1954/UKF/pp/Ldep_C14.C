/// @file Ldep_C14.C
/// @brief Is the luminosity really a constant? The check that decides where to normalise.
///
/// The absolute normalisation takes one number, L = (measured elastic) / (calculated elastic), and
/// applies it to every excited state. L is a beam flux times a target thickness, so it CANNOT
/// depend on scattering angle. Plotting it against angle is therefore a direct test of the
/// normalisation, and it also says where in angle the normalisation may legitimately be taken.
///
/// It is not flat everywhere, and the structure is informative rather than random:
///
///   * 62-72 deg: L spikes to ~110. The calculation has a deep diffraction minimum there
///     (1.5 mb/sr) which the measurement, with finite resolution, fills in. Dividing a shallow
///     measured minimum by a deep calculated one manufactures a spike. Nothing to do with flux.
///   * ~112 deg: L dips to ~9, the same effect inverted, where the data have a local minimum the
///     smooth calculation does not.
///   * 40-60 deg: L is flat at 28-33, because there the measured and calculated SHAPES agree.
///
/// The normalisation is therefore taken over 40-60 deg. Choosing instead a region where the
/// ACCEPTANCE is flat (75-120 deg) is the wrong criterion and was the original mistake here: that
/// range sits on the secondary maximum, where the measured peak (~82 deg) and the calculated one
/// (~93 deg) are displaced, so the ratio slides monotonically through it and the "systematic"
/// came out as -73%/+134% instead of -13%/+1%.
///
///   root -b -q 'Ldep_C14.C()'

void Ldep_C14(TString gsFile = "plots/elastic_sideband_gs_gf.root", TString gsCache = "plots/proton_kin_300gfx_nc.root",
              TString accDir = "/mnt/f/a1954_C14_acc_gf_z10_400/", TString frDir = "", Double_t zMin = 10.0,
              Double_t zMax = 400.0, Double_t gsLo = -0.5, Double_t gsHi = 0.7, Double_t Ladopted = 32.31,
              Double_t normLo = 40.0, Double_t normHi = 60.0, TString tag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (frDir.IsNull())
      frDir = here + "/../fresco/outputs/";

   TFile *fy = TFile::Open(here + "/" + gsFile);
   TFile *fc = TFile::Open(here + "/" + gsCache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   if (!fy || fy->IsZombie() || !fc || fc->IsZombie() || !fa || fa->IsZombie()) {
      printf("\033[1;31mmissing the g.s. yield, the cache or the acceptance\033[0m\n");
      return;
   }
   auto *Y = (TH1D *)fy->Get("yield_sideband");
   auto *t = (TTree *)fc->Get("pk");
   auto *A = (TH1D *)fa->Get("hAcc_gs_sum");
   if (!Y || !t || !A)
      return;

   TGraph fr;
   {
      std::ifstream in((frDir + "p14C_el_161_dsdo.dat").Data());
      double a, x;
      int n = 0;
      while (in >> a >> x)
         fr.SetPoint(n++, a, x);
      if (fr.GetN() == 0) {
         printf("\033[1;31mcannot read the KD03 elastic\033[0m\n");
         return;
      }
   }
   auto dOmega = [](double a, double b) {
      return 2 * TMath::Pi() * (std::cos(a * TMath::DegToRad()) - std::cos(b * TMath::DegToRad()));
   };

   auto *gL = new TGraphErrors();
   auto *gD = new TGraph();
   auto *gT = new TGraph();
   int n = 0, m = 0;
   std::vector<double> flat;
   printf("\n  theta_cm |    L    | L / %.1f\n", Ladopted);
   for (int b = 1; b <= Y->GetNbinsX(); ++b) {
      double lo = Y->GetBinLowEdge(b), w = Y->GetBinWidth(b), c = Y->GetBinCenter(b);
      double y = Y->GetBinContent(b), ey = Y->GetBinError(b);
      if (y <= 0)
         continue;
      TString base = TString::Format("ex>%g&&ex<%g&&thcm>=%g&&thcm<%g", gsLo, gsHi, lo, lo + w);
      double nAll = t->GetEntries(base);
      double nIn = t->GetEntries(base + TString::Format("&&vertexz>%g&&vertexz<%g", zMin, zMax));
      if (nAll < 10)
         continue;
      double f = nIn / nAll, a = A->GetBinContent(A->FindBin(c)), s = fr.Eval(c);
      if (a <= 0.05 || s <= 0)
         continue;
      double dat = y * f / a / dOmega(lo, lo + w);
      gL->SetPoint(n, c, dat / s);
      gL->SetPointError(n, 0, ey * f / a / dOmega(lo, lo + w) / s);
      ++n;
      gD->SetPoint(m, c, dat);
      gT->SetPoint(m, c, s);
      ++m;
      printf("  %6.0f   | %7.1f | %6.2f%s\n", c, dat / s, dat / s / Ladopted,
             (c >= normLo && c <= normHi) ? "   <- normalisation region" : "");
      if (c >= normLo && c <= normHi)
         flat.push_back(dat / s);
   }
   if (!flat.empty()) {
      std::sort(flat.begin(), flat.end());
      printf("\n  over %.0f-%.0f deg: median %.2f, spread %.2f to %.2f (%+.0f%% / %+.0f%%)\n", normLo, normHi,
             flat[flat.size() / 2], flat.front(), flat.back(), 100 * (flat.front() / flat[flat.size() / 2] - 1),
             100 * (flat.back() / flat[flat.size() / 2] - 1));
   }

   TCanvas *cv = new TCanvas("cL", "L vs angle", 1300, 540);
   cv->Divide(2, 1);
   cv->cd(1);
   gPad->SetLogy();
   gPad->SetGridy();
   auto *f1 = new TH1D("fL1", "elastic: measured vs KD03;#theta_{cm} [deg];d#sigma/d#Omega  [arb  /  mb sr^{-1}]", 1,
                       20, 140);
   f1->SetMinimum(0.3);
   f1->SetMaximum(3e4);
   f1->Draw();
   gD->SetMarkerStyle(20);
   gD->SetLineWidth(2);
   gD->Draw("LP");
   gT->SetLineColor(kRed + 1);
   gT->SetLineWidth(3);
   gT->Draw("L");
   auto *lg = new TLegend(0.50, 0.72, 0.88, 0.88);
   lg->AddEntry(gD, "data (arbitrary units)", "lp");
   lg->AddEntry(gT, "KD03 [mb/sr]", "l");
   lg->Draw();

   cv->cd(2);
   gPad->SetLogy();
   gPad->SetGridy();
   auto *f2 = new TH1D("fL2", "L = data / KD03   (must be CONSTANT);#theta_{cm} [deg];L [counts/mb]", 1, 20, 140);
   f2->SetMinimum(3);
   f2->SetMaximum(400);
   f2->Draw();
   gL->SetMarkerStyle(20);
   gL->SetMarkerSize(1.2);
   gL->SetMarkerColor(kAzure + 2);
   gL->SetLineColor(kAzure + 2);
   gL->SetLineWidth(2);
   gL->Draw("LP");
   auto *ln = new TLine(20, Ladopted, 140, Ladopted);
   ln->SetLineColor(kRed + 1);
   ln->SetLineWidth(3);
   ln->SetLineStyle(2);
   ln->Draw();
   auto *bx = new TBox(normLo, 3, normHi, 400);
   bx->SetFillColorAlpha(kGreen + 1, 0.15);
   bx->Draw();
   gL->Draw("LP");
   auto *lg2 = new TLegend(0.42, 0.74, 0.88, 0.88);
   lg2->AddEntry(gL, "L(#theta)", "lp");
   lg2->AddEntry(ln, TString::Format("adopted L = %.1f", Ladopted), "l");
   lg2->AddEntry(bx, TString::Format("normalisation %.0f-%.0f#circ", normLo, normHi), "f");
   lg2->Draw();

   TString png = here + "/plots/Ldep_C14" + tag + ".png";
   cv->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
