/// @file exz_view_C14.C
/// @brief The Ex versus vertex-z plane, which is what sets the usable target thickness.
///
/// This is the plot the excited-state normalisation rests on, so it is worth looking at directly
/// rather than inferring from summary numbers.
///
/// The elastic fills the whole chamber; the excited states do not. That is not physics -- the beam
/// loses only ~12 MeV crossing the full metre of H2 at 300 torr, which moves E_cm by less than
/// 1 MeV and closes no channel -- it is the trigger, which the simulation does not contain and
/// which therefore cannot be corrected for, only avoided.
///
/// It is avoided by working in the flat part of the distribution, where the trigger efficiency is
/// constant and cancels in an excited-to-elastic ratio. The fourth panel is the one that decides
/// where that is: the excited-state yield DIVIDED BY the elastic yield, per z bin. Wherever that
/// ratio is flat, the efficiency is common to both and drops out. Where it slopes, it does not.
///
///   root -b -q 'exz_view_C14.C()'

void exz_view_C14(TString cache = "plots/proton_kin_300gfx_ex.root", Double_t zLo = 10.0, Double_t zHi = 400.0,
                  TString tag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *f = TFile::Open(here + "/" + cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", cache.Data());
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t)
      return;

   const double ZMIN = -100, ZMAX = 1100;
   const int NZ = 48; // 25 mm bins

   auto proj = [&](const char *nm, const char *cut) {
      auto *h = new TH1D(nm, "", NZ, ZMIN, ZMAX);
      t->Draw(TString::Format("vertexz>>%s", nm), cut, "goff");
      h->SetDirectory(nullptr);
      return h;
   };
   TH1D *hEl = proj("zEl", "ex>-0.4&&ex<0.6");
   TH1D *hMu = proj("zMu", "ex>5.4&&ex<8.0");
   TH1D *h85 = proj("z85", "ex>8.2&&ex<8.9");
   TH1D *h93 = proj("z93", "ex>9.0&&ex<9.8");

   TCanvas *c = new TCanvas("cz", "Ex vs z", 1500, 950);
   c->Divide(2, 2);

   // ---- 1. the plane
   c->cd(1);
   gPad->SetLogz();
   auto *h2 = new TH2D("h2", "E_{x} vs vertex z;vertex z [mm];E_{x} [MeV]", 60, ZMIN, ZMAX, 110, -1.5, 11.0);
   t->Draw("ex:vertexz>>h2", "", "goff");
   h2->SetDirectory(nullptr);
   h2->Draw("colz");
   for (double z : {zLo, zHi}) {
      auto *l = new TLine(z, -1.5, z, 11.0);
      l->SetLineColor(kRed + 1);
      l->SetLineWidth(3);
      l->SetLineStyle(2);
      l->Draw();
   }

   // ---- 2. the z projections, raw
   c->cd(2);
   gPad->SetLogy();
   hEl->SetTitle("vertex z by E_{x} window;vertex z [mm];counts / 25 mm");
   int col[4] = {kBlack, kRed + 1, kAzure + 2, kGreen + 3};
   TH1D *hh[4] = {hEl, hMu, h85, h93};
   const char *nm[4] = {"elastic", "multiplet 5.4-8.0", "8.2-8.9", "9.0-9.8"};
   double mx = 0;
   for (auto *h : hh)
      mx = std::max(mx, h->GetMaximum());
   hEl->SetMaximum(mx * 3);
   hEl->SetMinimum(0.5);
   auto *lg = new TLegend(0.58, 0.68, 0.89, 0.88);
   for (int i = 0; i < 4; ++i) {
      hh[i]->SetLineColor(col[i]);
      hh[i]->SetLineWidth(2);
      hh[i]->Draw(i ? "hist same" : "hist");
      lg->AddEntry(hh[i], nm[i], "l");
   }
   lg->Draw();
   for (double z : {zLo, zHi}) {
      auto *l = new TLine(z, 0.5, z, mx * 3);
      l->SetLineColor(kRed + 1);
      l->SetLineWidth(3);
      l->SetLineStyle(2);
      l->Draw();
   }

   // ---- 3. each shape normalised to unit area: does the SHAPE differ, or only the rate?
   c->cd(3);
   auto *fr = new TH1D("frz", "same, each normalised to unit area;vertex z [mm];fraction / 25 mm", NZ, ZMIN, ZMAX);
   double my = 0;
   TH1D *hn[4];
   for (int i = 0; i < 4; ++i) {
      hn[i] = (TH1D *)hh[i]->Clone(TString::Format("n%d", i));
      hn[i]->SetDirectory(nullptr);
      if (hn[i]->Integral() > 0)
         hn[i]->Scale(1.0 / hn[i]->Integral());
      my = std::max(my, hn[i]->GetMaximum());
   }
   fr->SetMaximum(my * 1.35);
   fr->SetMinimum(0);
   fr->Draw();
   for (int i = 0; i < 4; ++i)
      hn[i]->Draw("hist same");
   lg->Draw();
   for (double z : {zLo, zHi}) {
      auto *l = new TLine(z, 0, z, my * 1.35);
      l->SetLineColor(kRed + 1);
      l->SetLineWidth(3);
      l->SetLineStyle(2);
      l->Draw();
   }

   // ---- 4. the decisive one: excited / elastic per z bin
   c->cd(4);
   auto *rMu = (TH1D *)hMu->Clone("rMu");
   rMu->SetDirectory(nullptr);
   rMu->Divide(hMu, hEl, 1, 1, "B");
   auto *r85 = (TH1D *)h85->Clone("r85");
   r85->SetDirectory(nullptr);
   r85->Divide(h85, hEl, 1, 1, "B");
   rMu->SetTitle("excited / elastic per z bin -- flat means the trigger cancels;"
                 "vertex z [mm];ratio");
   rMu->GetXaxis()->SetRangeUser(-50, 800);
   rMu->SetMaximum(0.6);
   rMu->SetMinimum(0);
   rMu->SetLineColor(kRed + 1);
   rMu->SetMarkerColor(kRed + 1);
   rMu->SetMarkerStyle(20);
   rMu->Draw("E1");
   r85->SetLineColor(kAzure + 2);
   r85->SetMarkerColor(kAzure + 2);
   r85->SetMarkerStyle(21);
   r85->Draw("E1 same");
   for (double z : {zLo, zHi}) {
      auto *l = new TLine(z, 0, z, 0.6);
      l->SetLineColor(kRed + 1);
      l->SetLineWidth(3);
      l->SetLineStyle(2);
      l->Draw();
   }
   auto *lg2 = new TLegend(0.55, 0.72, 0.89, 0.88);
   lg2->AddEntry(rMu, "multiplet / elastic", "lp");
   lg2->AddEntry(r85, "8.2-8.9 / elastic", "lp");
   lg2->Draw();

   // ---- the numbers behind the picture
   printf("\n  z [mm]     | elastic | multiplet |  8.2-8.9 | mult/elastic\n");
   for (int b = 1; b <= NZ; ++b) {
      double lo = hEl->GetBinLowEdge(b), e = hEl->GetBinContent(b), m = hMu->GetBinContent(b);
      if (lo < -50 || lo > 800)
         continue;
      printf("  %5.0f-%5.0f | %7.0f | %9.0f | %8.0f | %s%.3f%s\n", lo, lo + hEl->GetBinWidth(b), e, m,
             h85->GetBinContent(b), (lo >= zLo && lo < zHi) ? "\033[1m" : "", e > 0 ? m / e : 0,
             (lo >= zLo && lo < zHi) ? "\033[0m" : "");
   }
   auto frac = [&](TH1D *h) {
      double in = 0, all = h->Integral();
      for (int b = 1; b <= NZ; ++b) {
         double c = h->GetBinCenter(b);
         if (c >= zLo && c <= zHi)
            in += h->GetBinContent(b);
      }
      return all > 0 ? in / all : 0;
   };
   printf("\n  fraction inside the %.0f-%.0f mm window:  elastic %.3f   multiplet %.3f   8.2-8.9 %.3f\n", zLo, zHi,
          frac(hEl), frac(hMu), frac(h85));

   TString png = here + "/plots/exz_view_C14" + tag + ".png";
   c->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
