/// @file locus_fine_C15d.C
/// @brief The gate drawer's [Locus check], finely binned and written to a PNG.
///
///   root -b -q 'pid/locus_fine_C15d.C()'                          // proton gate, 0.5 deg x 0.005 T.m
///   root -b -q 'pid/locus_fine_C15d.C("pid/deuteron_C15d.json",1,2)'
///
/// The GUI panel is limited by the canvas it draws into; this writes the same comparison at
/// whatever binning is asked for, so a band can be judged against the curve rather than eyeballed.
///
/// ★ SAME PLANE AS PRODUCTION. IC [1045, 1225] with single pulse, which is what
/// gate_events_C15d.C applies before fitting. A locus check on any other plane is checking a gate
/// that will not be the one used.
///
/// ★ BOTH POLAR CONVENTIONS ARE DRAWN. The stored polar is either theta_lab or 180 - theta_lab
/// and the two are not distinguishable by eye on a band that lives at one end; drawing both means
/// the convention cannot be wrong silently. The one the band follows is the right one.

#include <fstream>

static double dpBrho(double thlab, double Eb, double exStar, double m3Amu, double m4Amu)
{
   const double u = 931.49401;
   const double m1 = 15.0105993 * u, m2 = 2.014102 * u;
   const double m3 = m3Amu * u, m4 = m4Amu * u + exStar;
   double E1 = Eb + m1, pb = std::sqrt(E1 * E1 - m1 * m1), thr = thlab * TMath::DegToRad(), best = -1;
   for (double ke = 0.2; ke < 170; ke += 0.02) {
      double E3 = ke + m3, p3 = std::sqrt(E3 * E3 - m3 * m3), E4 = E1 + m2 - E3;
      double px = p3 * std::sin(thr), pz = p3 * std::cos(thr);
      double m4x2 = E4 * E4 - (px * px + (pb - pz) * (pb - pz));
      if (m4x2 < 0) continue;
      if (std::fabs(std::sqrt(m4x2) - m4) < 0.06) best = ke;
   }
   if (best < 0) return -1;
   return std::sqrt(2 * m3 * best + best * best) / 299.792458;
}

static TGraph *locus(bool flip, double Eb, double ex, int col, double m3, double m4, int style = 1)
{
   auto *g = new TGraph();
   for (double th = 1; th <= 179; th += 0.5) {
      double b = dpBrho(th, Eb, ex, m3, m4);
      if (b > 0) g->SetPoint(g->GetN(), flip ? 180.0 - th : th, b);
   }
   g->SetLineColor(col);
   g->SetLineWidth(3);
   g->SetLineStyle(style);
   return g;
}

void locus_fine_C15d(TString gateFile = "pid/proton_C15d.json", int Z = 1, int A = 1,
                     double Ebeam = 207.0, double exGuide = 1.766,
                     // The stored polar is 180 - theta_lab on this dataset -- measured, see
                     // pid/open_gate_draw.sh. Solid = the convention in use, grey = the other one.
                     bool flipPolar = true,
                     // fine by default: 0.5 deg x 0.005 T.m
                     int nbx = 360, int nby = 400, double yhi = 2.0,
                     double icLo = 1045, double icHi = 1225,
                     TString pointsFile = "pid/points_C15d.root", TString outPng = "pid/plots/locus_fine_C15d.png")
{
   // ejectile -> residual: the compound is 15C + d, so fixing the ejectile fixes the residual
   double m3 = 1.00782503, m4 = 16.0147413;
   if (Z == 1 && A == 2) { m3 = 2.01410178; m4 = 15.0105993; }
   if (Z == 1 && A == 3) { m3 = 3.01604928; m4 = 14.0032420; }

   // --- the gate polygon -----------------------------------------------------------------------
   std::ifstream in(gateFile.Data());
   if (!in) { printf("\033[1;31mERROR: %s not found\033[0m\n", gateFile.Data()); return; }
   std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto *cut = new TCutG("gate", 0);
   size_t p = all.find("\"vertices\"");
   if (p == std::string::npos) p = all.find("\"points\"");
   int nv = 0;
   if (p != std::string::npos) {
      // Walk EVERY "[x , y]" pair from the key to the closing "]]" of the list. The earlier
      // bracket-hopping version dropped the last vertex, which silently changes the polygon.
      size_t end = all.find("]]", p);
      size_t q = all.find('[', p + 8);
      while (q != std::string::npos && (end == std::string::npos || q < end)) {
         double x, y;
         if (sscanf(all.c_str() + q + 1, " %lf , %lf", &x, &y) == 2) cut->SetPoint(nv++, x, y);
         q = all.find('[', q + 1);
      }
   }
   cut->Set(nv);
   printf("gate %s: %d vertices\n", gateFile.Data(), nv);
   if (nv < 3) { printf("\033[1;31mERROR: could not parse the polygon\033[0m\n"); return; }

   TFile *fp = TFile::Open(pointsFile);
   TTree *t = fp ? (TTree *)fp->Get("pts") : nullptr;
   if (!t) { printf("\033[1;31mERROR: no 'pts' in %s\033[0m\n", pointsFile.Data()); return; }

   Float_t sq, br, pol, ic;
   Int_t np;
   t->SetBranchAddress("sqrtdedx", &sq);
   t->SetBranchAddress("brho", &br);
   t->SetBranchAddress("polar", &pol);
   t->SetBranchAddress("ic", &ic);
   t->SetBranchAddress("npulse", &np);

   auto *h = new TH2D("hloc", Form("%s inside the gate, IC [%.0f,%.0f] single pulse"
                                   ";polar [deg];B#rho [T#upointm]",
                                   gateFile.Contains("proton") ? "15C(d,p) protons" : "gated tracks", icLo, icHi),
                      nbx, 0, 180, nby, 0, yhi);
   // The input TFile is the current directory here, so the histogram belongs to it and Close()
   // below would delete it -- the draw then segfaults on freed memory.
   h->SetDirectory(nullptr);
   Long64_t nAll = 0, nIn = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(ic >= icLo && ic <= icHi && np == 1)) continue;
      ++nAll;
      if (!cut->IsInside(sq, br)) continue;
      ++nIn;
      h->Fill(pol, br);
   }
   fp->Close();
   printf("IC-gated single pulse: %lld ; inside the gate: %lld (%.2f %%)\n", nAll, nIn, 100.0 * nIn / nAll);
   printf("binning: %d x %d  (%.2f deg x %.4f T.m)\n", nbx, nby, 180.0 / nbx, yhi / nby);

   gStyle->SetOptStat(0);
   auto *c1 = new TCanvas("c1", "locus", 1400, 950);
   c1->SetRightMargin(0.13);
   c1->SetLogz();
   h->Draw("colz");

   auto *g0 = locus(flipPolar, Ebeam, 0.0, kRed + 1, m3, m4);
   auto *gA = locus(!flipPolar, Ebeam, 0.0, kGray + 2, m3, m4);
   auto *gE = (exGuide > 0) ? locus(flipPolar, Ebeam, exGuide, kMagenta + 1, m3, m4, 7) : nullptr;
   gA->Draw("L same");
   g0->Draw("L same");
   if (gE) gE->Draw("L same");

   auto *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(g0, Form("g.s. locus, E_{b} = %.1f MeV", Ebeam), "l");
   if (gE) leg->AddEntry(gE, Form("E_{x} = %.3f MeV", exGuide), "l");
   leg->AddEntry(gA, "other polar convention", "l");
   leg->Draw();

   gSystem->mkdir(gSystem->DirName(outPng), kTRUE);
   c1->SaveAs(outPng);
   printf("\033[1;32mwrote\033[0m %s\n", outPng.Data());
}
