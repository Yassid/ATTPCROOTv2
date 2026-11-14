// cmp_prev_curr_dp.C -- PREVIOUS (Spyral-1.0 InterpSolver, ~/C17_dp_fits + Movilla notebook) vs
// CURRENT (genfit find805 + CATIMA) 16C(d,p)17C analysis, side by side.
//   (a) inclusive Ex, whole chamber, theta_cm 10-40, each analysis with ITS OWN selection + Ex recipe
//   (b) the same, vz 200-400 only
//   (c) bound-state yields per 5 deg bin (previous = fitted g.s.+0.217+0.335 sum; current = window - sideband)
//   (d) bound-state dsigma/dOmega, every version + Pereira-Lopez PLB 811 (17.2 MeV/u, other energy)
//   root -l -b -q cmp_prev_curr_dp.C
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TSystem.h>
#include <fstream>
#include <sstream>
#include <map>

namespace prev { // copied from ~/C17_dp_fits/Codes/C16_pd_AngBins.C (nominal cut set)
constexpr double U = 931.49401, Eb = 188.33;
constexpr double Mp = 1.007825 * U, Md = 2.0135532 * U, M16 = 16.0147 * U, M17 = 17.0226 * U;
constexpr double kP0 = 0.618451, kP1 = 0.0080509149, kOffset = 1.022, kShift = 0.208;
double om(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z); }
bool kine(double thl, double Ke, double &ex, double &tcm)
{
   double m1 = M16, m2 = Md, m3 = Mp, m4 = M17, Et1 = Eb + m1, Et3 = Ke + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1, u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * om(s, m1 * m1, m2 * m2) * om(u, m2 * m2, m3 * m3) - (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) + s + u - m2 * m2;
   if (a < 0) return false;
   double m4x = std::sqrt(a); ex = m4x - m4;
   double t = m2 * m2 + m4x * m4x - 2 * m2 * Et4;
   double c = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) + (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
              (om(s, m1 * m1, m2 * m2) * om(s, m3 * m3, m4x * m4x));
   c = std::max(-1.0, std::min(1.0, c));
   tcm = 180.0 - std::acos(c) * TMath::RadToDeg();
   return true;
}
} // namespace prev

using namespace exsel;
static double EbVz(double vz) { return 184.25 - 0.011030 * (vz - 500); }

struct Row { std::map<std::string, std::string> f; };
static std::vector<Row> readCsv(const char *fn)
{
   std::vector<Row> out; std::ifstream in(fn); std::string line; std::vector<std::string> hdr;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> tok; std::stringstream ss(line); std::string x;
      while (std::getline(ss, x, ',')) tok.push_back(x);
      if (hdr.empty()) { hdr = tok; continue; }
      Row r; for (size_t i = 0; i < hdr.size() && i < tok.size(); ++i) r.f[hdr[i]] = tok[i];
      out.push_back(r);
   }
   return out;
}
static double D(const Row &r, const char *k) { return std::atof(r.f.at(k).c_str()); }

void cmp_prev_curr_dp()
{
   gStyle->SetOptStat(0);
   const double exLo = -2, exHi = 8; const int nb = 100;
   TH1D *hP = new TH1D("hP", "", nb, exLo, exHi), *hC = new TH1D("hC", "", nb, exLo, exHi);
   TH1D *hPw = new TH1D("hPw", "", nb, exLo, exHi), *hCw = new TH1D("hCw", "", nb, exLo, exHi);
   long nPrev = 0, nCurr = 0;
   // ---- previous: Spyral InterpSolver, polar>=90, Ebeam 188.33, linear z correction
   {
      const char *dir = "/home/yassid/C16_dp/C16_dp/InterpSolver_root";
      void *d = gSystem->OpenDirectory(dir); const char *fn;
      while ((fn = gSystem->GetDirEntry(d))) {
         TString s(fn); if (!s.EndsWith("_1H.root")) continue;
         TFile f(TString(dir) + "/" + s); TTree *t = (TTree *)f.Get("parquettree");
         double pol, br, z;
         t->SetBranchAddress("polar", &pol); t->SetBranchAddress("brho", &br); t->SetBranchAddress("vertex_z", &z);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (pol * TMath::RadToDeg() < 90.0) continue;
            double p = br * 2.99792458e2, ke = std::sqrt(p * p + prev::Mp * prev::Mp) - prev::Mp, ex, tc;
            if (!prev::kine(pol, ke, ex, tc)) continue;
            double zcm = z * 100.0, exc = ex - prev::kP1 * zcm - prev::kP0 + prev::kOffset + prev::kShift;
            if (tc < 10 || tc >= 40) continue;
            hP->Fill(exc); ++nPrev;
            if (z * 1000 >= 200 && z * 1000 < 400) hPw->Fill(exc);
         }
      }
      gSystem->FreeDirectory(d);
   }
   // ---- current: genfit find805, adopted bound-AD selection (Ebeam(vz) pivot, keOff 0, kc off)
   {
      TFile f("/mnt/f/a1975/caches/.explorer/exp_dp_find805.root"); TTree *t = (TTree *)f.Get("pk");
      float ke, th, vz, c2;
      t->SetBranchAddress("ke", &ke); t->SetBranchAddress("theta", &th); t->SetBranchAddress("vertexz", &vz); t->SetBranchAddress("chi2ndf", &c2);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!(c2 <= 13.5) || !(th >= 89 && th <= 161) || !(ke >= 1.2) || !(vz >= -100 && vz <= 940)) continue;
         double ex, tc; if (!kine2b(EbVz(vz), th * TMath::DegToRad(), ke, ex, tc)) continue;
         if (tc < 10 || tc >= 40) continue;
         hC->Fill(ex); ++nCurr;
         if (vz >= 200 && vz < 400) hCw->Fill(ex);
      }
   }
   auto inwin = [](TH1D *h, double a, double b) { return h->Integral(h->FindBin(a + 1e-6), h->FindBin(b - 1e-6)); };
   printf("\nINCLUSIVE theta_cm 10-40   entries: previous %ld  current %ld\n", nPrev, nCurr);
   printf("  Ex -1.5..0.73  (bound)   prev %6.0f  curr %6.0f   | vz 200-400: prev %5.0f  curr %5.0f\n",
          inwin(hP, -1.5, 0.73), inwin(hC, -1.5, 0.73), inwin(hPw, -1.5, 0.73), inwin(hCw, -1.5, 0.73));
   printf("  Ex  0.73..7    (unbound) prev %6.0f  curr %6.0f   | vz 200-400: prev %5.0f  curr %5.0f\n",
          inwin(hP, 0.73, 7), inwin(hC, 0.73, 7), inwin(hPw, 0.73, 7), inwin(hCw, 0.73, 7));
   for (auto *h : {hP, hC}) {
      int ib = h->GetMaximumBin(); printf("  %s peak bin at Ex %+.2f MeV (max in -1..1: ", h == hP ? "prev" : "curr", h->GetBinCenter(ib));
      h->GetXaxis()->SetRangeUser(-1, 1); printf("%+.2f)\n", h->GetBinCenter(h->GetMaximumBin())); h->GetXaxis()->SetRange(0, 0);
   }

   // ---- bound angular distributions
   // previous repo: fine (2.5 deg) fitted, sum of the three bound states; merge pairs -> 5 deg
   auto fine = readCsv("/home/yassid/C17_dp_fits/results/fine/dsdo.csv");
   auto fy = readCsv("/home/yassid/C17_dp_fits/results/fine/yields.csv");
   auto isB = [](const Row &r) { auto s = r.f.at("state"); return s == "gs_3half+" || s == "ex_1half+_217" || s == "ex_5half+_335"; };
   double pX[12], pY[12], pE[12], y5[6] = {0}, y5e[6] = {0};
   for (int i = 0; i < 12; ++i) { pX[i] = 11.25 + 2.5 * i; pY[i] = 0; pE[i] = 0; }
   for (auto &r : fine) if (isB(r)) { int i = (int)std::lround((D(r, "bin_center") - 11.25) / 2.5); pY[i] += D(r, "dsdo_mbsr"); pE[i] += std::pow(D(r, "dsdo_err_mbsr"), 2); }
   for (int i = 0; i < 12; ++i) pE[i] = std::sqrt(pE[i]);
   for (auto &r : fy) if (isB(r)) { int i = (int)((D(r, "bin_lo") - 10) / 5); y5[i] += D(r, "counts"); y5e[i] += std::pow(D(r, "counts_err"), 2); }
   // Movilla notebook (far-end 0.9-1.0 m slice x10, Ex-only eps), its own printed numbers
   double mX[6] = {12.5, 17.5, 22.5, 27.5, 32.5, 37.5}, mY[6] = {23.1, 16.8, 3.67, 3.66, 4.16, 3.40}, mC[6] = {72, 70, 21, 25, 32, 27}, mE[6], z6[12] = {0};
   for (int i = 0; i < 6; ++i) mE[i] = mY[i] / std::sqrt(mC[i]);
   // current
   auto cur = readCsv("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/dp_bound_angdist_ADOPTED.csv");
   auto cfr = readCsv("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/dp_bound_angdist_vz0-915.csv");
   TGraphErrors *gA = new TGraphErrors, *gAx = new TGraphErrors, *gF = new TGraphErrors, *gYc = new TGraphErrors, *gYw = new TGraphErrors;
   for (auto &r : cur) {
      auto *g = D(r, "use") > 0 ? gA : gAx; int n = g->GetN();
      g->SetPoint(n, D(r, "thcm_mid") + 0.4, D(r, "dsdo_mb_sr")); g->SetPointError(n, 0, D(r, "err_mb_sr"));
      if (D(r, "thcm_mid") < 40) { int k = gYw->GetN(); gYw->SetPoint(k, D(r, "thcm_mid") + 0.6, D(r, "Y") - D(r, "bg")); gYw->SetPointError(k, 0, std::sqrt(D(r, "Y") + D(r, "bg"))); }
   }
   for (auto &r : cfr) if (D(r, "thcm_mid") < 40) {
      int n = gF->GetN(); gF->SetPoint(n, D(r, "thcm_mid") - 0.4, D(r, "dsdo_mb_sr")); gF->SetPointError(n, 0, D(r, "err_mb_sr"));
      int k = gYc->GetN(); gYc->SetPoint(k, D(r, "thcm_mid") + 0.2, D(r, "Y") - D(r, "bg")); gYc->SetPointError(k, 0, std::sqrt(D(r, "Y") + D(r, "bg")));
   }
   TGraphErrors *gY5 = new TGraphErrors;
   for (int i = 0; i < 6; ++i) { gY5->SetPoint(i, 12.5 + 5 * i - 0.6, y5[i]); gY5->SetPointError(i, 0, std::sqrt(y5e[i])); }
   // PLB 811 bound sum (17.2 MeV/u)
   TGraphErrors *gPLB = new TGraphErrors;
   { std::ifstream in("/home/yassid/fair_install/16C_dp/RCS/dp_AngularDistro_boundstates_Eex_no56_7May2018_CM.dat"); double a, b, c; while (in >> a >> b >> c) { int n = gPLB->GetN(); gPLB->SetPoint(n, a, b); gPLB->SetPointError(n, 0, c); } }

   printf("\nBOUND dsigma/dOmega (mb/sr), 5 deg bins 10-40\n  thcm   prev-repo(fit,5deg avg)  notebook(far-end x10)  CURRENT adopted(vz200-400)  current full(vz0-915)\n");
   for (int i = 0; i < 6; ++i) {
      double a = 0.5 * (pY[2 * i] + pY[2 * i + 1]), ae = 0.5 * std::hypot(pE[2 * i], pE[2 * i + 1]);
      printf("  %4.1f   %5.2f(%4.2f)               %5.2f(%4.2f)          %5.2f(%4.2f)%s           %5.2f(%4.2f)\n", mX[i], a, ae, mY[i], mE[i],
             D(cur[i], "dsdo_mb_sr"), D(cur[i], "err_mb_sr"), D(cur[i], "use") > 0 ? " " : "*", D(cfr[i], "dsdo_mb_sr"), D(cfr[i], "err_mb_sr"));
   }
   printf("  (* = excluded from the adopted result)\n\nBOUND YIELDS (counts)\n  thcm   prev-repo fit (full chamber)   curr window-bg vz0-915   curr vz200-400\n");
   for (int i = 0; i < 6; ++i)
      printf("  %4.1f   %6.1f(%4.1f)                  %6.1f                   %6.1f\n", mX[i], y5[i], std::sqrt(y5e[i]),
             D(cfr[i], "Y") - D(cfr[i], "bg"), D(cur[i], "Y") - D(cur[i], "bg"));

   // ---- drawing
   TCanvas *c = new TCanvas("c", "", 1500, 1100); c->Divide(2, 2);
   const double res[] = {2.763, 2.980, 3.661, 4.231, 4.841, 5.91, 6.30};
   auto drawEx = [&](int pad, TH1D *a, TH1D *b, const char *tt) {
      c->cd(pad); gPad->SetGrid();
      a->SetLineColor(kRed + 1); a->SetLineWidth(2); b->SetLineColor(kBlue + 1); b->SetLineWidth(2);
      a->SetTitle(TString::Format("%s;E_{x} (MeV);counts / 100 keV", tt));
      a->SetMaximum(1.15 * std::max(a->GetMaximum(), b->GetMaximum())); a->Draw("hist"); b->Draw("hist same");
      for (double e : res) { TLine *l = new TLine(e, 0, e, 0.12 * a->GetMaximum()); l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw(); }
      TLine *sn = new TLine(0.729, 0, 0.729, a->GetMaximum()); sn->SetLineColor(kGreen + 2); sn->SetLineStyle(7); sn->Draw();
      TLegend *L = new TLegend(0.50, 0.68, 0.89, 0.88);
      L->AddEntry(a, TString::Format("previous: Spyral, E_{b} 188.33 + z corr (%.0f)", a->Integral()), "l");
      L->AddEntry(b, TString::Format("current: genfit, E_{b}(vz) (%.0f)", b->Integral()), "l");
      L->AddEntry(sn, "S_{n} 0.729", "l"); L->Draw();
   };
   drawEx(1, hP, hC, "Inclusive, #theta_{cm} 10-40, whole chamber (prev resonances dashed)");
   drawEx(2, hPw, hCw, "Inclusive, #theta_{cm} 10-40, vz 200-400 mm");

   c->cd(3); gPad->SetGrid();
   TH1F *fr = gPad->DrawFrame(8, 0, 42, 1.15 * std::max({gYc->GetHistogram()->GetMaximum(), gY5->GetHistogram()->GetMaximum()}));
   fr->SetTitle("Bound-state yields (g.s.+0.217+0.335);#theta_{cm} (deg);counts per 5 deg");
   gY5->SetMarkerStyle(25); gY5->SetMarkerColor(kRed + 1); gY5->SetLineColor(kRed + 1); gY5->Draw("P");
   gYc->SetMarkerStyle(20); gYc->SetMarkerColor(kBlue + 1); gYc->SetLineColor(kBlue + 1); gYc->Draw("P");
   gYw->SetMarkerStyle(24); gYw->SetMarkerColor(kAzure + 7); gYw->SetLineColor(kAzure + 7); gYw->Draw("P");
   TLegend *L3 = new TLegend(0.12, 0.70, 0.62, 0.88);
   L3->AddEntry(gY5, "previous repo: fitted triplet, #sigma fixed 0.12", "p");
   L3->AddEntry(gYc, "current: window - sideband, vz 0-915", "p");
   L3->AddEntry(gYw, "current: window - sideband, vz 200-400", "p"); L3->Draw();

   c->cd(4); gPad->SetGrid(); gPad->SetLogy();
   TH1F *f4 = gPad->DrawFrame(5, 0.08, 45, 80);
   f4->SetTitle("Bound-state d#sigma/d#Omega;#theta_{cm} (deg);d#sigma/d#Omega (mb/sr)");
   TGraphErrors *gP = new TGraphErrors(12, pX, pY, z6, pE);
   for (int i = 6; i < 12; ++i) gP->SetPointError(i, 0, pE[i]);
   gPLB->SetMarkerStyle(21); gPLB->SetMarkerColor(kGray + 1); gPLB->SetLineColor(kGray + 1); gPLB->Draw("P");
   gP->SetMarkerStyle(25); gP->SetMarkerColor(kRed + 1); gP->SetLineColor(kRed + 1); gP->Draw("P");
   TGraphErrors *gM = new TGraphErrors(6, mX, mY, z6, mE);
   gM->SetMarkerStyle(26); gM->SetMarkerColor(kMagenta + 2); gM->SetLineColor(kMagenta + 2); gM->SetMarkerSize(1.3); gM->Draw("P");
   gF->SetMarkerStyle(24); gF->SetMarkerColor(kAzure + 7); gF->SetLineColor(kAzure + 7); gF->Draw("P");
   gAx->SetMarkerStyle(24); gAx->SetMarkerColor(kBlue + 1); gAx->SetLineColor(kBlue + 1); gAx->Draw("P");
   gA->SetMarkerStyle(20); gA->SetMarkerColor(kBlue + 1); gA->SetLineColor(kBlue + 1); gA->SetMarkerSize(1.3); gA->Draw("P");
   TLegend *L4 = new TLegend(0.40, 0.66, 0.89, 0.89);
   L4->AddEntry(gPLB, "Pereira-Lopez PLB 811 (17.2 MeV/u)", "p");
   L4->AddEntry(gM, "previous: notebook, far-end slice #times10", "p");
   L4->AddEntry(gP, "previous: C17_dp_fits repo, fitted triplet", "p");
   L4->AddEntry(gF, "current: vz 0-915, sim #varepsilon", "p");
   L4->AddEntry(gA, "current ADOPTED: vz 200-400 (open = excluded)", "p"); L4->Draw();
   c->SaveAs("cmp_prev_curr_dp.png");
}
