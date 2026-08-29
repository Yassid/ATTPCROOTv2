/// @file channel_leverage_C14.C
/// @brief Which channel would actually profit from finer pads, a higher field, or the vertex beam
/// energy -- computed from kinematics alone, before spending a single CPU hour on a simulation.
///
///   root -b -q 'channel_leverage_C14.C()'
///
/// THE QUESTION. On 14C(p,p') the field x pitch matrix bought almost nothing, because the recoil
/// protons the experiment accepts come out near theta_lab 77 deg where the excitation energy is
/// nearly insensitive to everything. That is a property of THAT channel's kinematics, not of the
/// detector. A transfer reaction whose ejectile goes BACKWARD in the laboratory samples a
/// completely different part of the same curves.
///
/// WHAT IS COMPUTED. For each channel, along the ground-state locus:
///     dEx/dE_beam   the leverage of the beam energy at the vertex
///     dEx/dKE       the leverage of the ejectile energy, i.e. of tracking
///     dEx/dtheta    the leverage of the ejectile angle
/// and then the excitation-energy resolution those imply when folded with the resolutions the
/// campaign MEASURED, for the worst and best cells of the matrix.
///
/// CAVEAT, and it is not small: sigma(KE) and sigma(theta) are carried over from the 14C(p,p')
/// simulation, where they were measured on protons of 0-35 MeV over theta_lab 20-90 deg. A
/// backward-going transfer proton of 2-8 MeV is a different track -- shorter, more curled, and at
/// 7 T possibly not reconstructed at all (the campaign's acceptance hole sits exactly there). So
/// what follows identifies WHERE the leverage is; it does not replace the simulation that would
/// confirm the detector can use it.

#include <vector>

static double le_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
/// Ex from the ejectile (KE, theta_lab); m4 is the GROUND-STATE residual mass
static double le_ex(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(th) * le_om2(s, m1 * m1, m2 * m2) * le_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                 (2 * m2 * m2) +
              s + uu - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

/// Forward kinematics: the ejectile (KE, theta_lab) at a given centre-of-mass angle.
/// theta_cm is the angle of the EJECTILE in the centre of mass.
static bool le_forward(double m1, double m2, double m3, double m4, double K, double thcm_deg, double &ke, double &thlab)
{
   double E1 = K + m1;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1;
   double rs = std::sqrt(s);
   if (rs < m3 + m4)
      return false;
   // ejectile momentum and energy in the CM
   double pcm = le_om2(s, m3 * m3, m4 * m4) / (2 * rs);
   double E3cm = std::sqrt(pcm * pcm + m3 * m3);
   // boost of the CM in the lab
   double plab1 = std::sqrt(E1 * E1 - m1 * m1);
   double beta = plab1 / (E1 + m2), gamma = 1.0 / std::sqrt(1 - beta * beta);
   double th = thcm_deg * TMath::DegToRad();
   double pz = gamma * (pcm * std::cos(th) + beta * E3cm);
   double pt = pcm * std::sin(th);
   double E3 = gamma * (E3cm + beta * pcm * std::cos(th));
   ke = E3 - m3;
   thlab = std::atan2(pt, pz) * TMath::RadToDeg();
   return ke > 0;
}

/// sigma of (KE_reco - KE_true) and of (theta_reco - theta_true) as a function of the TRUE proton
/// energy, measured on the campaign's own samples. Energy is the right matching variable: two
/// protons of the same energy make the same track whatever reaction produced them, whereas two
/// protons at the same lab angle from different channels do not.
struct KeRes {
   std::vector<double> keLo, keHi, sKE, sTh;
   double keMax{0};
   bool ok{false};
};
static double lv_q(std::vector<double> v, double p)
{
   if (v.size() < 25) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static KeRes lv_resolution(const TString &root, const char *cfg)
{
   KeRes r;
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/exres_gs_s*_%s.root 2>/dev/null | head -1", root.Data(), cfg, cfg));
   f = f.Strip(TString::kBoth);
   if (f.IsNull()) return r;
   TFile *fr = TFile::Open(f);
   if (!fr || fr->IsZombie()) return r;
   TTree *t = (TTree *)fr->Get("res");
   if (!t) { fr->Close(); return r; }
   double keTrue, keReco, thTrue, thReco;
   t->SetBranchAddress("keTrue", &keTrue);
   t->SetBranchAddress("keReco", &keReco);
   t->SetBranchAddress("thTrue", &thTrue);
   t->SetBranchAddress("thReco", &thReco);
   const int NB = 8;
   const double edge[NB + 1] = {0, 2, 4, 6, 9, 13, 18, 25, 40};
   std::vector<double> dk[NB], dt[NB];
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      for (int b = 0; b < NB; ++b)
         if (keTrue >= edge[b] && keTrue < edge[b + 1]) {
            dk[b].push_back(keReco - keTrue);
            dt[b].push_back(thReco - thTrue);
         }
   }
   for (int b = 0; b < NB; ++b) {
      double s1 = (lv_q(dk[b], .75) - lv_q(dk[b], .25)) / 1.349;
      double s2 = (lv_q(dt[b], .75) - lv_q(dt[b], .25)) / 1.349;
      if (std::isnan(s1)) continue;
      r.keLo.push_back(edge[b]);
      r.keHi.push_back(edge[b + 1]);
      r.sKE.push_back(s1);
      r.sTh.push_back(std::isnan(s2) ? 0.1 : s2);
      r.keMax = edge[b + 1];
   }
   r.ok = r.sKE.size() >= 3;
   fr->Close();
   return r;
}
/// look up the measured resolution at a proton energy; clamps at the ends and says so
static void lv_lookup(const KeRes &r, double ke, double &sKE, double &sTh, bool &extrap)
{
   extrap = false;
   if (!r.ok) { sKE = 0.1; sTh = 0.1; extrap = true; return; }
   for (size_t i = 0; i < r.sKE.size(); ++i)
      if (ke >= r.keLo[i] && ke < r.keHi[i]) { sKE = r.sKE[i]; sTh = r.sTh[i]; return; }
   extrap = true;
   size_t i = (ke < r.keLo.front()) ? 0 : r.sKE.size() - 1;
   sKE = r.sKE[i];
   sTh = r.sTh[i];
}

void channel_leverage_C14(Double_t EbeamPerU = 11.5, TString root = "/mnt/f/a1954_C14_hf", TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);
   const double u = 931.49401;
   const double m_p = 1.007825 * u, m_d = 2.014102 * u, m_t = 3.016049 * u;
   const double m_C14 = 14.003242 * u, m_C15 = 15.010599 * u, m_C13 = 13.003355 * u;

   struct Chan {
      const char *name;
      double m1, m2, m3, m4;
      int col, sty;
   };
   const double Ebeam = EbeamPerU * 14.0; // 14C beam in every case
   std::vector<Chan> ch = {
      {"14C(p,p')14C", m_C14, m_p, m_p, m_C14, kAzure + 2, 1},
      {"14C(d,p)15C", m_C14, m_d, m_p, m_C15, kRed + 1, 1},
      {"14C(d,t)13C", m_C14, m_d, m_t, m_C13, kOrange + 8, 1},
   };
   printf("\n  14C beam at %.1f MeV/u = %.1f MeV\n", EbeamPerU, Ebeam);
   for (auto &c : ch) {
      double Q = (c.m1 + c.m2 - c.m3 - c.m4);
      printf("  %-14s Q = %+7.3f MeV\n", c.name, Q);
   }

   // The resolutions come from the campaign's own samples, per proton ENERGY rather than as one
   // number: sigma(KE) runs from ~0.02 MeV for a 2 MeV proton to several hundred keV at 30 MeV,
   // and a single value would decide the answer by itself.
   KeRes resW = lv_resolution(root, "b285_attpc"); // today
   KeRes resB = lv_resolution(root, "b400_2mm");   // best cell of the matrix
   printf("\n  measured resolution vs proton energy (IQR/1.349)\n");
   printf("  %-12s", "KE [MeV]");
   for (size_t i = 0; i < resW.sKE.size(); ++i) printf(" %5.0f-%-4.0f", resW.keLo[i], resW.keHi[i]);
   printf("\n  %-12s", "2.85T AT-TPC");
   for (size_t i = 0; i < resW.sKE.size(); ++i) printf(" %10.3f", resW.sKE[i]);
   printf("\n  %-12s", "4T 2mm");
   for (size_t i = 0; i < resB.sKE.size(); ++i) printf(" %10.3f", resB.sKE[i]);
   printf("\n");
   const double sEbUncorr = 2.5; // MeV rms, from a vertex uniform over the drift
   const double sEbCorr = 0.32;  // MeV, beam-energy straggling about the profile

   std::vector<TGraph *> gLoc, gEb, gKe, gTh, gResW, gResB;
   for (size_t i = 0; i < ch.size(); ++i) {
      gLoc.push_back(new TGraph());
      gEb.push_back(new TGraph());
      gKe.push_back(new TGraph());
      gTh.push_back(new TGraph());
      gResW.push_back(new TGraph());
      gResB.push_back(new TGraph());
   }

   for (size_t i = 0; i < ch.size(); ++i) {
      auto &c = ch[i];
      printf("\n  ===== %s =====\n", c.name);
      printf("  %8s %9s %9s | %11s %10s %11s | %12s %12s\n", "theta_cm", "theta_lab", "KE [MeV]", "dEx/dEbeam",
             "dEx/dKE", "dEx/dth[/deg]", "2.85T AT-TPC", "4T 2mm+corr");
      // THE ANALYSIS CONVENTION. le_forward takes the EJECTILE angle in the centre of mass, which
      // is 180 deg minus the theta_cm every a1954/a1975 macro quotes (there, small theta_cm is the
      // forward-going HEAVY product, i.e. a backward-going light ejectile). Everything below is
      // reported in the analysis convention so it can be laid against a measured distribution.
      for (double tcmA = 4; tcmA <= 176; tcmA += 2) {
         double tcm = 180.0 - tcmA;
         double ke, thl;
         if (!le_forward(c.m1, c.m2, c.m3, c.m4, Ebeam, tcm, ke, thl))
            continue;
         if (ke < 0.3 || thl < 3 || thl > 177)
            continue;
         double th = thl * TMath::DegToRad();
         double dEb = le_ex(c.m1, c.m2, c.m3, c.m4, Ebeam + 0.5, th, ke) -
                      le_ex(c.m1, c.m2, c.m3, c.m4, Ebeam - 0.5, th, ke);
         double dKe = (le_ex(c.m1, c.m2, c.m3, c.m4, Ebeam, th, ke + 0.05) -
                       le_ex(c.m1, c.m2, c.m3, c.m4, Ebeam, th, ke - 0.05)) / 0.1;
         double dTh = (le_ex(c.m1, c.m2, c.m3, c.m4, Ebeam, th + 0.1 * TMath::DegToRad(), ke) -
                       le_ex(c.m1, c.m2, c.m3, c.m4, Ebeam, th - 0.1 * TMath::DegToRad(), ke)) / 0.2;
         if (std::isnan(dEb) || std::isnan(dKe) || std::isnan(dTh))
            continue;
         // the uncorrected arm carries the full vertex spread; the corrected one only straggling
         double sKEw, sThw, sKEb, sThb;
         bool exW, exB;
         lv_lookup(resW, ke, sKEw, sThw, exW);
         lv_lookup(resB, ke, sKEb, sThb, exB);
         double rW = std::sqrt(std::pow(dKe * sKEw, 2) + std::pow(dTh * sThw, 2) + std::pow(dEb * sEbUncorr, 2));
         double rB = std::sqrt(std::pow(dKe * sKEb, 2) + std::pow(dTh * sThb, 2) + std::pow(dEb * sEbCorr, 2));
         gLoc[i]->SetPoint(gLoc[i]->GetN(), thl, ke);
         gEb[i]->SetPoint(gEb[i]->GetN(), tcmA, std::fabs(dEb));
         gKe[i]->SetPoint(gKe[i]->GetN(), tcmA, std::fabs(dKe));
         gTh[i]->SetPoint(gTh[i]->GetN(), tcmA, std::fabs(dTh));
         gResW[i]->SetPoint(gResW[i]->GetN(), tcmA, rW);
         gResB[i]->SetPoint(gResB[i]->GetN(), tcmA, rB);
         if (((int)tcmA) % 20 == 0)
            printf("  %8.0f %9.1f %9.2f | %11.4f %10.3f %11.4f | %12.3f %12.3f %s\n", tcmA, thl, ke, dEb, dKe, dTh, rW,
                   rB, (thl > 90 ? "  <- theta_lab > 90, NOT covered by the (p,p') campaign" : (exW ? "  (KE outside the measured range)" : "")));
      }
   }

   auto *c1 = new TCanvas("lev", "lev", 1620, 1020);
   c1->Divide(2, 2);
   auto style = [&](TGraph *g, size_t i) {
      g->SetLineColor(ch[i].col);
      g->SetLineWidth(3);
      g->SetLineStyle(ch[i].sty);
   };
   auto panel = [&](int pad, std::vector<TGraph *> &g, const char *tit, const char *xt, const char *yt, double x0,
                    double x1, double y0, double y1, bool logy) {
      c1->cd(pad);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.13);
      if (logy) gPad->SetLogy();
      auto *f = new TH1F(Form("f%d", pad), Form(";%s;%s", xt, yt), 100, x0, x1);
      f->GetYaxis()->SetRangeUser(y0, y1);
      f->GetYaxis()->SetTitleOffset(1.25);
      f->GetYaxis()->SetTitleSize(0.045);
      f->GetXaxis()->SetTitleSize(0.045);
      f->SetLineColor(kWhite);
      f->Draw();
      auto *l = new TLegend(0.55, 0.70, 0.93, 0.88);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextSize(0.040);
      for (size_t i = 0; i < g.size(); ++i) {
         style(g[i], i);
         g[i]->Draw("l same");
         l->AddEntry(g[i], ch[i].name, "l");
      }
      l->Draw();
      auto *tx = new TLatex(0.17, 0.90, tit);
      tx->SetNDC();
      tx->SetTextSize(0.042);
      tx->Draw();
      return l;
   };

   panel(1, gLoc, "#bf{A}  where the ejectile goes", "#theta_{lab} [deg]", "ejectile KE [MeV]", 0, 180, 0.2, 60, kTRUE);
   {  // the (p,p') recoil cannot pass 90 deg; only a transfer ejectile reaches the backward half,
      // and that half is where the leverage in panels B and C lives
      auto *l90 = new TLine(90, 0.2, 90, 60);
      l90->SetLineStyle(2);
      l90->SetLineColor(kGray + 2);
      l90->Draw();
      auto *t90 = new TLatex(0.52, 0.20, "backward in the lab #rightarrow");
      t90->SetNDC();
      t90->SetTextSize(0.036);
      t90->SetTextColor(kGray + 3);
      t90->Draw();
   }
   panel(2, gEb, "#bf{B}  leverage of the beam energy", "#theta_{cm} [deg]", "|dE_{x}/dE_{beam}|", 0, 180, 1e-3, 2,
         kTRUE);
   panel(3, gKe, "#bf{C}  leverage of the ejectile energy", "#theta_{cm} [deg]", "|dE_{x}/dKE|", 0, 180, 0.02, 20,
         kTRUE);
   auto *legD = panel(4, gResW, "#bf{D}  implied #sigma(E_{x})", "#theta_{cm} [deg]", "#sigma(E_{x}) [MeV]", 0, 180,
                      0.005, 20, kTRUE);
   legD->SetX1NDC(0.60); legD->SetX2NDC(0.95);
   legD->SetY1NDC(0.17); legD->SetY2NDC(0.35);
   c1->cd(4);
   {  // where a transfer angular distribution actually has its yield
      auto *b = new TBox(8, 0.005, 45, 20);
      b->SetFillColorAlpha(kGray + 1, 0.16);
      b->SetLineWidth(0);
      b->Draw();
      auto *tb = new TLatex(0.20, 0.24, "#splitline{transfer}{peaks here}");
      tb->SetNDC();
      tb->SetTextSize(0.034);
      tb->SetTextColor(kGray + 3);
      tb->Draw();
      gPad->RedrawAxis();
   }
   for (size_t i = 0; i < gResB.size(); ++i) {
      gResB[i]->SetLineColor(ch[i].col);
      gResB[i]->SetLineWidth(3);
      gResB[i]->SetLineStyle(2);
      gResB[i]->Draw("l same");
   }
   auto *tx2 = new TLatex(0.17, 0.848, "solid: today.  dashed: 4 T + 2 mm, E_{beam}(z)");
   tx2->SetNDC();
   tx2->SetTextSize(0.038);
   tx2->Draw();
   auto *tx3 = new TLatex(0.17, 0.805, "steps are the KE bins of the measured lookup");
   tx3->SetNDC();
   tx3->SetTextSize(0.032);
   tx3->SetTextColor(kGray + 3);
   tx3->Draw();

   c1->SaveAs(outDir + "channel_leverage.png");
   printf("\n  wrote %schannel_leverage.png\n\n", outDir.Data());
}
