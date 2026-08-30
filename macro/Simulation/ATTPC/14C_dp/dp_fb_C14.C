/// @file dp_fb_C14.C
/// @brief 14C(d,p)15C kinematics and excitation energy, split FORWARD / BACKWARD in the lab.
///
///   root -b -q 'dp_fb_C14.C("/mnt/f/a1954_C14dp")'
///
/// The split is the point of this channel. A forward proton is fast (20-50 MeV) on a nearly
/// straight trajectory whose curvature is barely measured; a backward one is slow (2-8 MeV) and
/// spirals for metres, so its momentum is over-determined. They are two different measurements
/// sharing a detector, and averaging them hides both.
///
/// Writes three figures:
///   dp_fb_kinematics.png  the (theta_lab, KE) plane per configuration, with the two-body loci of
///                         every simulated level and the 90 deg divider
///   dp_fb_excitation.png  E_x per configuration, forward and backward overlaid, each normalised
///                         to itself so the WIDTHS compare rather than the yields
///   dp_fb_summary.png     sigma(E_x) and acceptance, forward vs backward, across the matrix
///
/// Missing samples are skipped and reported, so this can be run while a campaign is still going.

#include <algorithm>
#include <map>
#include <vector>

static const double U = 931.49401;
static const double M1 = 14.003242 * U, M2 = 2.0141018 * U, M3 = 1.007825 * U, M4 = 15.0105993 * U;
static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const char *CLB[NC] = {"2.85 T, AT-TPC", "2.85 T, 2 mm", "4 T, AT-TPC",
                              "4 T, 2 mm",      "7 T, AT-TPC", "7 T, 2 mm"};
static const int NL = 3;
static const char *LVL[NL] = {"gs", "ex0740", "ex3103"};
static const double LEX[NL] = {0.0, 0.740, 3.103};
static const int LCOL[NL] = {kBlack, kAzure + 2, kGreen + 3};

static double fb_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static bool fb_fwd(double m4, double K, double tcm, double &ke, double &thl)
{
   double E1 = K + M1, s = M1 * M1 + M2 * M2 + 2 * M2 * E1, rs = std::sqrt(s);
   if (rs < M3 + m4) return false;
   double pcm = fb_om2(s, M3 * M3, m4 * m4) / (2 * rs), E3cm = std::sqrt(pcm * pcm + M3 * M3);
   double p1 = std::sqrt(E1 * E1 - M1 * M1), beta = p1 / (E1 + M2), g = 1.0 / std::sqrt(1 - beta * beta);
   double th = tcm * TMath::DegToRad();
   ke = g * (E3cm + beta * pcm * std::cos(th)) - M3;
   thl = std::atan2(pcm * std::sin(th), g * (pcm * std::cos(th) + beta * E3cm)) * TMath::RadToDeg();
   return ke > 0;
}
static double fb_q(std::vector<double> v, double p)
{
   if (v.size() < 15) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}
static TString fb_find(const TString &d, const char *cfg, const char *lvl, const char *pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", d.Data(), cfg, pre, lvl, cfg));
   return f.Strip(TString::kBoth);
}

/// @param thSplit  the dividing lab angle. 90 deg is the physical forward/backward boundary.
void dp_fb_C14(TString root = "/mnt/f/a1954_C14dp", Double_t thSplit = 90.0, Double_t Ebeam = 155.9,
               TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   // [cfg][level][0=fwd,1=back]
   std::map<TString, std::map<TString, std::array<TH1D *, 2>>> hEx;
   std::map<TString, TH2D *> hKin;
   std::map<TString, std::array<std::vector<double>, 2>> resid; // Ex residual, per region, all levels
   std::map<TString, std::array<long, 2>> nTrk;

   for (int c = 0; c < NC; ++c) {
      hKin[CFG[c]] = new TH2D(Form("hk_%d", c), TString(CLB[c]) + ";#theta_{lab} [deg];proton KE [MeV]", 180, 0, 180,
                              240, 0, 60);
      for (int l = 0; l < NL; ++l)
         for (int r = 0; r < 2; ++r)
            hEx[CFG[c]][LVL[l]][r] =
               new TH1D(Form("he_%d_%d_%d", c, l, r), "", 200, -2.5, 6.0);
      nTrk[CFG[c]] = {0, 0};
      for (int l = 0; l < NL; ++l) {
         TString f = fb_find(root, CFG[c], LVL[l], "exres");
         if (f.IsNull()) { printf("  missing: %-12s %s\n", CFG[c], LVL[l]); continue; }
         TFile *fr = TFile::Open(f);
         if (!fr || fr->IsZombie()) continue;
         TTree *t = (TTree *)fr->Get("res");
         if (!t) { fr->Close(); continue; }
         double exReco, thTrue, thReco, keReco;
         t->SetBranchAddress("exReco", &exReco);
         t->SetBranchAddress("thTrue", &thTrue);
         t->SetBranchAddress("thReco", &thReco);
         t->SetBranchAddress("keReco", &keReco);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            int r = (thTrue >= thSplit) ? 1 : 0;
            hEx[CFG[c]][LVL[l]][r]->Fill(exReco);
            hKin[CFG[c]]->Fill(thReco, keReco);
            resid[CFG[c]][r].push_back(exReco - LEX[l]);
            ++nTrk[CFG[c]][r];
         }
         fr->Close();
      }
   }

   // ---- figure 1: the kinematic plane per configuration --------------------------------------
   {
      auto *cv = new TCanvas("fbk", "fbk", 1700, 1000);
      cv->Divide(3, 2);
      std::vector<TGraph *> loci;
      for (int l = 0; l < NL; ++l) {
         auto *g = new TGraph();
         for (double tcm = 1; tcm <= 179; tcm += 0.4) {
            double ke, thl;
            if (fb_fwd(M4 + LEX[l], Ebeam, tcm, ke, thl) && ke < 60) g->SetPoint(g->GetN(), thl, ke);
         }
         g->SetLineColor(LCOL[l]);
         g->SetLineWidth(2);
         loci.push_back(g);
      }
      for (int c = 0; c < NC; ++c) {
         cv->cd(c + 1);
         gPad->SetLogz();
         gPad->SetLeftMargin(0.13);
         hKin[CFG[c]]->Draw("colz");
         for (auto *g : loci) g->Draw("l same");
         auto *l = new TLine(thSplit, 0, thSplit, 60);
         l->SetLineStyle(2); l->SetLineColor(kGray + 2); l->Draw();
         auto *tf = new TLatex(0.17, 0.30, TString::Format("forward: %ld", nTrk[CFG[c]][0]));
         tf->SetNDC(); tf->SetTextSize(0.040); tf->Draw();
         auto *tb = new TLatex(0.17, 0.24, TString::Format("backward: %ld", nTrk[CFG[c]][1]));
         tb->SetNDC(); tb->SetTextSize(0.040); tb->SetTextColor(kRed + 1); tb->Draw();
      }
      cv->SaveAs(outDir + "dp_fb_kinematics.png");
   }

   // ---- figure 2: E_x, forward vs backward, per configuration ---------------------------------
   {
      auto *cv = new TCanvas("fbe", "fbe", 1700, 1000);
      cv->Divide(3, 2);
      printf("\n  sigma(E_x) [MeV], IQR/1.349 over all simulated levels\n");
      printf("  %-14s %10s %8s | %10s %8s\n", "config", "forward", "n", "backward", "n");
      for (int c = 0; c < NC; ++c) {
         cv->cd(c + 1);
         gPad->SetLeftMargin(0.13);
         // sum the levels, then normalise each region to itself: the comparison is of WIDTH
         TH1D *sum[2] = {nullptr, nullptr};
         for (int r = 0; r < 2; ++r)
            for (int l = 0; l < NL; ++l) {
               if (!hEx[CFG[c]][LVL[l]][r]) continue;
               if (!sum[r]) { sum[r] = (TH1D *)hEx[CFG[c]][LVL[l]][r]->Clone(Form("s_%d_%d", c, r));
                              sum[r]->SetDirectory(nullptr); }
               else sum[r]->Add(hEx[CFG[c]][LVL[l]][r]);
            }
         double ymax = 0;
         for (int r = 0; r < 2; ++r)
            if (sum[r] && sum[r]->Integral() > 0) { sum[r]->Scale(1.0 / sum[r]->Integral()); ymax = std::max(ymax, sum[r]->GetMaximum()); }
         bool first = true;
         const int col[2] = {kGray + 3, kRed + 1};
         for (int r = 0; r < 2; ++r) {
            if (!sum[r] || sum[r]->Integral() <= 0) continue;
            sum[r]->SetLineColor(col[r]);
            sum[r]->SetLineWidth(2);
            sum[r]->SetTitle(TString(CLB[c]) + ";E_{x} [MeV];normalised");
            sum[r]->GetYaxis()->SetRangeUser(0, 1.25 * ymax);
            sum[r]->Draw(first ? "hist" : "hist same");
            first = false;
         }
         for (int l = 0; l < NL; ++l) {
            auto *ln = new TLine(LEX[l], 0, LEX[l], 1.25 * ymax);
            ln->SetLineStyle(2); ln->SetLineColor(kGray + 1); ln->Draw();
         }
         if (c == 0 && !first) {
            auto *lg = new TLegend(0.55, 0.70, 0.90, 0.87);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
            if (sum[0]) lg->AddEntry(sum[0], "forward (#theta_{lab} < 90#circ)", "l");
            if (sum[1]) lg->AddEntry(sum[1], "backward (> 90#circ)", "l");
            lg->Draw();
         }
         double sf = (fb_q(resid[CFG[c]][0], .75) - fb_q(resid[CFG[c]][0], .25)) / 1.349;
         double sb = (fb_q(resid[CFG[c]][1], .75) - fb_q(resid[CFG[c]][1], .25)) / 1.349;
         printf("  %-14s %10.3f %8ld | %10.3f %8ld\n", CFG[c], sf, nTrk[CFG[c]][0], sb, nTrk[CFG[c]][1]);
         auto *tx = new TLatex(0.17, 0.84, TString::Format("#sigma  fwd %.3f   back %.3f MeV", sf, sb));
         tx->SetNDC(); tx->SetTextSize(0.040); tx->Draw();
      }
      cv->SaveAs(outDir + "dp_fb_excitation.png");
   }

   // ---- figure 3: the summary across the matrix ------------------------------------------------
   {
      auto *cv = new TCanvas("fbs", "fbs", 900, 620);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.16);
      gPad->SetGridy();
      gPad->SetLogy();
      auto *fr = new TH1F("fbsf", ";field [T] / pad pitch [mm];#sigma(E_{x}) [MeV]", NC, 0, NC);
      const char *sh[NC] = {"2.85 / 8#times12", "2.85 / 2mm", "4 / 8#times12", "4 / 2mm", "7 / 8#times12", "7 / 2mm"};
      for (int c = 0; c < NC; ++c) fr->GetXaxis()->SetBinLabel(c + 1, sh[c]);
      fr->GetYaxis()->SetRangeUser(0.01, 5);
      fr->SetLineColor(kWhite);
      fr->Draw();
      auto *gf = new TGraph(), *gb = new TGraph();
      for (int c = 0; c < NC; ++c) {
         double sf = (fb_q(resid[CFG[c]][0], .75) - fb_q(resid[CFG[c]][0], .25)) / 1.349;
         double sb = (fb_q(resid[CFG[c]][1], .75) - fb_q(resid[CFG[c]][1], .25)) / 1.349;
         if (!std::isnan(sf)) gf->SetPoint(gf->GetN(), c + 0.5, sf);
         if (!std::isnan(sb)) gb->SetPoint(gb->GetN(), c + 0.5, sb);
      }
      gf->SetMarkerStyle(20); gf->SetMarkerColor(kGray + 3); gf->SetLineColor(kGray + 3); gf->SetLineWidth(3);
      gb->SetMarkerStyle(21); gb->SetMarkerColor(kRed + 1); gb->SetLineColor(kRed + 1); gb->SetLineWidth(3);
      gf->SetMarkerSize(1.6); gb->SetMarkerSize(1.6);
      gf->Draw("pl same"); gb->Draw("pl same");
      auto *lg = new TLegend(0.55, 0.74, 0.92, 0.88);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
      lg->AddEntry(gf, "forward (#theta_{lab} < 90#circ)", "pl");
      lg->AddEntry(gb, "backward (> 90#circ)", "pl");
      lg->Draw();
      cv->SaveAs(outDir + "dp_fb_summary.png");
   }
   printf("\n  wrote %sdp_fb_{kinematics,excitation,summary}.png\n\n", outDir.Data());
}
