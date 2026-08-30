/// @file dp_perbin_C14.C
/// @brief Tracking efficiency, the implied measurement error, and the resolution -- all against
/// the TRUE lab angle in 10 degree bins, for every configuration on one set of axes.
///
///   root -b -q 'dp_perbin_C14.C("/mnt/f/a1954_C14dp")'
///
///   A  efficiency      truth-matched good fits / truth protons in the bin
///   B  implied sigma   measSigma * sqrt(chi2/ndf): what the hits are ACTUALLY scattered by. If it
///                      differs from the measSigma the fit was given, the quality metric is
///                      mis-scaled there -- and it varies with pad pitch, field and angle, so one
///                      global value cannot be right everywhere.
///   C  resolution      sigma of (KE_reco - KE_true)/KE_true
///
/// Binning in the TRUE angle throughout, so a track that is reconstructed at the wrong angle is
/// counted where it came from and shows up as an efficiency loss rather than migrating.

#include <algorithm>
#include <vector>

static const double U = 931.49401, MP = 1.007825 * U;
static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const char *CLB[NC] = {"2.85 T, AT-TPC", "2.85 T, 2 mm", "4 T, AT-TPC",
                              "4 T, 2 mm",      "7 T, AT-TPC", "7 T, 2 mm"};
static const int COL[NC] = {kGray + 3, kGray + 1, kAzure + 2, kAzure - 4, kRed + 1, kOrange + 7};
static const int MRK[NC] = {20, 24, 21, 25, 22, 26};
/// the measSigma each configuration was FITTED with (0.6 for the AT-TPC planes, 0.35 for 2 mm)
static const double MSUSED[NC] = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6};

static double pb_q(std::vector<double> v, double p)
{
   if (v.size() < 12) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_perbin_C14(TString root = "/mnt/f/a1954_C14dp", TString level = "gs", TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   const int NB = 18; // 10 degree bins, 0-180
   std::vector<TGraph *> gEff, gSig, gRes;

   for (int c = 0; c < NC; ++c) {
      TString fe = gSystem->GetFromPipe(
         TString::Format("ls %s/%s/exres_%s_s*_%s.root 2>/dev/null | head -1", root.Data(), CFG[c], level.Data(), CFG[c]));
      TString fg = gSystem->GetFromPipe(
         TString::Format("ls %s/%s/%s_s*_%s_genfit.root 2>/dev/null | head -1", root.Data(), CFG[c], level.Data(), CFG[c]));
      TString fs = gSystem->GetFromPipe(
         TString::Format("ls %s/sims_%.4s/%s_s*_sim.root 2>/dev/null | head -1", root.Data(), CFG[c], level.Data()));
      fe = fe.Strip(TString::kBoth); fg = fg.Strip(TString::kBoth); fs = fs.Strip(TString::kBoth);
      auto *ge = new TGraph(), *gs = new TGraph(), *gr = new TGraph();

      // denominator: truth protons per bin
      long nGen[NB] = {0};
      if (!fs.IsNull()) {
         TFile *f = TFile::Open(fs);
         TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
         TClonesArray *mc = nullptr;
         if (t) {
            t->SetBranchAddress("MCTrack", &mc);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               if (!mc) continue;
               for (int k = 0; k < mc->GetEntriesFast(); ++k) {
                  auto *p = (AtMCTrack *)mc->At(k);
                  if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1) continue;
                  double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
                  double pm = std::sqrt(px * px + py * py + pz * pz);
                  if (pm <= 0) continue;
                  int b = (int)(std::acos(pz / pm) * TMath::RadToDeg() / 10.0);
                  if (b >= 0 && b < NB) ++nGen[b];
                  break;
               }
            }
         }
         if (f) f->Close();
      }
      // numerator + resolution, from the exres tree (already truth-matched)
      long nOK[NB] = {0};
      std::vector<double> res[NB];
      if (!fe.IsNull()) {
         TFile *f = TFile::Open(fe);
         TTree *t = f ? (TTree *)f->Get("res") : nullptr;
         if (t) {
            double thTrue, keTrue, keReco;
            t->SetBranchAddress("thTrue", &thTrue);
            t->SetBranchAddress("keTrue", &keTrue);
            t->SetBranchAddress("keReco", &keReco);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               int b = (int)(thTrue / 10.0);
               if (b < 0 || b >= NB || keTrue <= 0) continue;
               ++nOK[b];
               res[b].push_back(100.0 * (keReco - keTrue) / keTrue);
            }
         }
         if (f) f->Close();
      }
      // implied residual, from chi2/ndf of every fit in the bin (no truth match: this is a
      // property of the FIT, not of the selection)
      std::vector<double> c2[NB];
      if (!fg.IsNull()) {
         TFile *f = TFile::Open(fg);
         TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
         TClonesArray *te = nullptr;
         if (t) {
            t->SetBranchAddress("AtTrackingEvent", &te);
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               if (!te || te->GetEntriesFast() == 0) continue;
               auto *ev = (AtTrackingEvent *)te->At(0);
               if (!ev) continue;
               for (auto &ft : ev->GetFittedTracks()) {
                  if (!ft) continue;
                  const auto &md = ft->GetTrackMetadata();
                  double ndf = md ? md->GetNdf() : 0, ch = md ? md->GetChi2() : 0;
                  if (ndf <= 0) continue;
                  int b = (int)(ft->GetKinematicsXtr().theta * TMath::RadToDeg() / 10.0);
                  if (b >= 0 && b < NB) c2[b].push_back(ch / ndf);
                  break;
               }
            }
         }
         if (f) f->Close();
      }
      for (int b = 0; b < NB; ++b) {
         double x = 10 * b + 5;
         if (nGen[b] >= 20) ge->SetPoint(ge->GetN(), x, (double)nOK[b] / nGen[b]);
         double m = pb_q(c2[b], .5);
         if (!std::isnan(m)) gs->SetPoint(gs->GetN(), x, MSUSED[c] * std::sqrt(m));
         double s = (pb_q(res[b], .75) - pb_q(res[b], .25)) / 1.349;
         if (!std::isnan(s)) gr->SetPoint(gr->GetN(), x, s);
      }
      gEff.push_back(ge); gSig.push_back(gs); gRes.push_back(gr);
   }

   auto *cv = new TCanvas("pb", "pb", 1650, 560);
   cv->Divide(3, 1);
   const char *yt[3] = {"tracking efficiency", "implied hit #sigma [mm]", "#sigma(KE)/KE [%]"};
   const char *tt[3] = {"A  efficiency", "B  measSigma implied by the fit", "C  energy resolution"};
   const double y0[3] = {0, 0.05, 0.05}, y1[3] = {1.05, 20, 60};
   const bool lg[3] = {false, true, true};
   for (int p = 0; p < 3; ++p) {
      cv->cd(p + 1);
      gPad->SetLeftMargin(0.14);
      gPad->SetBottomMargin(0.14);
      gPad->SetGridy();
      if (lg[p]) gPad->SetLogy();
      auto *fr = new TH1F(Form("f%d", p), TString(";#theta_{lab} true [deg];") + yt[p], 18, 0, 180);
      fr->GetYaxis()->SetRangeUser(y0[p], y1[p]);
      fr->GetYaxis()->SetTitleSize(0.048);
      fr->GetYaxis()->SetTitleOffset(1.35);
      fr->GetXaxis()->SetTitleSize(0.048);
      fr->SetLineColor(kWhite);
      fr->Draw();
      for (int c = 0; c < NC; ++c) {
         auto *g = (p == 0) ? gEff[c] : (p == 1 ? gSig[c] : gRes[c]);
         g->SetLineColor(COL[c]); g->SetMarkerColor(COL[c]); g->SetMarkerStyle(MRK[c]);
         g->SetMarkerSize(1.2); g->SetLineWidth(2);
         g->Draw("pl same");
      }
      auto *l90 = new TLine(90, y0[p], 90, y1[p]);
      l90->SetLineStyle(2); l90->SetLineColor(kGray + 2); l90->Draw();
      if (p == 1) { // the value actually handed to the fit
         auto *lms = new TLine(0, 0.6, 180, 0.6);
         lms->SetLineStyle(3); lms->SetLineColor(kGreen + 3); lms->SetLineWidth(2); lms->Draw();
         auto *t = new TLatex(0.17, 0.30, "green: the measSigma the fit was given (0.6 mm)");
         t->SetNDC(); t->SetTextSize(0.036); t->SetTextColor(kGreen + 3); t->Draw();
      }
      if (p == 0) {
         auto *leg = new TLegend(0.17, 0.16, 0.62, 0.44);
         leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
         for (int c = 0; c < NC; ++c) leg->AddEntry(gEff[c], CLB[c], "pl");
         leg->Draw();
      }
      auto *tx = new TLatex(0.17, 0.91, tt[p]);
      tx->SetNDC(); tx->SetTextSize(0.048); tx->Draw();
   }
   cv->SaveAs(outDir + "dp_perbin_" + level + ".png");
   printf("\n  wrote %sdp_perbin_%s.png\n\n", outDir.Data(), level.Data());
}
