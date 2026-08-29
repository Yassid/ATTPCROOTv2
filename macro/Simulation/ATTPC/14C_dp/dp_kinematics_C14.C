/// @file dp_kinematics_C14.C
/// @brief The 14C(d,p)15C kinematic plane as each configuration of the matrix reconstructs it.
///
///   root -b -q 'dp_kinematics_C14.C()'
///
/// One panel per configuration: the reconstructed proton (KE, theta_lab) of both simulated levels
/// together, with the two-body loci drawn over them. This is the first thing to look at in any of
/// these simulations -- before any resolution number, the question is whether the reconstructed
/// points sit on the curve the reaction says they must.
///
/// THE LOCI ARE PARAMETRIC IN theta_cm, not solved for KE at fixed angle. In inverse kinematics a
/// transfer locus can be double-valued in theta_lab -- two centre-of-mass angles landing at the
/// same lab angle with different energies -- and a bisection in KE silently returns one of them.
/// Stepping theta_cm and computing (theta_lab, KE) draws the whole curve including any turn-back.
///
/// THREE CURVES PER LEVEL, not one. The beam loses ~10 MeV crossing the chamber and the vertex is
/// uniform in z, so a single locus is wrong for every event except those at the mean vertex. The
/// solid line is the mid-chamber beam energy and the dashed pair are the entrance and far-end
/// values: the band between them is the kinematic spread the constant-beam-energy analysis has to
/// live with, drawn to scale.

#include <vector>

static const int NC = 6;
static const char *CFG[NC] = {"b285_attpc", "b285_2mm", "b400_attpc", "b400_2mm", "b700_attpc", "b700_2mm"};
static const char *CLB[NC] = {"2.85 T, AT-TPC pads", "2.85 T, 2 mm pads", "4 T, AT-TPC pads",
                              "4 T, 2 mm pads",      "7 T, AT-TPC pads", "7 T, 2 mm pads"};

static const double U = 931.49401;
static const double M1 = 14.003242 * U;   // 14C beam
static const double M2 = 2.0141018 * U;   // d target
static const double M3 = 1.007825 * U;    // p ejectile
static const double M4 = 15.0105993 * U;  // 15C residual, ground state

static double kn_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
/// the ejectile (KE, theta_lab) at a given centre-of-mass angle of the EJECTILE
static bool kn_forward(double m4, double K, double thcm_deg, double &ke, double &thlab)
{
   double E1 = K + M1;
   double s = M1 * M1 + M2 * M2 + 2 * M2 * E1;
   double rs = std::sqrt(s);
   if (rs < M3 + m4)
      return false;
   double pcm = kn_om2(s, M3 * M3, m4 * m4) / (2 * rs);
   double E3cm = std::sqrt(pcm * pcm + M3 * M3);
   double plab1 = std::sqrt(E1 * E1 - M1 * M1);
   double beta = plab1 / (E1 + M2), gamma = 1.0 / std::sqrt(1 - beta * beta);
   double th = thcm_deg * TMath::DegToRad();
   double pz = gamma * (pcm * std::cos(th) + beta * E3cm);
   double pt = pcm * std::sin(th);
   double E3 = gamma * (E3cm + beta * pcm * std::cos(th));
   ke = E3 - M3;
   thlab = std::atan2(pt, pz) * TMath::RadToDeg();
   return ke > 0;
}

/// @param eLo/eMid/eHi  beam energy at the far end / mid-chamber / entrance [MeV]
/// @param thLo..keHi    axis window. The default is the whole plane; pass (80,180,0,15) for the
///                      BACKWARD region, where the transfer cross section actually has its yield
///                      and where the full-range view compresses everything into one corner.
/// @param dTh, dKe      bin widths, degrees and MeV. Finer bins turn the plane into something
///                      closer to a scatter plot -- with ~5000 reconstructed protons over the full
///                      range, 1 deg x 0.2 MeV already puts most bins at 0 or 1 count, which shows
///                      individual tracks rather than a density. That is a legitimate view, but do
///                      not read a colour scale off it.
/// @param suffix        appended to the output file name, so the two views do not overwrite
void dp_kinematics_C14(TString root = "/mnt/f/a1954_C14dp_hf", Double_t eLo = 150.8, Double_t eMid = 155.9,
                       Double_t eHi = 161.0, TString outDir = "", Double_t thLo = 0, Double_t thHi = 180,
                       Double_t keLo = 0, Double_t keHi = 60, TString suffix = "", Double_t dTh = 1.0,
                       Double_t dKe = 0.2)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   const int NL = 2;
   const char *LVL[NL] = {"gs", "ex0740"};
   const double LEX[NL] = {0.0, 0.740};
   const int LCOL[NL] = {kBlack, kRed + 1};

   // the loci: [level][beam energy], parametric in theta_cm
   std::vector<TGraph *> loci;
   std::vector<int> lcol, lsty;
   const double eb[3] = {eHi, eMid, eLo};
   for (int l = 0; l < NL; ++l)
      for (int e = 0; e < 3; ++e) {
         auto *g = new TGraph();
         for (double tcm = 1; tcm <= 179; tcm += 0.5) {
            double ke, thl;
            if (!kn_forward(M4 + LEX[l], eb[e], tcm, ke, thl))
               continue;
            if (ke <= keLo || ke > keHi)
               continue;
            g->SetPoint(g->GetN(), thl, ke);
         }
         loci.push_back(g);
         lcol.push_back(LCOL[l]);
         lsty.push_back(e == 1 ? 1 : 2); // mid-chamber solid, the two ends dashed
      }

   auto *cv = new TCanvas("dpkin" + suffix, "dpkin", 1700, 1000);
   cv->Divide(3, 2);
   for (int c = 0; c < NC; ++c) {
      cv->cd(c + 1);
      gPad->SetLogz();
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.13);
      gPad->SetBottomMargin(0.13);
      auto *h = new TH2D(Form("hk%d%s", c, suffix.Data()), TString(CLB[c]) + ";#theta_{lab} [deg];proton KE [MeV]",
                         (int)std::max(10.0, (thHi - thLo) / dTh), thLo, thHi,
                         (int)std::max(10.0, (keHi - keLo) / dKe), keLo, keHi);
      h->GetXaxis()->SetTitleSize(0.045);
      h->GetYaxis()->SetTitleSize(0.045);
      h->GetYaxis()->SetTitleOffset(1.25);
      long n = 0, nIn = 0;
      for (int l = 0; l < NL; ++l) {
         TString f = gSystem->GetFromPipe(TString::Format("ls %s/%s/exres_%s_s*_%s.root 2>/dev/null | head -1",
                                                          root.Data(), CFG[c], LVL[l], CFG[c]));
         f = f.Strip(TString::kBoth);
         if (f.IsNull())
            continue;
         TFile *fr = TFile::Open(f);
         if (!fr || fr->IsZombie())
            continue;
         TTree *t = (TTree *)fr->Get("res");
         if (!t) { fr->Close(); continue; }
         double keReco, thReco;
         t->SetBranchAddress("keReco", &keReco);
         t->SetBranchAddress("thReco", &thReco);
         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            h->Fill(thReco, keReco);
            ++n;
            if (thReco >= thLo && thReco < thHi && keReco >= keLo && keReco < keHi)
               ++nIn;
         }
         fr->Close();
      }
      h->Draw("colz");
      for (size_t k = 0; k < loci.size(); ++k) {
         loci[k]->SetLineColor(lcol[k]);
         loci[k]->SetLineStyle(lsty[k]);
         loci[k]->SetLineWidth(lsty[k] == 1 ? 3 : 1);
         loci[k]->Draw("l same");
      }
      auto *tx = new TLatex(0.17, 0.86,
                            (thLo > 0 || thHi < 180 || keHi < 60)
                               ? TString::Format("%ld of %ld protons in view", nIn, n)
                               : TString::Format("%ld reconstructed protons", n));
      tx->SetNDC();
      tx->SetTextSize(0.040);
      tx->Draw();
      if (c == 0) {
         auto *lg = new TLegend(0.45, 0.62, 0.86, 0.86);
         lg->SetBorderSize(0);
         lg->SetFillStyle(0);
         lg->SetTextSize(0.038);
         lg->AddEntry(loci[1], "^{15}C g.s.", "l");
         lg->AddEntry(loci[4], "^{15}C 0.740 MeV", "l");
         lg->AddEntry(loci[0], "beam at entrance / far end", "l");
         lg->Draw();
      }
   }
   cv->SaveAs(outDir + "dp_kinematics" + suffix + ".png");
   printf("\n  wrote %sdp_kinematics%s.png\n", outDir.Data(), suffix.Data());
   printf("  loci drawn at E_beam = %.1f (entrance), %.1f (mid), %.1f (far end) MeV\n", eHi, eMid, eLo);
   printf("  bins %.2f deg x %.2f MeV\n\n", dTh, dKe);
}
