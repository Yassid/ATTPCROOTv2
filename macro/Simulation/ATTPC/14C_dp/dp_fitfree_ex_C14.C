/// @file dp_fitfree_ex_C14.C
/// @brief Excitation energy of 15C from hit geometry alone -- no helix fit, no GENFIT.
///
///   root -b -q 'dp_fitfree_ex_C14.C("/mnt/f/a1954_C14dp_sm","/mnt/f/a1954_C14dp/sims_b285",2.85,"2p85T")'
///
/// This is the payoff test for the SpecMAT/HELIOS route. dp_specmat_C14.C and dp_helios_C14.C
/// showed that the two geometric observables close on truth; the question that decides whether
/// the route is USEFUL is whether the resulting excitation energy is good enough to separate
/// states -- and in particular whether it recovers the 7 T column, where GENFIT diverges
/// (chi2/ndf ~1300 backward, only ~21 % of tracks passing the quality cut).
///
/// The chain is: d_max -> R; the apex identity psi = pi - 2*asin(d0/2R) with the z advance to the
/// apex -> theta_lab; then p_T = 0.3*B*R and p = p_T/sin(theta) -> KE; then the same two-body
/// inversion used everywhere else in this analysis -> Ex. Every input is a hit position.
///
/// The benchmark to compare against is the GENFIT result on the same channel, sigma(Ex) backward:
/// 0.197 MeV at 2.85 T, 0.241 at 4 T, and at 7 T a number that cannot be quoted because only
/// 625 of ~3000 tracks survive.
///
/// CAVEAT STATED UP FRONT: the beam energy is held constant here. The campaign showed that
/// evaluating it at the reconstructed vertex takes backward sigma(Ex) from 0.178 to 0.064 MeV,
/// so the numbers below are the CONSTANT-beam-energy figure and are directly comparable to the
/// 0.197 / 0.241 above, not to the vertex-corrected ones.

#include <algorithm>
#include <vector>

static double ffx_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics, byte for byte the one in ex_res_C14_hf.C / acceptance_C14.C / pp/ex_C14.C
static std::tuple<double, double> ffx_kine(double m1, double m2, double m3, double m4, double K_proj,
                                           double thetalab, double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * ffx_om2(s, m1 * m1, m2 * m2) * ffx_om2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (ffx_om2(s, m1 * m1, m2 * m2) * ffx_om2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

static double ffx_q(std::vector<double> v, double p)
{
   if (v.size() < 5)
      return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_fitfree_ex_C14(TString smDir = "/mnt/f/a1954_C14dp_sm", TString simDir = "/mnt/f/a1954_C14dp/sims_b285",
                       Double_t bField = 2.85, TString tag = "2p85T", Double_t Ebeam = 155.9,
                       Double_t thMin = 92.0, Int_t minHits = 120, Double_t activeR = 250.0, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   const double U = 931.49401;
   const double M1 = 14.003242 * U, M2 = 2.0141018 * U, M3 = 1.007825 * U, M4 = 15.0105993 * U;
   const double mp = M3, ZPAD = 1000.0, R2D = 57.29577951;
   const char *levs[3] = {"gs", "ex0740", "ex3103"};
   const double levEx[3] = {0.0, 0.740, 3.103};
   const int col[3] = {kAzure + 2, kOrange + 7, kGreen + 3};

   std::vector<double> exG[3], exT[3], cmG[3];

   for (int L = 0; L < 3; ++L) {
      TString rf = smDir + "/" + levs[L] + "/reco.root";
      TString sf = gSystem->GetFromPipe("ls " + simDir + "/" + TString(levs[L]) + "_*_sim.root 2>/dev/null | head -1");
      sf = sf.Strip(TString::kBoth);
      TFile *fr = TFile::Open(rf);
      TFile *fs = sf.IsNull() ? nullptr : TFile::Open(sf);
      if (!fr || fr->IsZombie() || !fs || fs->IsZombie()) {
         printf("  MISSING %s or %s\n", rf.Data(), sf.Data());
         continue;
      }
      TTree *tr = (TTree *)fr->Get("cbmsim");
      TTree *ts = (TTree *)fs->Get("cbmsim");
      TClonesArray *pe = nullptr, *mc = nullptr;
      tr->SetBranchAddress("AtPatternEvent", &pe);
      ts->SetBranchAddress("MCTrack", &mc);

      Long64_t N = std::min(tr->GetEntries(), ts->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         if (!mc)
            continue;
         double thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pm = std::sqrt(px * px + py * py + pz * pz);
            if (pm <= 0)
               break;
            thT = std::acos(pz / pm) * R2D;
            break;
         }
         if (thT < thMin)
            continue;

         tr->GetEntry(i);
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPatternEvent *)pe->At(0);
         if (!ev)
            continue;
         AtTrack *best = nullptr;
         size_t nb = 0;
         for (auto &t : ev->GetTrackCand())
            if (t.GetHitArray().size() > nb) {
               nb = t.GetHitArray().size();
               best = const_cast<AtTrack *>(&t);
            }
         if (!best || (int)nb < minHits)
            continue;
         auto *cl = best->GetHitClusterArray();
         if (!cl || cl->size() < 12)
            continue;
         std::vector<double> X, Y, Z;
         for (auto &c : *cl) {
            auto p = c.GetPosition();
            X.push_back(p.X());
            Y.push_back(p.Y());
            Z.push_back(ZPAD - p.Z());
         }
         if (std::hypot(X.back(), Y.back()) < std::hypot(X.front(), Y.front())) {
            std::reverse(X.begin(), X.end());
            std::reverse(Y.begin(), Y.end());
            std::reverse(Z.begin(), Z.end());
         }
         size_t iApex = 0;
         double dmax = -1;
         for (size_t k = 0; k < X.size(); ++k) {
            double d = std::hypot(X[k], Y[k]);
            if (d > dmax) { dmax = d; iApex = k; }
         }
         // QUALITY CUTS. a1975 established that this method needs them and I ran without any:
         // the apex fraction must sit strictly inside the track. If the apex is the LAST cluster
         // the track was truncated -- by the chamber wall or by the finder -- and d_max is then a
         // lower limit, not the diameter. That biases R low, hence KE low, hence Ex HIGH, and it
         // hits the g.s. hardest because its proton is the most energetic and has the largest
         // radius. Exactly the signature seen without the cut: sigma 1.084 (g.s.) vs 0.504
         // (3.103 MeV), all levels offset +0.45 MeV.
         double apexFrac = (double)iApex / (double)X.size();
         if (apexFrac < 0.10 || apexFrac > 0.75)
            continue;
         // and the apex must not sit against the active radius
         if (dmax > 0.95 * activeR)
            continue;
         double d0 = std::hypot(X[0], Y[0]);
         double R = 0.5 * dmax;
         if (!(R > 0) || d0 >= 2 * R)
            continue;
         double psi = TMath::Pi() - 2.0 * std::asin(d0 / (2 * R));
         double dz = Z[iApex] - Z[0];
         if (std::fabs(dz) < 1e-6)
            continue;
         double thG = std::atan2(psi * R, dz) * R2D;
         if (thG < 0)
            thG += 180.0;
         // REJECT HEMISPHERE-REFLECTED TRACKS. dp_specmat_C14.C has this cut and this macro did
         // not, which is the entire reason the 7 T Ex looked catastrophic (sigma 1.1-2.8 MeV)
         // while d_max and theta closed on truth at +-0.4 % and +-0.8 % on the same tracks. When
         // the z advance to the apex is small its SIGN is noisy, so theta comes out in the
         // forward hemisphere; cos(theta) enters the two-body inversion directly, so a reflected
         // angle does not merely add scatter, it puts Ex somewhere else entirely. The pitch
         // collapses with field, so reflections multiply at 7 T -- the same mechanism that
         // scrambles GENFIT's z-ordering near 90 deg.
         if (thG < thMin)
            continue;
         double sth = std::sin(thG / R2D);
         if (sth < 1e-3)
            continue;
         double pTot = (0.299792458 * bField * R) / sth;
         double keG = std::sqrt(pTot * pTot + mp * mp) - mp;
         if (!(keG > 0) || keG > 80)
            continue;

         auto [ex, tcm] = ffx_kine(M1, M2, M3, M4, Ebeam, thG / R2D, keG);
         if (!std::isfinite(ex) || std::fabs(ex) > 12)
            continue;
         exG[L].push_back(ex);
         exT[L].push_back(levEx[L]);
         cmG[L].push_back(tcm);
      }
      fr->Close();
      fs->Close();
   }

   auto *c = new TCanvas("cffx", "fitfree ex", 1500, 620);
   c->Divide(2, 1);

   // A: the spectrum. Frame on the data, fine bins -- never a sqrt(N) rule capped at some value.
   c->cd(1);
   gPad->SetGrid();
   double lo = 1e9, hi = -1e9;
   for (int L = 0; L < 3; ++L)
      for (double v : exG[L]) { lo = std::min(lo, v); hi = std::max(hi, v); }
   if (!(hi > lo)) { lo = -2; hi = 5; }
   lo = std::max(lo, -3.0);
   hi = std::min(hi, 8.0);
   auto *hs = new THStack("hs", TString::Format("A  E_{x} from hit geometry only, %s;E_{x} [MeV];tracks", tag.Data()));
   auto *lg = new TLegend(0.58, 0.66, 0.98, 0.90);
   lg->SetBorderSize(0);
   lg->SetFillStyle(0);
   printf("\n  === fit-free E_x, %s (B = %.2f T), backward theta_lab > %.0f deg ===\n", tag.Data(), bField, thMin);
   printf("  %-8s %6s  %9s  %9s  %9s\n", "level", "n", "true Ex", "median", "sigma");
   for (int L = 0; L < 3; ++L) {
      if (exG[L].size() < 5)
         continue;
      auto *h = new TH1D(TString::Format("hx%d", L), "", 180, lo, hi);
      for (double v : exG[L]) h->Fill(v);
      h->SetLineColor(col[L]);
      h->SetLineWidth(2);
      hs->Add(h);
      lg->AddEntry(h, TString::Format("E_{x} = %.3f  (n=%zu)", levEx[L], exG[L].size()), "l");
      double med = ffx_q(exG[L], .5), sig = (ffx_q(exG[L], .75) - ffx_q(exG[L], .25)) / 1.349;
      printf("  %-8s %6zu  %+9.3f  %+9.3f  %9.3f\n", levs[L], exG[L].size(), levEx[L], med, sig);
   }
   hs->Draw("nostack hist");
   lg->Draw();

   // B: Ex against theta_cm. A discrete state must be FLAT here; a slope is a kinematic
   // systematic that smears the integrated spectrum without ever looking like a bad fit.
   c->cd(2);
   gPad->SetGrid();
   double c0 = 1e9, c1 = -1e9;
   for (int L = 0; L < 3; ++L)
      for (double v : cmG[L]) { c0 = std::min(c0, v); c1 = std::max(c1, v); }
   if (!(c1 > c0)) { c0 = 0; c1 = 180; }
   auto *fB = gPad->DrawFrame(c0 - 2, lo, c1 + 2, hi,
                              "B  E_{x} vs #theta_{cm} -- a discrete state must be FLAT"
                              ";#theta_{cm} [deg];E_{x} [MeV]");
   fB->GetXaxis()->SetTitleSize(0.042);
   fB->GetYaxis()->SetTitleSize(0.042);
   for (int L = 0; L < 3; ++L) {
      if (exG[L].empty())
         continue;
      auto *g = new TGraph(exG[L].size(), cmG[L].data(), exG[L].data());
      g->SetMarkerStyle(20);
      g->SetMarkerSize(0.3);
      g->SetMarkerColor(col[L]);
      g->Draw("P same");
   }

   TString out = outDir + "dp_fitfree_ex_" + tag + ".png";
   c->SaveAs(out);
   printf("  GENFIT benchmark on this channel, backward: 0.197 MeV at 2.85 T, 0.241 at 4 T;\n"
          "  at 7 T GENFIT keeps only 625 of ~3000 tracks so its sigma is not quotable.\n");
   printf("  wrote %s\n\n", out.Data());
}
