/// @file acceptance_C16dt.C
/// @brief Detector acceptance vs theta_cm for 16C(d,t)15C, from MC truth, per state.
///
///   acceptance(theta_cm) = (generated reactions that end as a gated, Spyral-valid triton)
///                        / (reactions generated)
///
/// DENOMINATOR IS TRUTH. For every generated event the triton is taken from the MCTrack list by
/// SPECIES (pdg 1000010030), never by index, and its true (KE, theta_lab) is pushed through the
/// same two-body expressions the analysis uses to get the true theta_cm. So no reconstruction
/// inefficiency can touch it -- which is the entire point of a denominator.
///
/// NUMERATOR IS BINNED AT THE TRUE ANGLE, NOT THE RECONSTRUCTED ONE. This is the subtle one.
/// Binning the numerator by reconstructed theta_cm would fold the angular resolution into the
/// acceptance; the yields this acceptance corrects are ALREADY binned by reconstructed angle, so
/// the migration would then be applied twice. Truth on both sides, and the ratio is a pure
/// efficiency.
///
/// NO FIT IS REQUIRED, and that is deliberate. The numerator is a pattern track that (a) yields a
/// valid AtSpyralPID estimate and (b) falls inside the PID gate -- exactly the chain the data
/// passes through before any fitting. Requiring a converged genfit fit as well would fold the
/// fitter's convergence into a number that is supposed to describe the detector, and the fitter
/// is a moving target (its collapse rate went 2.51% -> 1.08% in one day).
///
/// PER STATE, because the kinematics differ: at Ex = 0 the tritons already sit at the 56 deg lab
/// limit, so the higher states are not a small correction to the ground state's acceptance.
///
/// The theta_lab window is applied to TRUTH on both numerator and denominator when it is enabled,
/// so it describes the data's cut rather than silently deleting the denominator's tails.
///
///   root -b -q 'acceptance_C16dt.C()'
#include <tuple>

static double dt_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics -> {Ex, theta_cm [deg]}, the same form the (p,d) acceptance uses
static std::tuple<double, double> dt_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * dt_omega2(s, m1 * m1, m2 * m2) * dt_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (dt_omega2(s, m1 * m1, m2 * m2) * dt_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

void acceptance_C16dt(TString simDir = "/mnt/f/a1975_C16_dt_sim/",
                      TString statesCSV = "gs_s3001,ex1_s3001,ex2_s3001,ex3_s3001,ex4_s3001",
                      TString exCSV = "0.0,0.740,3.103,4.657,6.358",
                      TString gateFile = "triton_dt_sim.json", Double_t Ebeam = 184.25, Int_t nBins = 28,
                      Double_t cmMax = 70.0, Double_t thLabMin = 0.0, Double_t thLabMax = 180.0,
                      // Vertex window from the explorer selection of 2026-08-18: the range over
                      // which the per-state z distributions are homogeneous.
                      Double_t vzLo = 10.0, Double_t vzHi = 500.0,
                      // THE SIMULATION'S DRIFT z IS MIRRORED against the experiment: measured on
                      // gs_s3001, corr(z_true, z_reco) = -1.000 with z_true + z_reco = 981.8 mm.
                      // Applying the data's vz window to the raw reconstructed z would therefore
                      // select the OPPOSITE half of the target and silently invert the correction.
                      // The numerator un-mirrors before cutting; 0 disables (real data).
                      Double_t zMirror = 981.8,
                      TString outFile = "data/acceptance_dt.root", TString outPng = "data/acceptance_dt.png")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);

   const double u = 931.49401;
   // 16C(d,t)15C: DETECTED particle is the triton, residual is 15C.
   const double m_C16 = 16.0147013 * u, m_d = 2.0135532 * u;
   const double m_t = 3.01550072 * u, m_C15 = 15.0105993 * u;

   // ---- the PID gate, read straight from its JSON ------------------------------------------
   std::ifstream in(gateFile.Data());
   if (!in) { printf("cannot open gate %s\n", gateFile.Data()); return; }
   std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   std::vector<double> vx, vy;
   size_t vp = all.find("\"vertices\"");
   size_t p = (vp == std::string::npos) ? std::string::npos : all.find('[', vp + 10);
   while (p != std::string::npos) {
      size_t a = all.find('[', p + 1);
      if (a == std::string::npos) break;
      size_t b = all.find(']', a);
      if (b == std::string::npos) break;
      double gx, gy;
      if (sscanf(all.substr(a + 1, b - a - 1).c_str(), "%lf , %lf", &gx, &gy) == 2) { vx.push_back(gx); vy.push_back(gy); }
      p = b;
      size_t nxt = all.find_first_not_of(" \n\r\t,", b + 1);
      if (nxt == std::string::npos || all[nxt] == ']') break;
   }
   if (vx.size() < 3) { printf("gate has only %zu vertices\n", vx.size()); return; }
   TCutG gate("gate", vx.size());
   for (size_t i = 0; i < vx.size(); ++i) gate.SetPoint(i, vx[i], vy[i]);
   printf("gate %s: %zu vertices\n", gateFile.Data(), vx.size());
   printf("theta_lab window %.0f-%.0f deg applied to TRUTH on BOTH numerator and denominator\n\n",
          thLabMin, thLabMax);

   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   TObjArray *ta = statesCSV.Tokenize(",");
   TObjArray *te = exCSV.Tokenize(",");
   TFile *fo = TFile::Open(outFile, "RECREATE");
   auto *c = new TCanvas("cacc", "acceptance", 1100, 750);
   TLegend *leg = new TLegend(0.14, 0.14, 0.44, 0.36);
   leg->SetBorderSize(0);
   const int col[5] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1};

   for (int is = 0; is < ta->GetEntries(); ++is) {
      TString tg = ((TObjString *)ta->At(is))->GetString().Strip(TString::kBoth);
      double resEx = (is < te->GetEntries()) ? ((TObjString *)te->At(is))->GetString().Atof() : 0.0;

      TString fs = simDir + tg + "_sim.root", fr = simDir + tg + "_reco.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) { printf("skip %s\n", tg.Data()); continue; }
      TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
      TTree *ts = Fs ? (TTree *)Fs->Get("cbmsim") : nullptr;
      TTree *tr = Fr ? (TTree *)Fr->Get("cbmsim") : nullptr;
      if (!ts || !tr) { printf("skip %s (no tree)\n", tg.Data()); continue; }

      // ONE ENTRY PER GENERATED EVENT ON BOTH SIDES, or the event-by-event match is a guess.
      if (ts->GetEntries() != tr->GetEntries()) {
         printf("  *** %s: sim has %lld entries, reco %lld -- REFUSING to match by index ***\n",
                tg.Data(), ts->GetEntries(), tr->GetEntries());
         Fs->Close(); Fr->Close();
         continue;
      }

      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pe);

      // SetDirectory(nullptr) is NOT tidiness. TFile::Open above made the SIM file the current
      // directory, so a histogram booked here would be OWNED BY IT -- and Fs->Close() at the end
      // of this iteration would delete it, leaving the canvas holding freed memory and crashing
      // at SaveAs after the whole analysis had already run. Memory-resident, written explicitly.
      auto *hGen = new TH1D("hGen_" + tg, "", nBins, 0, cmMax);
      auto *hRec = new TH1D("hRec_" + tg, "", nBins, 0, cmMax);
      hGen->SetDirectory(nullptr); hRec->SetDirectory(nullptr);
      hGen->Sumw2(); hRec->Sumw2();

      long nGen = 0, nRec = 0;
      for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
         ts->GetEntry(i);
         // --- truth: find the triton by SPECIES, never by index ---
         double cmT = -1, thLabT = -1, zTrue = -1e9;
         for (int k = 0; k < mc->GetEntries(); ++k) {
            auto *t = (AtMCTrack *)mc->At(k);
            if (!t || t->GetPdgCode() != 1000010030) continue; // triton
            zTrue = t->GetStartZ() * 10.0; // cm -> mm
            double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000; // GeV -> MeV
            double pm = std::sqrt(px * px + py * py + pz * pz);
            if (pm <= 0) continue;
            double KE = std::sqrt(pm * pm + m_t * m_t) - m_t;
            double th = std::acos(pz / pm);
            auto [ex, cm] = dt_kine(m_C16, m_d, m_t, m_C15, Ebeam, th, KE);
            if (!std::isfinite(cm)) continue;
            cmT = cm;
            thLabT = th * TMath::RadToDeg();
            break;
         }
         if (cmT < 0) continue;
         // the data's theta_lab window, on TRUTH, so numerator and denominator see one selection
         if (thLabMax > thLabMin && (thLabT < thLabMin || thLabT > thLabMax)) continue;
         // DENOMINATOR: reactions that really happened inside the vertex window, on TRUTH.
         if (vzHi > vzLo && (zTrue < vzLo || zTrue > vzHi)) continue;
         ++nGen;
         hGen->Fill(cmT);

         // --- numerator: a gated, Spyral-valid triton in THIS event, binned at the TRUE angle ---
         tr->GetEntry(i);
         if (!pe || pe->GetEntries() == 0) continue;
         auto *P = (AtPatternEvent *)pe->At(0);
         if (!P) continue;
         bool found = false;
         for (auto &trk : P->GetTrackCand()) {
            AtTrack &T = const_cast<AtTrack &>(trk);
            auto r = spy.Estimate(T);
            if (!r.valid) continue;
            if (!gate.IsInside(r.sqrtdEdx, r.brho)) continue;
            // NUMERATOR: the data's vertex cut, applied to the RECONSTRUCTED z after un-mirroring,
            // so the numerator carries the z resolution the data's cut also suffers.
            if (vzHi > vzLo) {
               double zr = (zMirror > 0) ? (zMirror - r.vertex.Z()) : r.vertex.Z();
               if (zr < vzLo || zr > vzHi) continue;
            }
            // truth-match the track itself, so a gated contaminant cannot count as a success
            int nT = 0, nTot = 0;
            for (const auto &hit : T.GetHitArray()) {
               if (!hit) continue;
               const auto &mcs = hit->GetMCSimPointArray();
               if (mcs.empty()) continue;
               ++nTot;
               if (mcs[0].A == 3 && mcs[0].Z == 1) ++nT;
            }
            if (nTot > 0 && nT > nTot / 2) { found = true; break; }
         }
         if (found) { ++nRec; hRec->Fill(cmT); }
      }

      auto *hAcc = (TH1D *)hRec->Clone("hAcc_" + tg);
      hAcc->SetDirectory(nullptr);
      hAcc->Divide(hRec, hGen, 1, 1, "B"); // binomial errors
      hAcc->SetTitle(TString::Format("16C(d,t)15C acceptance;#theta_{cm} [deg];acceptance"));
      hAcc->SetMinimum(0); hAcc->SetMaximum(1.05);
      hAcc->SetLineColor(col[is % 5]); hAcc->SetMarkerColor(col[is % 5]);
      hAcc->SetMarkerStyle(20); hAcc->SetMarkerSize(0.8);

      printf("=== %-10s Ex = %.3f MeV   generated %ld, accepted %ld  -> overall %.3f\n",
             tg.Data(), resEx, nGen, nRec, nGen ? (double)nRec / nGen : 0.0);
      for (int b = 1; b <= hGen->GetNbinsX(); ++b) {
         double g = hGen->GetBinContent(b);
         if (g < 1) continue;
         printf("   %5.1f-%5.1f  gen %6.0f  acc %6.0f   %.3f +- %.3f\n", hGen->GetBinLowEdge(b),
                hGen->GetBinLowEdge(b) + hGen->GetBinWidth(b), g, hRec->GetBinContent(b),
                hAcc->GetBinContent(b), hAcc->GetBinError(b));
      }
      printf("\n");

      fo->cd();
      hGen->Write(); hRec->Write(); hAcc->Write();
      c->cd();
      hAcc->Draw(is == 0 ? "E1" : "E1 SAME");
      leg->AddEntry(hAcc, TString::Format("Ex = %.3f MeV", resEx), "lp");
      Fs->Close(); Fr->Close();
   }
   c->cd(); leg->Draw();
   c->SaveAs(outPng);
   fo->Close();
   printf("acceptance -> %s\nplot       -> %s\n", outFile.Data(), outPng.Data());
}
