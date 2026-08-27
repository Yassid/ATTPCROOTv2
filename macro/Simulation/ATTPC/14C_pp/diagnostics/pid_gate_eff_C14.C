/// @file pid_gate_eff_C14.C
/// @brief Measure the efficiency of the proton PID gate + backward cut ON THE SIMULATION.
///
/// WHY THIS IS THE MISSING NUMBER. Every 14C cross section is normalised on a luminosity fitted
/// to the elastic against a CALCULATED elastic cross section. The FRIBDAQ scalers give a second
/// luminosity that owes nothing to an optical model, and the two disagree:
///
///     eps = L_elastic / L_scaler = 0.59 (KD03), 0.80 (Menet), 0.86 (Becchetti-Greenlees)
///
/// eps is the detection efficiency that is NOT in the simulated acceptance. And it is not zero by
/// construction, because acceptance_C14.C applies exactly ONE selection -- chi2/ndf < chi2Cut --
/// while the data additionally pass, in pipeline/gate_events_C14.C:
///     * the hand-drawn proton polygon proton_14C.json, on (sqrt(dE/dx), Brho)
///     * theta_lab > 90 deg on the track
///     * the hardware trigger                              <-- NOT measurable here
/// The first two ARE measurable here: run the identical AtSpyralPID estimate and the identical
/// TCutG on simulated protons and count what survives. What comes out is an upper bound on eps
/// (the trigger can only lower it further).
///
/// The test that matters: if eps_measured comes out near 0.86, Becchetti-Greenlees is corroborated
/// by a route sharing nothing with the dip position or the shape, and the scaler luminosity
/// becomes a genuine cross-check of the absolute scale. If it comes out near 1.0, then the 14-41 %
/// has to be found somewhere else and the elastic/scaler disagreement is unexplained.
///
/// The denominator is deliberately the SAME population acceptance_C14.C calls reconstructed:
/// pattern-level tracks in events that produced one. This isolates the gate, which is the thing
/// the acceptance omits -- it does not re-measure the reconstruction.
///
/// MEASURED on gs_s1001 (8000 events), 2026-08-27:
///     3295 events with a pattern, 3183 pattern tracks, 3051 with a valid PID estimate
///     2826 inside the polygon (0.9263 of valid);  theta_lab > 90 keeps all of them
///     eps(gate, per track)  = 0.9263
///     eps(gate, per EVENT)  = 0.8577    <-- the number to compare
/// The per-event value is the right one: the acceptance is defined per generated reaction and the
/// luminosity counts beam particles. Against eps = L_elastic/L_scaler = 0.855 for
/// Becchetti-Greenlees, that agrees to 0.3 %, so the whole elastic/scaler deficit is the PID gate
/// and the trigger is consistent with unity. It also excludes Koning-Delaroche quantitatively:
/// its eps = 0.59 would need the trigger to cost a further factor 0.69 on top of the gate.
///
///   root -b -q 'pid_gate_eff_C14.C("/mnt/f/a1954_C14_acc_catima/gs_s1001_reco.root")'
/// *** THE FIELD SIGN IS +2.85 FOR SIMULATION, -2.85 FOR DATA. *** AtSpyralPID returns a SIGNED
/// Brho, so passing the data value here puts every simulated track at negative Brho while the
/// hand-drawn gate spans +0.040 to +0.964 -- and the efficiency comes back as a clean 0.0000,
/// which reads as a result rather than as a bug. The default below is the SIMULATION value
/// because that is what this macro is for.
void pid_gate_eff_C14(TString recoFile, Double_t bField = +2.85, Double_t thMin = 90.0,
                      TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/"
                                           "Unpack_HDF5/a1954/UKF/pid/proton_14C.json",
                      TString tag = "gs_s1001")
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   // THE GATE, read with the SAME loader gate_events_C14.C uses. Do not hand-roll this: the
   // json begins {"name":"proton_14C", "Z":1, "A":1, ...} and a naive "pull every number" parser
   // silently prepends 14, 1 and 1 to the vertex list, deforming the polygon. Anchor on
   // "vertices" and stop at the matching bracket.
   TCutG *pcut = nullptr;
   {
      std::ifstream in(protonGate.Data());
      if (!in) { printf("\033[1;31mno gate %s\033[0m\n", protonGate.Data()); return; }
      std::string str((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      auto pos = str.find('[', str.find("vertices"));
      if (pos == std::string::npos) { printf("\033[1;31mno vertices in %s\033[0m\n", protonGate.Data()); return; }
      std::vector<double> nums;
      const char *q = str.c_str() + pos, *end = str.c_str() + str.size();
      int depth = 0;
      while (q < end) {
         if (*q == '[') depth++;
         if (*q == ']') { depth--; if (depth <= 0) { ++q; break; } }
         char *np = nullptr; double v = strtod(q, &np);
         if (np != q) { nums.push_back(v); q = np; } else ++q;
      }
      pcut = new TCutG("pcut", nums.size() / 2);
      for (size_t i = 0, k = 0; i + 1 < nums.size(); i += 2, ++k) pcut->SetPoint(k, nums[i], nums[i + 1]);
      printf("\n  gate %s: %d vertices, x %.3f-%.3f, y %.3f-%.3f\n", protonGate.Data(), pcut->GetN(),
             *std::min_element(pcut->GetX(), pcut->GetX() + pcut->GetN()),
             *std::max_element(pcut->GetX(), pcut->GetX() + pcut->GetN()),
             *std::min_element(pcut->GetY(), pcut->GetY() + pcut->GetN()),
             *std::max_element(pcut->GetY(), pcut->GetY() + pcut->GetN()));
   }

   TFile *f = TFile::Open(recoFile);
   TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("cbmsim") : nullptr;
   if (!t) { printf("\033[1;31mcannot open %s\033[0m\n", recoFile.Data()); return; }
   TClonesArray *pe = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   auto *hAll = new TH2D("hAll", "simulated pattern tracks;#sqrt{dE/dx};B#rho [Tm]", 220, 0, 22, 200, 0, 1.3);
   auto *hIn = new TH2D("hIn", "", 220, 0, 22, 200, 0, 1.3);
   auto *hThAll = new TH1D("hThAll", "polar angle;#theta_{lab} [deg];tracks", 90, 0, 180);
   auto *hThIn = new TH1D("hThIn", "", 90, 0, 180);
   long long nTrk = 0, nValid = 0, nBack = 0, nGate = 0, nBoth = 0, nEvt = 0, nEvtKept = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!pe || pe->GetEntries() == 0) continue;
      auto *p = (AtPatternEvent *)pe->At(0);
      if (!p) continue;
      ++nEvt;
      bool kept = false;
      for (auto &trk : p->GetTrackCand()) {
         AtTrack &tr = const_cast<AtTrack &>(trk);
         ++nTrk;
         auto r = spy.Estimate(tr);
         if (!r.valid) continue;
         ++nValid;
         double th = r.polar * TMath::RadToDeg();
         hAll->Fill(r.sqrtdEdx, r.brho);
         hThAll->Fill(th);
         bool inG = pcut->IsInside(r.sqrtdEdx, r.brho), inB = th > thMin;
         if (inG) { ++nGate; hIn->Fill(r.sqrtdEdx, r.brho); }
         if (inB) ++nBack;
         if (inG && inB) { ++nBoth; hThIn->Fill(th); kept = true; }
      }
      if (kept) ++nEvtKept;
   }
   printf("\n  ===== proton gate efficiency on the simulation (%s) =====\n", tag.Data());
   printf("    events with a pattern             %10lld\n", nEvt);
   printf("    pattern tracks                    %10lld\n", nTrk);
   printf("    with a valid PID estimate         %10lld   (%.4f of tracks)\n", nValid, nTrk ? (double)nValid / nTrk : 0);
   printf("    inside the proton polygon         %10lld   (%.4f of valid)\n", nGate, nValid ? (double)nGate / nValid : 0);
   printf("    theta_lab > %.0f                    %10lld   (%.4f of valid)\n", thMin, nBack, nValid ? (double)nBack / nValid : 0);
   printf("    BOTH  (what the data require)     %10lld   (%.4f of valid)\n", nBoth, nValid ? (double)nBoth / nValid : 0);
   printf("    events keeping >=1 such track     %10lld   (%.4f of events)\n", nEvtKept, nEvt ? (double)nEvtKept / nEvt : 0);
   double eps = nValid ? (double)nBoth / nValid : 0;
   double epsEvt = nEvt ? (double)nEvtKept / nEvt : 0;
   printf("\n    eps(gate, per track)  = %.4f\n    eps(gate, per event)  = %.4f\n", eps, epsEvt);
   printf("\n    for comparison, eps = L_elastic/L_scaler:  KD03 0.59, Menet 0.80, Becchetti-Greenlees 0.86\n");
   printf("    this is an UPPER BOUND on eps -- the trigger is not simulated and can only lower it.\n");

   auto *c = new TCanvas("cpg", "", 1300, 560); c->Divide(2, 1);
   c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.13);
   hAll->Draw("colz"); pcut->SetLineColor(kRed + 1); pcut->SetLineWidth(3); pcut->Draw("l same");
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);
   tx.DrawLatex(0.15, 0.86, Form("simulated protons, gate keeps %.3f", eps));
   c->cd(2); gPad->SetGridx(); gPad->SetGridy();
   hThAll->SetLineColor(kGray + 2); hThAll->SetLineWidth(2); hThAll->Draw("hist");
   hThIn->SetLineColor(kAzure + 2); hThIn->SetLineWidth(2); hThIn->SetFillColorAlpha(kAzure + 2, 0.35);
   hThIn->Draw("hist same");
   auto *lg = new TLegend(0.15, 0.72, 0.58, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hThAll, "all valid PID estimates", "l");
   lg->AddEntry(hThIn, "gate + backward", "f"); lg->Draw();
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c->SaveAs(out + "12_pid_gate_eff_" + tag + ".png");
   printf("\n  wrote %s12_pid_gate_eff_%s.png\n\n", out.Data(), tag.Data());
}
