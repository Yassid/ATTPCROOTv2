/// @file order_score_p_dp.C
/// @brief Backward (d,p) protons on data (2026-11-02): proton-band tracks of the circle PID, runs 0016 0026 0058 0088, fitted with the
///        proton hypothesis by order_fit_dp.C: mode 0 (no stray-hit removal) and mode 9 (removal with the large-pad guard); plus the adopted
///        gf_dp_find805 production read back for the same (event, trackID). Adopted selection: theta 89-161, KE > 1.2 MeV, chi2/ndf < 13.5,
///        vz -100..940 mm; Ex from 16C(d,p)17C (Q = -1.496) with Ebeam(vz) = 184.25 - 0.01103 (vz - 500) and keOff +0.25 MeV (the kinematic
///        correction of the explorer is not applied). Reports: selected, bound window -1.5..1.0, junk Ex < -2 (no physical state), off-axis
///        vertices > 20 mm, peak width; writes the spectra.
///   root -l -b -q 'order_score_p_dp.C("plots/order_dp_ex.png")'
namespace ordp {
double exDP(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613, mp = 938.272, M4 = mC + md - mp + 1.496;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ep = ke + mp, pp = std::sqrt(Ep * Ep - mp * mp);
   const double E4 = Eb + md - Ep, px = -pp * std::sin(th), pz = pb - pp * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - M4;
}
struct Tally { long sel{0}, bound{0}, junk{0}, off{0}; };
void add(Tally &t, TH1D *h, double ke, double th, double chi, double vz, double vr)
{
   const double deg = th * 57.29578;
   if (deg < 89 || deg > 161 || ke <= 1.2 || chi >= 13.5 || vz < -100 || vz > 940) return;
   const double ex = exDP(ke + 0.25, th, 184.25 - 0.01103 * (vz - 500));
   ++t.sel;
   h->Fill(ex);
   if (ex > -1.5 && ex < 1.0) ++t.bound;
   if (ex < -2) ++t.junk;
   if (vr > 20) ++t.off;
}
}
void order_score_p_dp(TString png, TString O = "/mnt/f/a1975/order_test/", TString F = "/mnt/f/a1975/gf_dp_find805/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   gStyle->SetOptStat(0);
   const char *runs[4] = {"0016", "0026", "0058", "0088"};
   TH1D *h[3];
   const char *nm[3] = {"production clusters", "stray hits dropped, large-pad guard", "adopted find805"};
   ordp::Tally T[3];
   long nAll = 0, nFit[3] = {0}, nGuardCh = 0;
   for (int k = 0; k < 3; ++k) h[k] = new TH1D(Form("h%d", k), ";E_{x} [MeV];tracks / 0.2 MeV", 70, -5, 9);
   for (int ir = 0; ir < 4; ++ir) {
      TChain o("o");
      o.Add(O + "p_" + runs[ir] + "_0.root"); o.Add(O + "p_" + runs[ir] + "_1.root");
      Long64_t ev; int tid, nDrop, ok[2];
      double ke[2], th[2], c2[2], nd[2], vr[2], vz[2];
      o.SetBranchAddress("ev", &ev); o.SetBranchAddress("tid", &tid); o.SetBranchAddress("nDrop", &nDrop);
      const int mode[2] = {0, 9};
      for (int m = 0; m < 2; ++m) {
         o.SetBranchAddress(Form("ok%d", mode[m]), &ok[m]); o.SetBranchAddress(Form("ke%d", mode[m]), &ke[m]); o.SetBranchAddress(Form("th%d", mode[m]), &th[m]);
         o.SetBranchAddress(Form("c2%d", mode[m]), &c2[m]); o.SetBranchAddress(Form("nd%d", mode[m]), &nd[m]); o.SetBranchAddress(Form("vr%d", mode[m]), &vr[m]); o.SetBranchAddress(Form("vz%d", mode[m]), &vz[m]);
      }
      TFile ff(F + "run_" + runs[ir] + "_multifit_genfitter_p.root");
      auto *tf = (TTree *)ff.Get("cbmsim");
      TClonesArray *te = nullptr, *pe = nullptr;
      tf->SetBranchAddress("AtTrackingEvent", &te); tf->SetBranchAddress("AtPIDEvent", &pe);
      for (Long64_t i = 0; i < o.GetEntries(); ++i) {
         o.GetEntry(i);
         ++nAll;
         for (int m = 0; m < 2; ++m) {
            if (!ok[m]) continue;
            ++nFit[m];
            ordp::add(T[m], h[m], ke[m], th[m], c2[m] / nd[m], vz[m], vr[m]);
         }
         if (ev >= tf->GetEntries()) continue;
         tf->GetEntry(ev);
         if (!pe->GetEntries() || !te->GetEntries()) continue;
         auto *pid = (AtPIDEvent *)pe->At(0);
         if (!pid || pid->GetSpyral().empty()) continue;
         auto *tev = (AtTrackingEvent *)te->At(0);
         for (auto &ft : tev->GetFittedTracks()) {
            if (!ft || ft->GetTrackID() != tid) continue;
            auto &kx = ft->GetKinematicsXtr();
            const double ndf = ft->GetTrackMetadata()->GetNdf();
            if (kx.kineticEnergy <= 0 || ndf <= 0) break;
            ++nFit[2];
            auto v = ft->GetVertex();
            ordp::add(T[2], h[2], kx.kineticEnergy, kx.theta, ft->GetTrackMetadata()->GetChi2() / ndf, v.Z(), std::hypot(v.X(), v.Y()));
            break;
         }
      }
   }
   printf("proton-band tracks, 4 runs: %ld\n", nAll);
   const int col[3] = {kGray + 2, kBlue + 2, kRed + 1};
   auto *c = new TCanvas("c", "c", 1100, 700);
   gPad->SetGrid();
   double ymax = 0;
   for (int k = 0; k < 3; ++k) ymax = std::max(ymax, h[k]->GetMaximum());
   auto *lg = new TLegend(0.55, 0.72, 0.89, 0.89);
   for (int k = 0; k < 3; ++k) {
      TF1 g(Form("g%d", k), "gaus", -1.5, 1.0);
      h[k]->Fit(&g, "QNR");
      printf("  %-36s fit %5.1f %%  selected %5ld  bound(-1.5,1) %5ld  junk(Ex<-2) %5ld (%.1f %%)  off-axis>20 %5ld (%.1f %%)  peak %+.2f, sigma %.2f MeV\n", nm[k], 100. * nFit[k] / nAll, T[k].sel, T[k].bound, T[k].junk, 100. * T[k].junk / std::max(1L, T[k].sel), T[k].off, 100. * T[k].off / std::max(1L, T[k].sel), g.GetParameter(1), g.GetParameter(2));
      h[k]->SetLineColor(col[k]); h[k]->SetLineWidth(2); h[k]->SetMaximum(1.15 * ymax);
      h[k]->Draw(k ? "hist same" : "hist");
      lg->AddEntry(h[k], nm[k], "l");
   }
   lg->Draw();
   c->SaveAs(png);
}
