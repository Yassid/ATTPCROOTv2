/// @file dedxeom_calib_dp.C
/// @brief Per-run charge -> MeV/mm constant k for the species-blind EOM (pid_dedxeom_dp.C), from the deuteron band
///        (2026-11-01). Tracks with the hit-refit/arc-angle Brho bA in [bLo, bHi] whose ln(hs * bA) lies within
///        +-win of the band mode (deuterons dominate) give k_i = meanS(p = bA, m_d, lseg) / hs^2 (CATIMA, velocity-
///        scaled proton table); k = median. species (default 1 = d) selects the band mass, e.g. 2 for a triton-only
///        simulation (closure test). eomPrefix: momentum from the deuteron-hypothesis EOM fit (pid_eom_dp.C trees,
///        <eomPrefix><run>.root, same track order: p = e0b1 * sin(hth) / sin(ta)) instead of the circle bA -- the circle
///        momentum is biased low by deceleration and inflated k by ~16 % (Bragg-energy calibration, 2026-11-01). Also checks the velocity scaling against the d and t tables.
///   root -l -b -q 'dedxeom_calib_dp.C("/mnt/f/a1975/pid3d_dedxeom/A_","0016,0026,0058,0088")'
#include "eom_brho_dp.h"
void dedxeom_calib_dp(TString prefix, TString runsCSV, double bLo = 0.2, double bHi = 0.4, double win = 0.08, int species = 1,
                      TString eomPrefix = "")
{
   gSystem->Load("libAtTools.so");
   eom::Table tab[3];
   eom::makeTables(tab);
   printf("velocity scaling check (MeV/mm): KE 2 / 5 / 10 MeV\n");
   for (int k = 1; k < 3; ++k) {
      printf("  %s table:", k == 1 ? "d" : "t");
      for (double ke : {2.0, 5.0, 10.0}) printf("  %.5f (scaled p %.5f)", tab[k].dedx(ke), eom::sVel(tab[0], ke, eom::massMeV(k)));
      printf("\n");
   }
   for (auto *o : *runsCSV.Tokenize(",")) {
      TString r = ((TObjString *)o)->GetString();
      TFile f(prefix + r + ".root");
      auto *w = (TTree *)f.Get("w");
      TString sel = Form("va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&lseg>5&&bA>%g&&bA<%g", bLo, bHi);
      auto *h = new TH1D(Form("h%s", r.Data()), "", 800, 0, 12);
      w->Draw(Form("log(hs*bA*100)>>h%s", r.Data()), sel, "goff");
      h->Smooth(2);
      double mode = h->GetBinCenter(h->GetMaximumBin());
      h->SetDirectory(nullptr);
      double hs, bA, lseg;
      int va, v3, ok1, nh;
      double sp3;
      w->SetBranchAddress("hs", &hs); w->SetBranchAddress("bA", &bA); w->SetBranchAddress("lseg", &lseg);
      w->SetBranchAddress("va", &va); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("ok1", &ok1);
      w->SetBranchAddress("nh", &nh); w->SetBranchAddress("sp3", &sp3);
      std::vector<double> kk;
      TFile *fe = eomPrefix.Length() ? TFile::Open(eomPrefix + r + ".root") : nullptr;
      TTree *we = fe ? (TTree *)fe->Get("w") : nullptr;
      double e0b1 = -1, hthE = -1, ta = -1, brA = 0, brE = 0;
      int nhE = 0;
      long bad = 0;
      if (we) {
         we->SetBranchAddress("e0b1", &e0b1); we->SetBranchAddress("hth", &hthE); we->SetBranchAddress("nh", &nhE);
         we->SetBranchAddress("br", &brE);
         w->SetBranchAddress("ta", &ta); w->SetBranchAddress("br", &brA);
      }
      for (Long64_t i = 0; i < w->GetEntries(); ++i) {
         w->GetEntry(i);
         if (we) { we->GetEntry(i); if (nhE != nh || std::abs(brE - brA) > 1e-9) { ++bad; continue; } }
         if (!va || !v3 || !ok1 || sp3 >= 6 || nh < 30 || hs <= 0 || lseg <= 5 || bA < bLo || bA > bHi) continue;
         if (std::abs(std::log(hs * bA * 100) - mode) > win) continue;
         double p = bA * 299.792458;
         if (we) {
            if (e0b1 <= 0 || hthE <= 0 || ta <= 0) continue;
            p = e0b1 * 299.792458 * std::sin(hthE) / std::sin(ta);
         }
         kk.push_back(eom::meanS(tab[0], p, eom::massMeV(species), lseg) / (hs * hs));
      }
      std::sort(kk.begin(), kk.end());
      if (kk.size() < 20) { printf("run %s: only %zu band tracks (mode %.3f)\n", r.Data(), kk.size(), mode); continue; }
      if (we) printf("  (EOM momentum; alignment mismatches %ld)\n", bad);
      printf("run %s: %zu band tracks, k = %.6g MeV/mm per (charge/mm)  (p16 %.4g p84 %.4g)\n", r.Data(), kk.size(),
             kk[kk.size() / 2], kk[kk.size() * 16 / 100], kk[kk.size() * 84 / 100]);
   }
}
