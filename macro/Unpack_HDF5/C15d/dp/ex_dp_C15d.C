/// @file ex_dp_C15d.C
/// @brief 15C(d,p)16C excitation energy from the fitted protons, and an Ebeam scan.
///
///   root -b -q 'dp/ex_dp_C15d.C()'                 // Ebeam from the (d,d') ridge
///   root -b -q 'dp/ex_dp_C15d.C(95)'               // try another beam energy
///
/// Two-body inversion, nothing fitted beyond the proton itself:
///   E3 = KE + m3,  p3 = sqrt(E3^2 - m3^2)
///   E4 = E1 + m2 - E3,  p4 = p1 - p3  (vector)
///   Ex = sqrt(E4^2 - p4^2) - m4
///
/// ★ Ebeam IS NOT INDEPENDENTLY KNOWN HERE. 90 MeV comes from this analysis's own 15C(d,d) elastic
/// ridge. The ground state sitting at Ex = 0 is therefore a CONSISTENCY CHECK between two channels
/// rather than a calibration -- if the g.s. lands away from zero, either the beam energy or the
/// energy scale is wrong, and the scan below says which beam energy would put it at zero.
///
/// ★ FORWARD AND BACKWARD ARE KEPT SEPARATE. In inverse kinematics the two hemispheres have very
/// different resolution -- backward protons are slow and lose a large fraction of their energy in
/// the gas -- so a combined spectrum hides which half carries the structure. They are also the two
/// halves that a seeding or ordering error would move relative to each other, which is exactly the
/// failure the upstream "fit backward tracks from the vertex end" commit fixed.

namespace {
const double kU = 931.49410242;
const double kM1 = 15.0105993 * kU;  // 15C beam
const double kM2 = 2.0141018 * kU;   // d target
const double kM3 = 1.0078250 * kU;   // p ejectile
const double kM4 = 16.0147013 * kU;  // 16C residual

/// Ex of the residual from the ejectile's lab kinetic energy and angle.
double ExFrom(double keLab, double thDeg, double Ebeam)
{
   const double E1 = Ebeam + kM1;
   const double p1 = std::sqrt(E1 * E1 - kM1 * kM1);
   const double E3 = keLab + kM3;
   const double p3 = std::sqrt(std::max(0.0, E3 * E3 - kM3 * kM3));
   const double th = thDeg * TMath::DegToRad();
   const double E4 = E1 + kM2 - E3;
   const double p4z = p1 - p3 * std::cos(th);
   const double p4t = p3 * std::sin(th);
   const double m4sq = E4 * E4 - (p4z * p4z + p4t * p4t);
   if (m4sq <= 0)
      return -999;
   return std::sqrt(m4sq) - kM4;
}
} // namespace

void ex_dp_C15d(Double_t Ebeam = 90.0, TString fitDir = "/home/yassid/C15d_fit/",
                TString outDir = "dp/plots/", Double_t chi2Cut = 5.0, Double_t exLo = -6, Double_t exHi = 14)
{
   gSystem->mkdir(outDir, kTRUE);
   TChain ch("kin");
   const int nf = ch.Add(fitDir + "*_kin_p.root");
   Int_t run, event, track, ndf, fwd;
   Double_t ke, th, keX, thX, vz, c2;
   ch.SetBranchAddress("run", &run);
   ch.SetBranchAddress("keXtr", &keX);
   ch.SetBranchAddress("thetaXtr", &thX);
   ch.SetBranchAddress("vz", &vz);
   ch.SetBranchAddress("chi2ndf", &c2);
   ch.SetBranchAddress("dirFwd", &fwd);

   auto *hAll = new TH1D("exAll", Form("15C(d,p)16C at E_{beam}=%.0f MeV;E_{x}(^{16}C) [MeV];counts", Ebeam),
                         200, exLo, exHi);
   auto *hF = new TH1D("exF", "forward;E_{x} [MeV];counts", 200, exLo, exHi);
   auto *hB = new TH1D("exB", "backward;E_{x} [MeV];counts", 200, exLo, exHi);
   auto *hKT = new TH2D("hKT", "protons;#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 200, 0, 40);

   std::vector<float> K, T;
   Long64_t n = 0;
   for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
      ch.GetEntry(i);
      if (!(keX > 0) || c2 > chi2Cut)
         continue;
      ++n;
      K.push_back(keX);
      T.push_back(thX);
      hKT->Fill(thX, keX);
      const double ex = ExFrom(keX, thX, Ebeam);
      if (ex < -900)
         continue;
      hAll->Fill(ex);
      (thX > 90 ? hB : hF)->Fill(ex);
   }

   std::cout << "\033[1;33m=== 15C(d,p)16C excitation energy ===\033[0m\n"
             << "  kin files : " << nf << "\n"
             << "  protons   : " << n << "  (chi2/ndf < " << chi2Cut << ")\n"
             << "  forward   : " << hF->GetEntries() << "   backward : " << hB->GetEntries() << "\n"
             << "  Ebeam     : " << Ebeam << " MeV  (from the (d,d) elastic ridge -- NOT independent)\n";

   auto peak = [](TH1D *h, double lo, double hi) {
      int b1 = h->FindBin(lo), b2 = h->FindBin(hi), ib = -1;
      double best = -1;
      for (int b = b1; b <= b2; ++b)
         if (h->GetBinContent(b) > best) { best = h->GetBinContent(b); ib = b; }
      return ib > 0 ? h->GetBinCenter(ib) : -999;
   };
   printf("  peak nearest 0 : all %.2f   forward %.2f   backward %.2f MeV\n",
          peak(hAll, -3, 3), peak(hF, -3, 3), peak(hB, -3, 3));

   // ---- what beam energy puts the ground state at zero? --------------------------------------
   printf("\n  %-10s %10s %10s %10s\n", "Ebeam", "peak(all)", "peak(fwd)", "peak(bwd)");
   for (double E = Ebeam - 20; E <= Ebeam + 20; E += 5) {
      TH1D a("a", "", 200, exLo, exHi), fh("fh", "", 200, exLo, exHi), bh("bh", "", 200, exLo, exHi);
      for (size_t i = 0; i < K.size(); ++i) {
         const double ex = ExFrom(K[i], T[i], E);
         if (ex < -900) continue;
         a.Fill(ex);
         (T[i] > 90 ? bh : fh).Fill(ex);
      }
      printf("  %-10.0f %10.2f %10.2f %10.2f\n", E, peak(&a, -3, 3), peak(&fh, -3, 3), peak(&bh, -3, 3));
   }
   printf("\n  16C bound/low-lying states for reference: 1.766 (2+), 3.027 (2+), 3.986, 4.088\n");

   TFile fo(outDir + "ex_dp_C15d.root", "RECREATE");
   hAll->Write(); hF->Write(); hB->Write(); hKT->Write();
   auto *c = new TCanvas("cex", "ex", 1200, 900);
   c->Divide(2, 2);
   c->cd(1); hKT->Draw("colz"); gPad->SetLogz();
   c->cd(2); hAll->Draw("hist");
   c->cd(3); hF->SetLineColor(kAzure + 1); hF->Draw("hist");
   c->cd(4); hB->SetLineColor(kRed + 1); hB->Draw("hist");
   c->SaveAs(outDir + "ex_dp_C15d.png");
   fo.Close();
   std::cout << "  \033[1;32mwrote\033[0m " << outDir << "ex_dp_C15d.{root,png}\n";
}
