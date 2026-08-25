/// @file phasespace_1n_C14.C
/// @brief Generate the 1-neutron-emission phase space for 14C(p,p')14C* -> 13C + n, in the SAME
///        reconstructed-Ex variable the data uses.
///
/// Sn(14C) = 8.176 MeV, so everything above that in the (p,p') spectrum can be continuum from
/// neutron emission -- including the 8.317 state, which is unbound by 141 keV. A fit over 6-10 MeV
/// that has no continuum term has to absorb it into the peak areas or a straight line, and that is
/// what makes the upper half of the spectrum fit badly.
///
/// Ported from Downloads/PhaseSpace_16Cpp16C_test.C (the 16C(p,p') version, same structure).
/// The proton is generated in a 3-body final state p' + 13C + n and its Ex is then computed with
/// the TWO-BODY expression, exactly as the data analysis does -- that mismatch is the point: it is
/// what gives the continuum its shape in the reconstructed variable.
///
///   root -b -q 'phasespace_1n_C14.C(200000)'
void phasespace_1n_C14(Int_t nEvents = 200000, Double_t EbeamTot = 159.75, Double_t sigEx = 0.15,
                       Double_t cmMin = 20.0, Double_t cmMax = 140.0,
                       TString out = "plots/phasespace_1n_C14.root")
{
   if (!gROOT->GetClass("TGenPhaseSpace")) gSystem->Load("libPhysics");
   const double u = 931.49410242;
   const double m1 = 14.0032420 * u;    // 14C beam
   const double m2 = 1.00782503 * u;    // p target
   const double m3 = 1.00782503 * u;    // p' ejectile
   const double m4 = 14.0032420 * u;    // 14C recoil, for the two-body Ex
   const double m5 = 13.00335484 * u;   // 13C after 1n emission
   const double mn = 1.00866492 * u;    // neutron
   const double T1 = EbeamTot;          // beam KE, MeV (total, not per nucleon)
   const double p1 = TMath::Sqrt(T1 * (T1 + 2 * m1));

   TLorentzVector target(0, 0, 0, m2), beam(0, 0, p1, T1 + m1);
   TLorentzVector W = beam + target;
   Double_t masses[3] = {m5, m3, mn};   // 13C + p' + n
   TGenPhaseSpace ev;
   if (!ev.SetDecay(W, 3, masses)) { printf("\033[1;31mSetDecay failed -- not enough energy\033[0m\n"); return; }

   // two-body Ex from the proton, the SAME expression ex_C14.C uses
   auto ExOf = [&](double K, double thLab) {
      double Et1 = T1 + m1, Et3 = K + m3;
      double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
      double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
      auto om = [](double x, double y, double z) { return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); };
      double a = (std::cos(thLab) * om(s, m1*m1, m2*m2) * om(uu, m2*m2, m3*m3)
                  - (s - m1*m1 - m2*m2) * (m2*m2 + m3*m3 - uu)) / (2*m2*m2) + s + uu - m2*m2;
      return a < 0 ? -99.0 : std::sqrt(a) - m4;
   };
   // theta_cm of the proton, for the same angular window as the data
   double etot = T1 + m1 + m2, etotcm = TMath::Sqrt(m1*m1 + m2*m2 + 2*m2*(T1+m1));
   double gam = etot / etotcm, beta = TMath::Sqrt(1 - 1/gam/gam);

   auto *h = new TH1D("hPS", "1n phase space, reconstructed E_{x};E_{x} [MeV];weight", 300, 4, 16);
   TRandom3 rnd(20260826);
   double wsum = 0;
   for (Int_t i = 0; i < nEvents; ++i) {
      double w = ev.Generate();
      TLorentzVector *pp = ev.GetDecay(1);           // the proton
      double K = pp->E() - m3, th = pp->Theta();
      // proton theta_cm, boosting along the beam
      TLorentzVector q = *pp;
      q.Boost(0, 0, -beta);
      double thcm = q.Theta() * TMath::RadToDeg();
      if (thcm < cmMin || thcm > cmMax) continue;
      double ex = ExOf(K, th);
      if (ex < -90) continue;
      h->Fill(ex + rnd.Gaus(0, sigEx), w);           // smeared with the experimental resolution
      wsum += w;
   }
   h->Smooth(1);
   TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile fo(dir + "/" + out, "RECREATE");
   h->Write();
   fo.Close();
   printf("  generated %d events, %.0f weighted in theta_cm %g-%g\n", nEvents, wsum, cmMin, cmMax);
   printf("  Ex range with weight: %.2f to %.2f MeV   (Sn = 8.176)\n",
          h->GetBinCenter(h->FindFirstBinAbove(0)), h->GetBinCenter(h->FindLastBinAbove(0)));
   printf("  wrote %s/%s\n", dir.Data(), out.Data());
}
