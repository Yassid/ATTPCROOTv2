/// @file tp_kinematics_Be10.C
/// @brief The two-body kinematics of 10Be(t,p)12Be, and the LEVERAGE factors that decide whether a
///        better detector would help -- computed analytically, with 14C(d,p)15C alongside.
///
///   root -b -q 'tp_kinematics_Be10.C()'
///
/// This needs no simulation: it is exact relativistic two-body kinematics. It exists because the
/// (d,p) campaign's whole framing was that dEx/dKE and dEx/dE_beam -- how much an error in the
/// measured proton energy, or in the beam energy, costs in excitation energy -- decide what a
/// field or pad upgrade can buy. Those are properties of the reaction, not of the detector, so
/// they should be known before a single event is generated.

static double kn_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// Ex from (K_beam, theta_lab, K_ejectile) -- the expression the analysis inverts
static double kn_ex(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(th) * kn_om2(s, m1 * m1, m2 * m2) * kn_om2(uu, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) / (2 * m2 * m2) + s + uu - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

/// forward kinematics: theta_cm -> (theta_lab, KE) of the ejectile, for residual mass m4.
///
/// CONVENTION, and it is not a detail. `thcm` here is the ANALYSIS/DWBA angle, the one
/// acceptance_C14.C and every angular distribution use: the angle between the outgoing light
/// particle and the INCOMING light particle, i.e. 180 deg minus the ejectile's CM angle measured
/// from the beam. In inverse kinematics that flip is the whole story -- SMALL theta_cm, where a
/// transfer angular distribution carries its yield, is the BACKWARD-lab, LOW-energy proton, not
/// the forward 50 MeV one. Getting it the wrong way round makes the channel look easy.
static bool kn_fwd(double m1, double m2, double m3, double m4, double K, double thcm_analysis, double &thlab,
                   double &Ke)
{
   const double thcm = TMath::Pi() - thcm_analysis;
   double E1 = K + m1;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * E1;
   double rs = std::sqrt(s);
   if (rs < m3 + m4) return false;
   double pcm = kn_om2(s, m3 * m3, m4 * m4) / (2 * rs);
   double Ecm3 = std::sqrt(pcm * pcm + m3 * m3);
   // boost of the CM in the lab
   double plab = std::sqrt(E1 * E1 - m1 * m1);
   double beta = plab / (E1 + m2);
   double gam = 1.0 / std::sqrt(1 - beta * beta);
   double pz = gam * (pcm * std::cos(thcm) + beta * Ecm3);
   double pt = pcm * std::sin(thcm);
   double E3 = gam * (Ecm3 + beta * pcm * std::cos(thcm));
   Ke = E3 - m3;
   thlab = std::atan2(pt, pz);
   return true;
}

void tp_kinematics_Be10(Double_t Ebeam = 112.20)
{
   const double U = 931.49401;
   struct Chan { const char *name; double m1, m2, m3, m4, E; double ex[4]; int nex; };
   Chan ch[2] = {
      {"10Be(t,p)12Be", 10.0135341 * U, 3.0160493 * U, 1.007825 * U, 12.0269221 * U, Ebeam,
       {0.0, 2.109, 2.251, 2.715}, 4},
      {"14C(d,p)15C  ", 14.003242 * U, 2.0141018 * U, 1.007825 * U, 15.0105993 * U, 159.75,
       {0.0, 0.740, 0.0, 0.0}, 2},
   };

   for (int c = 0; c < 2; ++c) {
      auto &k = ch[c];
      double Q = k.m1 + k.m2 - k.m3 - k.m4;
      printf("\n================ %s  at %.2f MeV (%.2f MeV/u) ================\n", k.name, k.E,
             k.E / std::round(k.m1 / U));
      printf("Q(g.s.) = %+.3f MeV\n", Q);
      printf("\n  theta_cm is the ANALYSIS/DWBA angle: small theta_cm = backward lab (see kn_fwd).\n");
      printf("\n  %8s %10s %10s | %10s %12s %14s\n", "theta_cm", "theta_lab", "KE_p [MeV]", "dEx/dKE",
             "dEx/dE_beam", "dEx/dtheta[/deg]");
      for (double thcm = 5; thcm <= 175; thcm += (thcm < 60 ? 5 : 15)) {
         double thlab, Ke;
         if (!kn_fwd(k.m1, k.m2, k.m3, k.m4, k.E, thcm * TMath::DegToRad(), thlab, Ke)) continue;
         if (Ke <= 0.2) continue;
         double h = 0.01;
         double dKE = (kn_ex(k.m1, k.m2, k.m3, k.m4, k.E, thlab, Ke + h) -
                       kn_ex(k.m1, k.m2, k.m3, k.m4, k.E, thlab, Ke - h)) / (2 * h);
         double dEb = (kn_ex(k.m1, k.m2, k.m3, k.m4, k.E + h, thlab, Ke) -
                       kn_ex(k.m1, k.m2, k.m3, k.m4, k.E - h, thlab, Ke)) / (2 * h);
         double ht = 0.001;
         double dTh = (kn_ex(k.m1, k.m2, k.m3, k.m4, k.E, thlab + ht, Ke) -
                       kn_ex(k.m1, k.m2, k.m3, k.m4, k.E, thlab - ht, Ke)) / (2 * ht) * TMath::DegToRad();
         printf("  %8.0f %10.1f %10.2f | %10.3f %12.4f %14.4f\n", thcm, thlab * TMath::RadToDeg(), Ke, std::fabs(dKE),
                std::fabs(dEb), std::fabs(dTh));
      }
      // where each level lands
      printf("\n  level-by-level lab angle and proton energy, at theta_cm = 10 / 30 / 60 deg:\n");
      for (int l = 0; l < k.nex; ++l) {
         printf("    Ex = %.3f :", k.ex[l]);
         for (double thcm : {10.0, 30.0, 60.0}) {
            double thlab, Ke;
            if (kn_fwd(k.m1, k.m2, k.m3, k.m4 + k.ex[l], k.E, thcm * TMath::DegToRad(), thlab, Ke))
               printf("   %6.1f deg / %6.2f MeV", thlab * TMath::RadToDeg(), Ke);
            else
               printf("   %21s", "closed");
         }
         printf("\n");
      }
   }
   printf("\nREADING THIS. dEx/dKE is how many MeV of excitation-energy error one MeV of proton\n"
          "energy error buys. A LARGE value means the measurement is hard but also that improving\n"
          "the proton energy pays; a small one means the channel is insensitive either way.\n"
          "dEx/dE_beam is the same for the beam, and it is what the vertex-by-vertex beam-energy\n"
          "correction removes -- which is free, in software.\n");
   printf("\nkinematics done\n\n");
}
