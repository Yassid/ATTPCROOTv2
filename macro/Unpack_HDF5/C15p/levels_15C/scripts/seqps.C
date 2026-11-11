#include "psmodel.h"
#include "exlib.h"
// ============================================================================================
//  The continuum template done properly: SEQUENTIAL decay, with the detector acceptance folded in.
//
//  WHY THE OLD TEMPLATE WAS WRONG.  `psalt.C` used TGenPhaseSpace to make d + 14C + n as a
//  uniform, uncorrelated three-body final state.  That throws the deuteron FORWARD:
//
//     theta        0-10 10-20 20-30 30-40 40-50 50-60 60-70 70-80 80-90
//     uniform 3b    7.8  20.3  25.4  22.6  14.9   6.9   1.9   0.1   0.0  %
//     real data     0.0   1.5   2.8   3.4   4.3   6.5   7.2  35.4  38.8  %
//
//  The two distributions are disjoint.  But the physical process is SEQUENTIAL,
//      15C + d -> d' + 15C*(M),   then   15C* -> 14C + n,
//  so the deuteron has TWO-BODY kinematics against a recoil of mass M = M(15C) + Ex and goes
//  BACKWARD, into exactly the angular region where the data and the discrete states live.
//  Uniform three-body phase space is simply the wrong kinematics for this channel -- which is why
//  the fit kept throwing the template away.  It was never evidence against a 14C+n continuum.
//
//  Here: for each Ex, generate the two-body d' + 15C*(Ex) kinematics over CM angle, keep the
//  deuterons that land inside the detector's (theta_lab, KE) acceptance, and weight by the
//  14C+n density of states.  acceptance(Ex) x dens(Ex) is the template.
// ============================================================================================
static double ACC_TH_LO = 40, ACC_TH_HI = 90, ACC_KE_LO = 0.5, ACC_KE_HI = 18.0;

// fraction of CM solid angle (isotropic) that lands inside the acceptance, for a recoil of mass M
double acc_of(double Ex, double Kb, double* meanTh = nullptr, double* meanKE = nullptr)
{
   double M4 = PS::M15C + Ex, M3 = PS::MD;
   double Eb = Kb + PS::M15C, pb = std::sqrt(Kb * Kb + 2 * Kb * PS::M15C);
   TLorentzVector W(0, 0, pb, Eb + PS::MD);
   double s = W.M2(), rs = std::sqrt(s);
   if (rs <= M3 + M4) return 0.0;
   double pcm = PS::omega(s, M3 * M3, M4 * M4) / (2 * rs);
   double E3 = std::sqrt(pcm * pcm + M3 * M3);
   TVector3 b = W.BoostVector();
   const int N = 2000; double acc = 0, tot = 0, sTh = 0, sKE = 0;
   for (int i = 0; i < N; i++) {
      double c = -1.0 + 2.0 * (i + 0.5) / N;             // uniform in cos(theta_cm)
      double st = std::sqrt(std::max(0.0, 1 - c * c));
      TLorentzVector d(pcm * st, 0, pcm * c, E3);
      d.Boost(b);
      double th = d.Theta() * 180.0 / M_PI, ke = d.E() - M3;
      tot += 1;
      if (th >= ACC_TH_LO && th <= ACC_TH_HI && ke >= ACC_KE_LO && ke <= ACC_KE_HI) {
         acc += 1; sTh += th; sKE += ke;
      }
   }
   if (meanTh) *meanTh = acc > 0 ? sTh / acc : 0;
   if (meanKE) *meanKE = acc > 0 ? sKE / acc : 0;
   return acc / tot;
}
void seqps()
{
   double Kb = 191.25;
   printf("Sequential 15C(d,d')15C*(Ex), 15C* -> 14C+n.  Acceptance = fraction of CM solid angle\n");
   printf("landing in theta_lab %.0f-%.0f deg and KE %.1f-%.1f MeV.\n\n", ACC_TH_LO, ACC_TH_HI, ACC_KE_LO, ACC_KE_HI);
   printf("   Ex     acc      <theta>  <KE>    dens(Ex)   dens*acc   (dens*acc normalised at 9)\n");
   double n9 = 0; { double a = acc_of(9.0, Kb); n9 = PS::dens(9.0, Kb) * a; }
   for (double X = 2; X <= 20.01; X += 2) {
      double th, ke, a = acc_of(X, Kb, &th, &ke);
      double d = PS::dens(X, Kb);
      printf("  %5.1f   %6.3f   %6.1f  %5.1f   %9.3g  %9.3g      %6.2f\n",
             X, a, th, ke, d, d * a, n9 > 0 ? d * a / n9 : 0);
   }
   // and what the DATA looks like, same normalisation
   TFile f("dd_ex.root");
   for (const char* nm : {"ATTPCROOT", "Spyral"}) {
      TH1D* h = (TH1D*)f.Get(nm);
      auto avg = [&](double a, double b) { double s = 0; int n = 0;
         for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i);
            if (x > a && x < b) { s += h->GetBinContent(i); n++; } } return n ? s / n : 0.0; };
      double d9 = avg(8, 10);
      printf("\n  %-10s data (normalised at 8-10):", nm);
      for (double X = 2; X <= 20.01; X += 2) printf("  %4.1f:%5.2f", X, avg(X - 1, X + 1) / d9);
      printf("\n");
   }
}
