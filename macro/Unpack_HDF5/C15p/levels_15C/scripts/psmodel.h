// Three-body phase space for  15C + d -> d' + (14C + n).
//
// We detect the deuteron and reconstruct the missing mass M of the recoil system; Ex = M - m(15C).
// A RESONANCE in 15C* decays to 14C+n afterwards, so its Ex is still the true excitation and it
// appears as a peak. DIRECT breakup with no intermediate resonance populates the invariant mass
// over phase space instead -- that is what should describe the continuum under the peaks.
//
//   dN/dM  proportional to  p3*(M) * q(M)
//     p3*(M) = omega(s, m3^2, M^2) / (2 sqrt(s))     deuteron momentum in the overall CM
//     q(M)   = omega(M^2, m14^2, mn^2) / (2 M)       relative momentum inside the 14C+n system
//
// q ~ sqrt(Ex - Sn) at threshold (the sqrt rise); p3* -> 0 at the kinematic limit (the turnover).
#include <cmath>
namespace PS {
  const double U = 931.49401;
  const double M15C = 15.0105993*U, MD = 2.0141018*U;
  const double M14C = 14.0032420*U, MN  = 1.0086649*U;
  const double SN   = (M14C + MN) - M15C;              // 15C neutron separation energy, MeV
  inline double omega(double x,double y,double z){
    double a = x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z;
    return a>0 ? std::sqrt(a) : 0.0; }
  // s from the entrance channel: 15C beam of kinetic energy Kb on a deuteron at rest
  inline double sOf(double Kb){ return M15C*M15C + MD*MD + 2*MD*(Kb + M15C); }
  // phase-space density in Ex (unnormalised)
  inline double dens(double Ex,double Kb){
    double M = M15C + Ex;
    if (Ex <= SN) return 0.0;                          // below the 14C+n threshold
    double s = sOf(Kb), rs = std::sqrt(s);
    if (rs <= MD + M) return 0.0;                      // above the kinematic limit
    double p3 = omega(s, MD*MD, M*M)/(2*rs);
    double q  = omega(M*M, M14C*M14C, MN*MN)/(2*M);
    return p3*q; }
}
