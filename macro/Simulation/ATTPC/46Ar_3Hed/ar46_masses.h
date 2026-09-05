/// @file ar46_masses.h
/// @brief THE masses and beam constants of 46Ar(3He,d)47K. One definition, no ROOT dependencies.
///
/// This file exists because hand-converting these to amu for the browser explorer got them wrong
/// on 2026-09-05: the beam was out by 2.66 MeV and the residual by 1.60, so the page computed
/// Q = +3.521 MeV instead of +7.777 and every excitation energy it displayed was shifted by
/// -4.26 MeV. The macros were fine -- they use the MeV values below directly -- so nothing
/// disagreed except the one place a human had retyped the numbers.
///
/// RULE: never retype these anywhere. Include this header and divide by kU if amu are needed.
/// It deliberately includes nothing, so both the physics macros (which need AtMCTrack and friends)
/// and the explorer generator (which needs none of that) can take it.
#ifndef AR46_MASSES_H
#define AR46_MASSES_H

namespace Ar46 {

constexpr double kU = 931.49401; ///< MeV per u

// NUCLEAR masses in MeV: 46Ar + 3He -> d + 47K.  Q(g.s.) = kMb + kMt - kMe - kMR = +7.777 MeV,
// which reproduces the proposal's value and the two-body kinematics checked on 2026-08-11
// (theta_cm 15 deg -> theta_lab 131.4 deg, T_d 5.1 MeV; 80 deg -> 59.8 deg, 55.3 MeV).
constexpr double kMb = 42809.757; ///< 46Ar
constexpr double kMt = 2808.392;  ///< 3He
constexpr double kMR = 43734.759; ///< 47K
constexpr double kMe = 1875.613;  ///< d

/// Beam KE at z = 0 and its loss per cm of drift. The 95.7 MeV lost across the metre is why any
/// viewer of this reaction has to take the beam energy at the VERTEX rather than as a constant.
constexpr double kTb0 = 598.0;
constexpr double kdEdz = 0.957;

constexpr double kQgs = kMb + kMt - kMe - kMR;

} // namespace Ar46
#endif
