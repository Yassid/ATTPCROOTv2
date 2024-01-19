#ifndef _PDG_UTILS_H_
#define _PDG_UTILS_H_

#include <utility>

namespace AtTools {

/**
 * Namespace to help manage PDG particle codes.
 * Note: We treat protons the same an H ions.
 */
namespace pdg {
bool IsPseudoParticle(int pdgc);
bool IsIon(int pdgc);

// Note: PDG codes for nuclear targets can be computed using pdg::IonPdgCode(A,Z)
// PDG2006 convention: 10LZZZAAAI
std::pair<int, int> PdgToIon(int pdgc); //< Returns <A,Z>.
constexpr int PdgToZ(int pdgc);
constexpr int PdgToA(int pdgc);
constexpr int IonPdgCode(int A, int Z, int L = 0, int I = 0);

constexpr bool IsProton(int pdgc);
constexpr bool IsNeutron(int pdgc);

constexpr int kElectron = 11;
constexpr int kPositron = -11;

constexpr int kProton = 2212;
constexpr int kNeutron = 2112;
constexpr int kDeuteron = 1*10000+2*10;
constexpr int kAlpha= 2*10000+4*10;


} // namespace pdg
} // namespace AtTools

#endif // _PDG_UTILS_H_