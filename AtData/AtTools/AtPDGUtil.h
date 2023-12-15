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
int PdgToZ(int pdgc);
int PdgToA(int pdgc);
int IonPdgCode(int A, int Z, int L = 0, int I = 0);

bool IsProton(int pdgc);
bool IsNeutron(int pdgc);

const int kElectron = 11;  
const int kPositron = -11; 

const int kProton = 2212;
const int kNeutron = 2112; 
const int kPdgDeuteron = 1000010020;

} // namespace pdg
} // namespace AtTools

#endif // _PDG_UTILS_H_