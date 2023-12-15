#include "AtPDGUtil.h"

namespace AtTools::pdg {

bool IsPseudoParticle(int pdgc)
{
   return pdgc;
}

bool IsIon(int pdgc)
{
   return (pdgc > 1000000000 && pdgc < 1999999999) || pdgc == kProton;
}

std::pair<int, int> PdgToIon(int pdgc)
{
   return std::pair<int, int>(PdgToA(pdgc), PdgToZ(pdgc));
}

int PdgToZ(int pdgc)
{
   if (pdgc == kProton)
      return 1;
   if (pdgc == kNeutron)
      return 0;
   int Z = (pdgc / 10000) - 1000 * (pdgc / 10000000); // don't factor out!
   return Z;
}

int PdgToA(int pdgc)
{
   if (pdgc == kProton || pdgc == kNeutron)
      return 1;
   int A = (pdgc / 10) - 1000 * (pdgc / 10000); // don't factor out!
   return A;
}

int IonPdgCode(int A, int Z, int L, int I)
{
   if (A == 1)
      if (Z == 1)
         return kProton;
      else
         return kNeutron;

   return 1000000000 + L * 100000000 + Z * 10000 + A * 10 + I;
}

bool IsProton(int pdgc)
{
   return pdgc == kProton || pdgc == IonPdgCode(1, 1);
}
bool IsNeutron(int pdgc)
{
   return pdgc == kNeutron || pdgc == IonPdgCode(1, 0);
}
} // namespace AtTools::pdg
