#include "OpticalPotentialLibrary.h"
#include <cstdio>
#include <string>
int main(int c,char**v){double E=atof(v[1]); for(std::string n:{"A","H","Z","D","C"}){auto p=callPotential(n,16,6,E,1);
printf("%s ok%d V %.3f r %.4f a %.4f Wv %.3f rw %.4f aw %.4f Ws %.3f rs %.4f as %.4f Vso %.3f rso %.4f aso %.4f rc %.4f\n",n.c_str(),(int)p.ok,p.v,p.r0,p.a,p.vi,p.ri0,p.ai,p.vsi,p.rsi0,p.asi,p.vso,p.rso0,p.aso,p.rc0);}}
