#include "psmodel.h"
// What the continuum actually looks like, with no fit model in the way.
// Background-only windows (the states A and B sit at 4.5-4.7 and 6.6-6.9, so skip 3.5-8.0).
void whynops(){
  TFile f("dd_ex.root");
  double mu = PS::M14C*PS::MN/(PS::M14C+PS::MN);              // 14C-n reduced mass, MeV
  double kap = std::sqrt(2*mu*PS::SN);                        // halo momentum scale, MeV/c
  printf("Sn = %.4f MeV   reduced mass %.1f MeV   kappa = sqrt(2*mu*Sn) = %.1f MeV/c\n", PS::SN, mu, kap);
  printf("  -> the 15C g.s. neutron momentum distribution has width %.0f MeV/c,\n", kap);
  printf("     i.e. a 14C-n relative energy scale of kappa^2/(2mu) = %.2f MeV (= Sn, by construction)\n\n", kap*kap/(2*mu));
  for(const char* nm : {"ATTPCROOT","Spyral"}){
    TH1D* h=(TH1D*)f.Get(nm);
    auto avg=[&](double a,double b){ double s=0; int n=0;
      for(int i=1;i<=h->GetNbinsX();i++){ double x=h->GetBinCenter(i);
        if(x>a&&x<b){s+=h->GetBinContent(i);n++;} } return n?s/n:0.0; };
    double sub=avg(-2.5,-1.5);
    printf("== %s   (counts per 100 keV bin)\n", nm);
    printf("   below threshold (Ex -2.5..-1.5): %6.1f   <- instrumental, no 15C decay can live here\n", sub);
    struct {double a,b;} W[]={{2.0,3.4},{8.0,10},{10,12},{12,14},{14,16},{16,18},{18,20}};
    printf("   Ex window    data   uniform PS   PS x halo FF   (both normalised at Ex 8-10)\n");
    double d0=avg(8,10);
    double psN=0, hfN=0;
    auto psAt=[&](double X){ return PS::dens(X,191.25); };
    auto hfAt=[&](double X){ double M=PS::M15C+X;
      if(X<=PS::SN) return 0.0;
      double q=PS::omega(M*M,PS::M14C*PS::M14C,PS::MN*PS::MN)/(2*M);
      double s=PS::sOf(191.25), rs=std::sqrt(s);
      double p3=PS::omega(s,PS::MD*PS::MD,M*M)/(2*rs);
      return p3*q/std::pow(q*q+kap*kap,2); };
    { double sp=0,sh=0; int n=0;
      for(double X=8;X<10;X+=0.1){ sp+=psAt(X); sh+=hfAt(X); n++; }
      psN=d0/(sp/n); hfN=d0/(sh/n); }
    for(auto&w:W){ double d=avg(w.a,w.b), sp=0, sh=0; int n=0;
      for(double X=w.a;X<w.b;X+=0.1){ sp+=psAt(X); sh+=hfAt(X); n++; }
      printf("   %5.1f-%-5.1f %6.1f   %8.1f     %8.1f\n", w.a,w.b,d, psN*sp/n, hfN*sh/n); }
    printf("\n");
  }
}
