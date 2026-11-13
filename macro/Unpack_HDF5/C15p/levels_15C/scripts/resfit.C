#pragma once
#include <Math/SpecFuncMathMore.h>
#include "binacc.h"
#include "psmodel.h"
#include "binbuild.C"
#include "lineshape.C"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
// =============================================================================================
//  RESONANCE-REGION fit of 15C(d,d'), Ex 3-9 MeV, 5 deg theta_CM bins out to 60 deg.
//
//  Three changes from the previous version, all of them simplifications:
//
//  1. **Ex restricted to 3-9 MeV.** The ground state, the 0.740 MeV 5/2+, the continuum turn-on
//     at Sn and the low-energy background turn-on at 2.5 are all BELOW the window. What remains
//     of them is a smooth tail, so the whole low-Ex complication -- which is where the
//     background/BW1 degeneracy lived -- is simply out of the fit.
//  2. **5 deg angular bins out to 60 deg**: 10 points instead of 5. 84 % of the data sits in
//     theta_CM 10-30 and was being squashed into two coarse bins; the channel closes near 65 deg
//     so there is nothing beyond.
//  3. **Only the resonances the data populates.** Forward selection decides; a resonance seen in
//     (p,p') need not be populated in (d,d'), and carrying one that is not just lets it absorb
//     whatever the continuum misses.
//
//  Model per bin:   smooth baseline  +  selected resonances
//     baseline(Ex) = C * dens(Ex)*acc(Ex)  +  Bg * exp(-kappa*Ex)     (kappa SHARED)
//     the second term is the tail of the low-energy background, whose turn-on is below the window
//  Resonances: barrier BW (Lane-Thomas) at the (p,p') energies, convolved with that bin's own
//  measured response.
// =============================================================================================
namespace RB {
   const double SN = 1.218, AMU = 931.494, HBARC = 197.327;
   const double R_FM = 1.4 * std::pow(15.0, 1.0/3.0), MU = AMU * 14.0 / 15.0;
   inline double P_ell(double Ed, int l)
   { if (Ed <= 0) return 0;
     double r = std::sqrt(2*MU*Ed)/HBARC*R_FM;
     double F = r*ROOT::Math::sph_bessel(l,r), G = -r*ROOT::Math::sph_neumann(l,r);
     double d = F*F+G*G; return d>0 ? r/d : 0; }
   inline double BW(double E,double Er,double Gr,int l)
   { double Ed=E-SN, Edr=Er-SN; if(Ed<=0||Edr<=0) return 0;
     double Pr=P_ell(Edr,l); if(Pr<=0) return 0;
     double h=0.5*Gr*P_ell(Ed,l)/Pr;
     return (h/TMath::Pi())/((E-Er)*(E-Er)+h*h); }
}
struct RT {
   std::vector<double> y; double x0=0,dx=1; int n=0;
   void build(double Er,double Gr,int l,double sig,double lo,double hi,int N=1200)
   { n=N; x0=lo; dx=(hi-lo)/(N-1); y.assign(N,0.0);
     const double dE=0.002,E0=RB::SN+1e-6,E1=23.0;
     int M=(int)((E1-E0)/dE)+1; std::vector<double> bw(M);
     for(int j=0;j<M;j++) bw[j]=RB::BW(E0+j*dE,Er,Gr,l);
     double cut=5*sig,nr=1.0/(std::sqrt(2*TMath::Pi())*sig);
     for(int i=0;i<N;i++){ double x=x0+i*dx,s=0;
       int jl=std::max(0,(int)((x-cut-E0)/dE)),jh=std::min(M-1,(int)((x+cut-E0)/dE));
       for(int j=jl;j<=jh;j++){ double u=(x-(E0+j*dE))/sig; s+=bw[j]*std::exp(-0.5*u*u); }
       y[i]=s*dE*nr; } }
   double Eval(double x) const
   { if(n<2) return 0; double t=(x-x0)/dx; int i=(int)t;
     if(i<0||i>=n-1) return 0; double f=t-i; return y[i]*(1-f)+y[i+1]*f; }
};
namespace RF {
   const int NR = 6;
   const double EPP6[NR] = {4.55,5.05,6.10,6.80,7.50,8.40};
   const double GPP6[NR] = {0.020,0.020,0.300,0.300,0.400,0.500};
   const int    LPP6[NR] = {1,1,2,2,2,2};
   const double PPGS6 = 0.119;
   double exLo = 3.0, exHi = 9.0;
   bool USE[NR] = {true,true,true,true,true,true};
   int  NUSE = 6;
   // TIE the 4.43/4.93 doublet: 0.5 MeV apart against a response FWHM of ~0.9, so they are not
   // separable and they visibly trade (ATTPCROOT 18 deg: 1.000/0.443, 22 deg: 0.016/1.000).
   // One amplitude for the pair with a shared internal ratio. 6.68 is left free -- its angular
   // distribution is already clean and reproducible between codes.
   bool TIE12 = false;
   struct Bin { TH1D* h=nullptr; std::vector<RT> T; TGraph* gc=nullptr;
                double gsTail=0, gsMu=0, sgL=0, sgR=0, alL=0, alR=0, gsAmp=0, e5Amp=0;
                double N=0, cmLo=0, cmHi=0; bool in[NR]; };
   std::vector<Bin> B;
   // reference response, taken from the highest-statistics bin, for bins whose own fit fails
   bool refValid=false; double refN=0, refSgL=0, refSgR=0, refAlL=0, refAlR=0, refMu=0;
   int nper() { return 2 + NUSE - (TIE12 ? 1 : 0); }   // C, Bg, then the amplitudes
   double model(int b,double X,const double* p)
   {
      const Bin& d = B[b];
      const double* q = p + 1 + b*nper();   // p[0] = shared kappa
      double v = std::max(0.0,q[0]) * (d.gc ? std::max(0.0,d.gc->Eval(X)) : 0.0)
               + std::max(0.0,q[1]) * std::exp(-std::max(0.0,p[0]) * X);
      // the g.s. and 5/2+ tails, FIXED from this bin's own peak fit
      v += d.gsAmp * shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
         + d.e5Amp * shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
      if(TIE12){
         // q[2] is the doublet amplitude; the LAST parameter overall is the internal ratio
         double r = std::max(0.0, p[1 + (int)B.size()*nper()]);
         double a = std::max(0.0,q[2]);
         if(d.in[0]) v += a * d.T[0].Eval(X);
         if(d.in[1]) v += a * r * d.T[1].Eval(X);
         int j=1;
         for(int i=2;i<NR;i++){ if(!USE[i]) continue;
           if(d.in[i]) v += std::max(0.0,q[2+j]) * d.T[i].Eval(X);
           j++; }
      } else {
         int j=0;
         for(int i=0;i<NR;i++){ if(!USE[i]) continue;
           if(d.in[i]) v += std::max(0.0,q[2+j]) * d.T[i].Eval(X);
           j++; }
      }
      return v;
   }
   double fcn(const double* p)
   {
      double s=0;
      for(int b=0;b<(int)B.size();b++){ TH1D* h=B[b].h; if(!h) continue;
        for(int i=1;i<=h->GetNbinsX();i++){ double x=h->GetBinCenter(i);
          if(x<exLo||x>exHi) continue;
          double n=h->GetBinContent(i), m=model(b,x,p);
          if(m<1e-9) m=1e-9;
          s += 2.0*(m-n+(n>0?n*std::log(n/m):0.0)); } }
      return s;
   }
}

// Prepare one bin: measure the response and the g.s./5/2+ amplitudes on the FULL spectrum (where
// those peaks live), then freeze them; build the resonance templates and the continuum shape.
void rfPrep(const Cfg& k,int ib)
{
   RF::Bin d;
   TFile f("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
   TH1D* h=(TH1D*)f.Get(Form("%s_cm%d",k.nm,ib));
   if(!h){ RF::B.push_back(d); return; }
   d.h=(TH1D*)h->Clone(Form("rf_%s_%d",k.nm,ib)); d.h->SetDirectory(nullptr);
   d.N=d.h->GetEntries(); d.cmLo=TCM0+ib*TCMSTEP; d.cmHi=d.cmLo+TCMSTEP;
   double s0=k.fw/2.355, pk=d.h->GetBinContent(d.h->GetMaximumBin());
   TF1 gt(Form("rfgt_%s_%d",k.nm,ib),f_GT,-2.0,2.4,8);
   gt.SetParameters(pk,0,s0,s0,1.5,1.5,5,0);
   gt.SetParLimits(1,-0.8,0.8); gt.SetParLimits(2,0.05,2.0); gt.SetParLimits(3,0.05,2.0);
   gt.SetParLimits(4,0.3,6.0);  gt.SetParLimits(5,0.3,6.0);
   d.h->Fit(&gt,"RQ0NL");
   d.gsAmp=gt.GetParameter(0); d.gsMu=gt.GetParameter(1);
   d.sgL=gt.GetParameter(2);   d.sgR=gt.GetParameter(3);
   d.alL=gt.GetParameter(4);   d.alR=gt.GetParameter(5);
   d.e5Amp = 0.08*d.gsAmp;
   // ⚠ SANITY CHECK ON THE RESPONSE FIT.  The elastic peak is only populated at FORWARD
   // theta_CM; in the large-angle bins there is essentially no ground state, the fit over
   // [-2, 2.4] has nothing to latch onto, and it returns a huge low shape -- ATTPCROOT bin 8 gave
   // FWHM_R = 4.69 MeV and Spyral bin 9 4.71, whose "tail" then contributed 22 and 26 counts at
   // Ex = 3 where the data has 6 and 1.  That alone produced chi2/ndf of 9.4 and 6.5.
   // If the fit is not credible, drop the g.s. contribution and borrow the response from the
   // reference bin (the one with the most counts), which is what the convolution needs anyway.
   bool ok = (2.355*d.sgR < 2.0) && (2.355*d.sgL < 2.0) && (d.gsAmp > 0.02*pk);
   if(!ok){
      d.gsAmp = 0; d.e5Amp = 0;
      if(RF::refValid){ d.sgL=RF::refSgL; d.sgR=RF::refSgR; d.alL=RF::refAlL; d.alR=RF::refAlR;
                        d.gsMu=RF::refMu; }
      printf("   [%s bin %d] response fit rejected (FWHM_R %.2f); g.s. dropped, response borrowed\n",
             k.nm,ib,2.355*gt.GetParameter(3));
   } else if(!RF::refValid || d.N > RF::refN){
      RF::refValid=true; RF::refN=d.N; RF::refSgL=d.sgL; RF::refSgR=d.sgR;
      RF::refAlL=d.alL; RF::refAlR=d.alR; RF::refMu=d.gsMu;
   }
   d.T.resize(RF::NR);
   for(int i=0;i<RF::NR;i++){ double mu=(RF::EPP6[i]-RF::PPGS6)+d.gsMu;
     d.in[i]=(mu>RF::exLo+d.sgR)&&(mu<RF::exHi-d.sgR);
     d.T[i].build(mu,RF::GPP6[i],RF::LPP6[i],d.sgR,RF::exLo-1.5,23.0); }
   Pars A=parsOf(k); Win w{14,86,k.keLo,k.keHi};
   d.gc=new TGraph();
   for(int i=0;i<=460;i++){ double X=i*0.05;
     d.gc->SetPoint(i,X,PS::dens(X,k.eb)*accBin(A,k.eb,X,d.cmLo,d.cmHi,w)); }
   RF::B.push_back(d);
}

double rfRun(const Cfg& k,int nbin,bool talk,std::vector<std::vector<double>>* out=nullptr)
{
   RF::B.clear();
   // two passes: the first establishes the reference response from the richest bin, the second
   // lets the failing bins borrow it
   RF::refValid=false; RF::refN=0;
   for(int ib=0;ib<nbin;ib++) rfPrep(k,ib);
   RF::B.clear();
   for(int ib=0;ib<nbin;ib++) rfPrep(k,ib);
   RF::NUSE=0; for(int i=0;i<RF::NR;i++) if(RF::USE[i]) RF::NUSE++;
   int npar=1+nbin*RF::nper()+(RF::TIE12?1:0);
   auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
   m->SetMaxFunctionCalls(2000000); m->SetMaxIterations(400000);
   m->SetTolerance(0.1); m->SetPrintLevel(0); m->SetStrategy(1);
   ROOT::Math::Functor fnc(&RF::fcn,npar); m->SetFunction(fnc);
   m->SetLimitedVariable(0,"kappa",0.25,0.005,0.0,2.0);
   for(int b=0;b<nbin;b++){ RF::Bin& d=RF::B[b]; if(!d.h) continue;
     double lvl=0;int nl=0;
     for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i);
       if(x>8&&x<9){ lvl+=d.h->GetBinContent(i); nl++; } }
     lvl=nl?lvl/nl:1; if(lvl<1) lvl=1;
     int o=1+b*RF::nper();
     double c9=d.gc?d.gc->Eval(6.0):1;
     m->SetLimitedVariable(o+0,Form("C_%d",b), c9>0?0.3*lvl/c9:1e-3, 1e-5, 0, 1e5);
     m->SetLimitedVariable(o+1,Form("Bg_%d",b), lvl*2, lvl*1e-2, 0, lvl*500);
     int j=0;
     for(int i=(RF::TIE12?1:0);i<RF::NR;i++){
       if(RF::TIE12 && i==1){ // the doublet slot, seeded on 4.43
         double sc=d.in[0]?d.T[0].Eval((RF::EPP6[0]-RF::PPGS6)+d.gsMu):0;
         double seed=(sc>0)?0.3*lvl/sc:0;
         m->SetLimitedVariable(o+2+j,Form("D12_%d",b),seed,std::max(1e-6,seed*1e-2),0,seed*300+1e3);
         j++; continue; }
       if(!RF::USE[i]) continue;
       double sc=d.in[i]?d.T[i].Eval((RF::EPP6[i]-RF::PPGS6)+d.gsMu):0;
       double seed=(sc>0)?0.3*lvl/sc:0;
       if(!d.in[i]) m->SetFixedVariable(o+2+j,Form("a%d_%d",i+1,b),0);
       else m->SetLimitedVariable(o+2+j,Form("a%d_%d",i+1,b),seed,std::max(1e-6,seed*1e-2),0,seed*300+1e3);
       j++; } }
   if(RF::TIE12) m->SetLimitedVariable(1+nbin*RF::nper(),"r_4.93/4.43",1.0,0.02,0.0,10.0);
   m->Minimize();
   const double* p=m->X();
   // Extract the per-bin amplitudes ALWAYS -- this used to live inside the `if(talk)` block, so
   // calling with talk = false returned an empty `out` and any figure built from it came out blank.
   for(int b=0;b<nbin;b++){ RF::Bin& d=RF::B[b]; if(!d.h) continue;
     const double* q=p+1+b*RF::nper();
     std::vector<double> row; int j=0;
     if(RF::TIE12){ row.push_back(d.in[0]?q[2]*d.T[0].Eval((RF::EPP6[0]-RF::PPGS6)+d.gsMu):0); j=1; }
     for(int i=(RF::TIE12?2:0);i<RF::NR;i++){ if(!RF::USE[i]) continue;
       row.push_back(d.in[i]?q[2+j]*d.T[i].Eval((RF::EPP6[i]-RF::PPGS6)+d.gsMu):0); j++; }
     if(out) out->push_back(row); }
   if(talk){
     printf("\n================ %s   Ex [%.1f,%.1f], %d x %.0f deg bins, %d resonances ================\n",
            k.nm,RF::exLo,RF::exHi,nbin,TCMSTEP,RF::NUSE);
     printf("   shared kappa = %.4f /MeV     -2lnL = %.1f   (%d parameters)\n",p[0],m->MinValue(),npar);
     printf("   bin  theta_CM     N   chi2/ndf |");
     if(RF::TIE12) printf("   D12  ");
     for(int i=(RF::TIE12?2:0);i<RF::NR;i++) if(RF::USE[i]) printf("  %5.2f ",RF::EPP6[i]-RF::PPGS6);
     printf("\n");
     for(int b=0;b<nbin;b++){ RF::Bin& d=RF::B[b]; if(!d.h) continue;
       double c2=0;int nd=0;
       for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
         if(x<RF::exLo||x>RF::exHi||n<5) continue;
         double mm=RF::model(b,x,p); c2+=(n-mm)*(n-mm)/n; nd++; }
       printf("   %2d   %2.0f-%2.0f  %7.0f  %6.2f  |",b,d.cmLo,d.cmHi,d.N,c2/std::max(1,nd-RF::nper()));
       if(out) for(double x:(*out)[b]) printf(" %6.1f ",x);
       printf("\n"); }
   }
   double v=m->MinValue(); delete m; return v;
}
void resfit(){ setCM(10,5,60); for(auto&k:CFG) rfRun(k,10,true); }
