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
struct RTD {
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
namespace RK {
   // A SEVENTH slot, not in the (p,p') scheme.  The residuals of the richest bin show strength
   // at Ex ~ 3.0-3.7 that nothing in the model covers, in BOTH codes -- the long-standing
   // "unexplained strength at Ex ~ 3.4".  A state seen in (d,d') but not (p,p') is allowed.
   const int NR = 7;
   double EPP6[NR] = {4.55,5.05,6.10,6.80,7.50,8.40, 3.52};
   double GPP6[NR] = {0.020,0.020,0.300,0.300,0.400,0.500, 0.300};
   int    LPP6[NR] = {1,1,2,2,2,2, 2};
   const double PPGS6 = 0.119;
   double exLo = 3.0, exHi = 9.0;
   bool USE[NR] = {true,true,true,true,true,true,false};
   int  NUSE = 6;
   // TIE the 4.43/4.93 doublet: 0.5 MeV apart against a response FWHM of ~0.9, so they are not
   // separable and they visibly trade (ATTPCROOT 18 deg: 1.000/0.443, 22 deg: 0.016/1.000).
   // One amplitude for the pair with a shared internal ratio. 6.68 is left free -- its angular
   // distribution is already clean and reproducible between codes.
   bool TIE12 = false;
   // ---- GROUPS -------------------------------------------------------------------------
   // Two pairs in this spectrum are below resolution and their individual amplitudes are the
   // degenerate direction of the fit:  4.43/4.93 (0.50 MeV apart) and 5.98/6.68 (0.70 MeV),
   // against a forward-bin response FWHM of 0.93-1.11.  Tie each pair to ONE amplitude with a
   // FIXED internal weight.  Fitting that internal ratio, as the old TIE12 did, is fitting the
   // degenerate direction itself -- which is why the per-bin yields swung by factors of 50.
   // The deliverable is then the PAIR's yield, which is what the data actually determine.
   int    GRP[NR];        // resonance -> amplitude slot, -1 = not in the fit
   double WGT[NR];        // weight within its group when RIDX < 0
   int    RIDX[NR];       // -> GLOBAL ratio parameter for this member; -1 = weight fixed at 1
   int    NGRP = 0, NRAT = 0;
   bool   PAIRS = false;  // true: 0+1 share a slot, 2+3 share a slot
   // NPERBIN: give the CONTINUUM its shape parameter per bin instead of globally.  That costs
   // exactly the 10 parameters a new resonance amplitude costs, so the two can be compared
   // head to head.  If the continuum does the same job, the "state" was continuum freedom.
   bool   NPERBIN = false;
   // ONE global shift applied to every state in every bin.  If our Ex scale is offset -- which the
   // NNDC comparison suggests by ~0.2 -- then every resonance is misplaced by the same amount and
   // no amount of per-state tuning fixes it.  This is the single parameter that tests that.
   double DELTA = 0.0;
   std::vector<double> SEED;     // non-empty and right-sized: used as the starting point
   std::vector<double> LASTP;    // the converged parameter vector, for seeding the next fit
   int   ONLYBIN = -1;           // >=0: fit just this bin
   double FIXRAT[4] = {1,1,1,1}; // pair ratios held at these values when ONLYBIN >= 0
   void setGroups(){
      for(int i=0;i<NR;i++){ GRP[i]=-1; WGT[i]=1.0; RIDX[i]=-1; }
      NGRP=0; NRAT=0;
      if(PAIRS){
         // One GLOBAL ratio per pair, shared by every angular bin.  Free per bin it is the
         // degenerate direction (yields swung by factors of 50); fixed at 1:1 it is too stiff and
         // the richest bin stays badly fitted.  One shared ratio is the middle that works.
         if(USE[0]||USE[1]){ int g=NGRP++; if(USE[0])GRP[0]=g;
            if(USE[1]){ GRP[1]=g; if(USE[0]) RIDX[1]=NRAT++; } }
         if(USE[2]||USE[3]){ int g=NGRP++; if(USE[2])GRP[2]=g;
            if(USE[3]){ GRP[3]=g; if(USE[2]) RIDX[3]=NRAT++; } }
         for(int i=4;i<NR;i++) if(USE[i]) GRP[i]=NGRP++;
      } else for(int i=0;i<NR;i++) if(USE[i]) GRP[i]=NGRP++;
   }
   double LASTP0 = 0;
   double LASTA[32][10]={{0}}; double LASTBG[32]={0}; double LASTR=1.0;
   // Background FORM.  The second baseline term used to be the empirical exp(-kappa*Ex).
   // 15C is an s-wave one-neutron halo (Sn = 1.218 MeV), and DIRECT breakup -- which is not at
   // all the same thing as the uniform three-body phase space tested earlier -- has a standard
   // shape peaked just above threshold and falling steeply:
   //       s-wave   dsigma/dEx  ~  sqrt(Ex-Sb) / Ex^2
   //       p-wave   dsigma/dEx  ~  (Ex-Sb)^1.5 / Ex^3
   // p[0] is the single SHARED shape parameter: kappa for EXP/LIN, the binding energy Sb for the
   // breakup forms.  Floating Sb turns the fit into a measurement of the halo binding energy.
   // PS/HALO/HALOF need the beam energy to evaluate the three-body density.
   double EB = 191.25;
   // ONETERM: use the chosen form as the SOLE baseline, dropping the separate C*dens*acc term.
   // Without this the breakup forms compete against PS::dens, which near threshold has almost
   // the SAME sqrt(Ex-Sn) shape -- the two share strength and the breakup parameter runs away.
   bool ONETERM = false;
   enum { EXP=0, SWAVE=1, PWAVE=2, LIN=3, PS=4, HALO=5, HALOF=6, HALON=7 };
   int  BKG  = EXP;
   bool SBFIX = false; double SBVAL = 1.218;
   inline double bshape(double X,double s)
   {  switch(BKG){
        case SWAVE: { double u=X-s; return u>0 ? std::sqrt(u)/(X*X) : 0.0; }
        case PWAVE: { double u=X-s; return u>0 ? u*std::sqrt(u)/(X*X*X) : 0.0; }
        case LIN:   { double v=1.0-s*X; return v>0 ? v : 0.0; }
        // uniform three-body phase space, p3*(M) q(M): zero free shape parameters
        case PS:    return PS::dens(X,EB);
        // DIRECT s-wave halo breakup = phase space x the halo momentum distribution.
        // |phi(q)|^2 ~ 1/(q^2+kappa^2)^2 with q^2/2mu = Ex-Sn and kappa^2/2mu = Sn, so the
        // denominator is simply Ex^2.  ZERO free shape parameters.
        case HALO:  { double d=PS::dens(X,EB); return d>0 ? d/(X*X) : 0.0; }
        // same, with the halo binding energy FLOATED: if it returns Sn the model is confirmed
        // the halo form factor with its POWER floated.  An s-wave halo predicts exactly 2:
        // whatever the fit returns is a direct, interpretable measure of how far off we are.
        case HALON: { double d=PS::dens(X,EB); return d>0 ? d*std::pow(X,-s) : 0.0; }
        case HALOF: { double d=PS::dens(X,EB); if(d<=0) return 0.0;
                      double den=X-PS::SN+std::max(1e-3,s); return d/(den*den); }
        default:    return std::exp(-std::max(0.0,s)*X); } }
   struct Bin { TH1D* h=nullptr; std::vector<RTD> T; TGraph* gc=nullptr;
                double gsTail=0, gsMu=0, sgL=0, sgR=0, alL=0, alR=0, gsAmp=0, e5Amp=0;
                double N=0, cmLo=0, cmHi=0; bool in[NR]; TGraph* ga=nullptr; };
   std::vector<Bin> B;
   // reference response, taken from the highest-statistics bin, for bins whose own fit fails
   bool refValid=false; double refN=0, refSgL=0, refSgR=0, refAlL=0, refAlR=0, refMu=0;
   int nper() { return 2 + (PAIRS ? NGRP : NUSE - (TIE12 ? 1 : 0)) + (NPERBIN?1:0); }
   double model(int b,double X,const double* p)
   {
      const Bin& d = B[b];
      const double* q = p + 1 + b*nper();   // p[0] = shared kappa
      // ⚠ The continuum model is ABSOLUTE in Ex: its denominator is Ex = E_rel + Sn, referenced
      // to the 15C ground state.  The resonance templates below are already placed relative to
      // this bin's MEASURED g.s. centroid (mu = EPP6 - PPGS6 + gsMu), but the continuum was being
      // evaluated at the raw X.  ATTPCROOT's g.s. sits at +0.027 and Spyral's at -0.287, so the
      // model was evaluated 0.31 MeV apart in the two codes -- which for 1/Ex^n over Ex 3-9
      // mimics a change in n and is why they disagreed by 3.5 sigma on the exponent.
      const double Xt = X - d.gsMu;          // true excitation energy
      const double shp = NPERBIN ? q[2 + (PAIRS?NGRP:0) + (PAIRS?0:0)] : p[0];
      double v = (ONETERM ? 0.0 : std::max(0.0,q[0]) * (d.gc ? std::max(0.0,d.gc->Eval(Xt)) : 0.0))
               + std::max(0.0,q[1]) * bshape(Xt,shp)
                                      * (d.ga ? std::max(0.0,d.ga->Eval(Xt)) : 1.0);
      // the g.s. and 5/2+ tails, FIXED from this bin's own peak fit
      v += d.gsAmp * shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
         + d.e5Amp * shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
      if(PAIRS){
         const double* rb = (ONLYBIN>=0) ? FIXRAT : p + 1 + (int)B.size()*nper();
         for(int i=0;i<NR;i++) if(GRP[i]>=0 && d.in[i]){
            double w = (RIDX[i]>=0) ? std::max(0.0,rb[RIDX[i]]) : WGT[i];
            v += std::max(0.0,q[2+GRP[i]]) * w * d.T[i].Eval(X); }
      } else if(TIE12){
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
void rkPrep(const Cfg& k,int ib)
{
   RK::Bin d;
   TFile f("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
   TH1D* h=(TH1D*)f.Get(Form("%s_cm%d",k.nm,ib));
   if(!h){ RK::B.push_back(d); return; }
   d.h=(TH1D*)h->Clone(Form("rkd_%s_%d",k.nm,ib)); d.h->SetDirectory(nullptr);
   d.N=d.h->GetEntries(); d.cmLo=TCM0+ib*TCMSTEP; d.cmHi=d.cmLo+TCMSTEP;
   double s0=k.fw/2.355, pk=d.h->GetBinContent(d.h->GetMaximumBin());
   TF1 gt(Form("rkdgt_%s_%d",k.nm,ib),f_GT,-2.0,2.4,8);
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
   // ⚠ The amplitude test compares against the BIN'S OWN maximum, which in a bin with no
   // ground state is just noise -- so it passes trivially.  Spyral 50-55 has gsAmp = 8 in a bin
   // with ~4 counts per bin, and the peak fit returned gsMu = +0.515 against a true -0.28: a
   // centroid 0.8 MeV wrong, which then placed EVERY resonance template 0.8 MeV too high
   // (mu = EPP6 - PPGS6 + gsMu) and pushed 8.28 out of the window.  chi2R 920, EDM 1.5e+05.
   // Also require the centroid to agree with the reference bin, which IS well measured.
   bool ok = (2.355*d.sgR < 2.0) && (2.355*d.sgL < 2.0) && (d.gsAmp > 0.02*pk)
             && (!RK::refValid || std::fabs(d.gsMu - RK::refMu) < 0.30);
   if(!ok){
      d.gsAmp = 0; d.e5Amp = 0;
      if(RK::refValid){ d.sgL=RK::refSgL; d.sgR=RK::refSgR; d.alL=RK::refAlL; d.alR=RK::refAlR;
                        d.gsMu=RK::refMu; }
      printf("   [%s bin %d] response fit rejected (FWHM_R %.2f); g.s. dropped, response borrowed\n",
             k.nm,ib,2.355*gt.GetParameter(3));
   } else if(!RK::refValid || d.N > RK::refN){
      RK::refValid=true; RK::refN=d.N; RK::refSgL=d.sgL; RK::refSgR=d.sgR;
      RK::refAlL=d.alL; RK::refAlR=d.alR; RK::refMu=d.gsMu;
   }
   d.T.resize(RK::NR);
   for(int i=0;i<RK::NR;i++){ double mu=(RK::EPP6[i]-RK::PPGS6)+d.gsMu+RK::DELTA;
     d.in[i]=(mu>RK::exLo+d.sgR)&&(mu<RK::exHi-d.sgR);
     d.T[i].build(mu,RK::GPP6[i],RK::LPP6[i],d.sgR,RK::exLo-1.5,23.0); }
   RK::EB=k.eb;
   Pars A=parsOf(k); Win w{14,86,k.keLo,k.keHi};
   d.gc=new TGraph();
   for(int i=0;i<=460;i++){ double X=i*0.05;
     d.gc->SetPoint(i,X,PS::dens(X,k.eb)*accBin(A,k.eb,X,d.cmLo,d.cmHi,w)); }
   // the breakup term is a physical yield, so it must pass through the acceptance too
   d.ga=new TGraph();
   for(int i=0;i<=460;i++){ double X=i*0.05;
     d.ga->SetPoint(i,X,accBin(A,k.eb,X,d.cmLo,d.cmHi,w)); }
   RK::B.push_back(d);
}

double rkRun(const Cfg& k,int nbin,bool talk,std::vector<std::vector<double>>* out=nullptr)
{
   RK::B.clear();
   // two passes: the first establishes the reference response from the richest bin, the second
   // lets the failing bins borrow it
   RK::refValid=false; RK::refN=0;
   for(int ib=0;ib<nbin;ib++) rkPrep(k,ib);
   RK::B.clear();
   for(int ib=0;ib<nbin;ib++) rkPrep(k,ib);
   RK::NUSE=0; for(int i=0;i<RK::NR;i++) if(RK::USE[i]) RK::NUSE++;
   RK::setGroups();
   int npar=1+nbin*RK::nper()+(RK::PAIRS?RK::NRAT:((RK::TIE12)?1:0));
   const bool useSeed = ((int)RK::SEED.size()==npar);
   auto SD=[&](int i,double v){ return useSeed ? RK::SEED[i] : v; };
   auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
   m->SetMaxFunctionCalls(2000000); m->SetMaxIterations(400000);
   m->SetTolerance(0.1); m->SetPrintLevel(0); m->SetStrategy(1);
   ROOT::Math::Functor fnc(&RK::fcn,npar); m->SetFunction(fnc);
   if(RK::NPERBIN) m->SetFixedVariable(0,"(per-bin)",0.0);
   else if(RK::BKG==RK::PS||RK::BKG==RK::HALO) m->SetFixedVariable(0,"(none)",0.0);
   else if(RK::BKG==RK::HALOF) m->SetLimitedVariable(0,"Sb",1.0,0.01,0.02,6.0);
   else if(RK::BKG==RK::HALON){ if(RK::SBFIX) m->SetFixedVariable(0,"n",RK::SBVAL);
                                else m->SetLimitedVariable(0,"n",2.0,0.02,0.0,9.0); }
   else if(RK::BKG==RK::SWAVE||RK::BKG==RK::PWAVE){
      if(RK::SBFIX) m->SetFixedVariable(0,"Sb",RK::SBVAL);
      else          m->SetLimitedVariable(0,"Sb",1.0,0.01,-3.0,2.9);
   } else if(RK::BKG==RK::LIN) m->SetLimitedVariable(0,"slope",0.08,0.002,0.0,0.33);
   else                        m->SetLimitedVariable(0,"kappa",0.25,0.005,0.0,2.0);
   for(int b=0;b<nbin;b++){ RK::Bin& d=RK::B[b]; if(!d.h) continue;
     double lvl=0;int nl=0;
     for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i);
       if(x>8&&x<9){ lvl+=d.h->GetBinContent(i); nl++; } }
     lvl=nl?lvl/nl:1; if(lvl<1) lvl=1;
     int o=1+b*RK::nper();
     double c9=d.gc?d.gc->Eval(6.0):1;
     if(RK::ONETERM) m->SetFixedVariable(o+0,Form("C_%d",b),0.0);
     else m->SetLimitedVariable(o+0,Form("C_%d",b), c9>0?0.3*lvl/c9:1e-3, 1e-5, 0, 1e5);
     { double f6=RK::bshape(6.0,(RK::BKG==RK::HALOF)?1.218:((RK::BKG==RK::HALON)?2.0:0.25));
       double sd=(f6>0)? lvl/f6 : lvl*2;
       m->SetLimitedVariable(o+1,Form("Bg_%d",b), SD(o+1,sd), std::max(1e-9,sd*1e-2), 0, sd*1e4); }
     if(RK::NPERBIN){ int ip=o+2+(RK::PAIRS?RK::NGRP:RK::NUSE);
       m->SetLimitedVariable(ip,Form("n_%d",b),SD(ip,2.0),0.02,0.0,9.0); }
     if(RK::PAIRS){
       for(int g=0;g<RK::NGRP;g++){
         // seed from the first member of the group that this bin can see
         double sc=0; bool any=false;
         for(int i=0;i<RK::NR;i++) if(RK::GRP[i]==g && d.in[i]){
           sc=std::max(sc,d.T[i].Eval((RK::EPP6[i]-RK::PPGS6)+d.gsMu)); any=true; }
         if(!any){ m->SetFixedVariable(o+2+g,Form("G%d_%d",g,b),0); continue; }
         double seed=(sc>0)?0.3*lvl/sc:0;
         m->SetLimitedVariable(o+2+g,Form("G%d_%d",g,b),SD(o+2+g,seed),
                               std::max(1e-6,seed*1e-2),0,seed*300+1e3); }
       continue; }
     int j=0;
     for(int i=(RK::TIE12?1:0);i<RK::NR;i++){
       if(RK::TIE12 && i==1){ // the doublet slot, seeded on 4.43
         double sc=d.in[0]?d.T[0].Eval((RK::EPP6[0]-RK::PPGS6)+d.gsMu):0;
         double seed=(sc>0)?0.3*lvl/sc:0;
         m->SetLimitedVariable(o+2+j,Form("D12_%d",b),seed,std::max(1e-6,seed*1e-2),0,seed*300+1e3);
         j++; continue; }
       if(!RK::USE[i]) continue;
       double sc=d.in[i]?d.T[i].Eval((RK::EPP6[i]-RK::PPGS6)+d.gsMu):0;
       double seed=(sc>0)?0.3*lvl/sc:0;
       if(!d.in[i]) m->SetFixedVariable(o+2+j,Form("a%d_%d",i+1,b),0);
       else m->SetLimitedVariable(o+2+j,Form("a%d_%d",i+1,b),seed,std::max(1e-6,seed*1e-2),0,seed*300+1e3);
       j++; } }
   if(RK::PAIRS){ for(int r=0;r<RK::NRAT;r++){ int ip=1+nbin*RK::nper()+r;
        m->SetLimitedVariable(ip,Form("ratio%d",r),SD(ip,1.0),0.02,0.0,20.0); } }
   else if(RK::TIE12) m->SetLimitedVariable(1+nbin*RK::nper(),"r_4.93/4.43",1.0,0.02,0.0,10.0);
   m->Minimize();
   const double* p=m->X(); RK::LASTP0=p[0];
   RK::LASTP.assign(p,p+npar);
   RK::LASTR = (RK::TIE12 && !RK::PAIRS) ? std::max(0.0,p[1+nbin*RK::nper()]) : 1.0;
   // Extract the per-bin amplitudes ALWAYS -- this used to live inside the `if(talk)` block, so
   // calling with talk = false returned an empty `out` and any figure built from it came out blank.
   for(int b=0;b<nbin;b++){ RK::Bin& d=RK::B[b]; if(!d.h) continue;
     const double* q=p+1+b*RK::nper();
     RK::LASTBG[b]=std::max(0.0,q[1]);
     for(int i=0;i<RK::NR;i++) RK::LASTA[b][i]=0;
     // ⚠ with TIE12 on, q[2] is the DOUBLET amplitude, so the free resonances start at q[3].
     // Starting jj at 0 here read every amplitude one slot low.
     if(RK::PAIRS){ const double* rb=p+1+nbin*RK::nper();
       for(int i=0;i<RK::NR;i++){
         double w=(RK::RIDX[i]>=0)?std::max(0.0,rb[RK::RIDX[i]]):RK::WGT[i];
         RK::LASTA[b][i] = (RK::GRP[i]>=0 && d.in[i]) ? std::max(0.0,q[2+RK::GRP[i]])*w : 0.0; } }
     else { int jj=(RK::TIE12?1:0);
       if(RK::TIE12){ RK::LASTA[b][0]=d.in[0]?std::max(0.0,q[2]):0.0; }
       for(int i=(RK::TIE12?2:0);i<RK::NR;i++){ if(!RK::USE[i]) continue;
         RK::LASTA[b][i]=d.in[i]?std::max(0.0,q[2+jj]):0.0; jj++; } }
     if(RK::PAIRS){
       std::vector<double> row;
       for(int g=0;g<RK::NGRP;g++){
         double tot=0;
         for(int i=0;i<RK::NR;i++) if(RK::GRP[i]==g && d.in[i]){
           // the pair's combined height at the HIGHER member's peak
           double mu=(RK::EPP6[i]-RK::PPGS6)+d.gsMu, sum=0;
           const double* rb=p+1+nbin*RK::nper();
           for(int k=0;k<RK::NR;k++) if(RK::GRP[k]==g && d.in[k]){
             double w=(RK::RIDX[k]>=0)?std::max(0.0,rb[RK::RIDX[k]]):RK::WGT[k];
             sum += w*d.T[k].Eval(mu); }
           tot=std::max(tot,std::max(0.0,q[2+g])*sum); }
         row.push_back(tot); }
       if(out) out->push_back(row);
       continue; }
     std::vector<double> row; int j=0;
     if(RK::TIE12){ row.push_back(d.in[0]?q[2]*d.T[0].Eval((RK::EPP6[0]-RK::PPGS6)+d.gsMu):0); j=1; }
     for(int i=(RK::TIE12?2:0);i<RK::NR;i++){ if(!RK::USE[i]) continue;
       row.push_back(d.in[i]?q[2+j]*d.T[i].Eval((RK::EPP6[i]-RK::PPGS6)+d.gsMu):0); j++; }
     if(out) out->push_back(row); }
   if(talk){
     printf("\n================ %s   Ex [%.1f,%.1f], %d x %.0f deg bins, %d resonances ================\n",
            k.nm,RK::exLo,RK::exHi,nbin,TCMSTEP,RK::NUSE);
     printf("   shared shape par = %.4f     -2lnL = %.1f   (%d parameters)\n",p[0],m->MinValue(),npar);
     printf("   bin  theta_CM     N   chi2/ndf |");
     if(RK::TIE12) printf("   D12  ");
     for(int i=(RK::TIE12?2:0);i<RK::NR;i++) if(RK::USE[i]) printf("  %5.2f ",RK::EPP6[i]-RK::PPGS6);
     printf("\n");
     for(int b=0;b<nbin;b++){ RK::Bin& d=RK::B[b]; if(!d.h) continue;
       double c2=0;int nd=0;
       for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
         if(x<RK::exLo||x>RK::exHi||n<5) continue;
         double mm=RK::model(b,x,p); c2+=(n-mm)*(n-mm)/n; nd++; }
       printf("   %2d   %2.0f-%2.0f  %7.0f  %6.2f  |",b,d.cmLo,d.cmHi,d.N,c2/std::max(1,nd-RK::nper()));
       if(out) for(double x:(*out)[b]) printf(" %6.1f ",x);
       printf("\n"); }
   }
   double v=m->MinValue(); delete m; return v;
}
void resfitD(){ setCM(10,5,60); for(auto&k:CFG) rkRun(k,10,true); }
