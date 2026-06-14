/// @file excorr_pd.C
/// @brief Flatten the Ex-vs-theta_cm dependence of the genfit 16C(p,d)15C data.
/// A discrete state must have Ex independent of theta_cm, but the reconstructed g.s.
/// band drifts with theta_cm (a kinematic systematic that smears the integrated
/// spectrum). This fits the g.s.-region Ex(theta_cm) profile with a quadratic and
/// subtracts it event-by-event, bringing every angle to the same Ex. Shows the
/// g.s.-centroid drift + the fit, and the before/after integrated spectrum with the
/// resolution gain. Ported from pd/excorr_pd_a1975.C (UKF) to the genfit cache.
/// Prints the pol2 coefficients so plot_pd_angbins.C can apply the same correction.
///
///   root -l -b -q 'pp/excorr_pd.C'

static double e_om2(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
static std::pair<double,double> e_kine(double Kp,double thl,double Ke){
   const double u=931.49401; double m1=16.0147*u,m2=1.007825*u,m3=2.014102*u,m4=15.010599*u;
   double Et1=Kp+m1,Et3=Ke+m3,Et4=Et1+m2-Et3; double s=m1*m1+m2*m2+2*m2*Et1,uu=m2*m2+m3*m3-2*m2*Et3;
   double a=(cos(thl)*e_om2(s,m1*m1,m2*m2)*e_om2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   if(a<0) return {NAN,NAN}; double m4x=std::sqrt(a); double ex=m4x-m4;
   double t=m2*m2+m4x*m4x-2*m2*Et4;
   double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4x*m4x)+(m1*m1-m2*m2)*(m3*m3-m4x*m4x))
                                    /(e_om2(s,m1*m1,m2*m2)*e_om2(s,m3*m3,m4x*m4x)));
   return {ex, tcm*TMath::RadToDeg()};
}
// g.s. centroid + FWHM from a doublet fit with the 0.74 separation fixed (stable)
static void gsFit(TH1 *h,double &gs,double &fwhm){
   double pk=h->GetBinCenter(h->GetMaximumBin());
   TF1 dg("dg","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1]-0.740)/[4])^2)",pk-1.1,pk+1.1);
   double mx=h->GetMaximum(); dg.SetParameters(0.6*mx,pk-0.4,0.30,mx,0.30);
   dg.SetParLimits(1,pk-1.1,pk+0.3); dg.SetParLimits(2,0.10,0.8); dg.SetParLimits(4,0.10,0.8);
   h->Fit(&dg,"QRN"); gs=dg.GetParameter(1); fwhm=2.3548*dg.GetParameter(2);
}

void excorr_pd(double exShift=-0.34, double Ebeam=192, double chi2max=5, double icMin=950, double icMax=1350,
               double tcLo=15, double tcHi=120, int nProf=21, TString cache="/tmp/pd_kin.root",
               TString out="/tmp/pd_excorr.png")
{
   gStyle->SetOptStat(0);
   TFile *f=TFile::Open(cache); if(!f||f->IsZombie()){ printf("cannot open %s\n",cache.Data()); return; }
   TTree *t=(TTree*)f->Get("pk"); if(!t){ printf("no tree pk\n"); return; }
   float ke,theta,vz,chi2ndf,ic; t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta);
   t->SetBranchAddress("vz",&vz); t->SetBranchAddress("chi2ndf",&chi2ndf); t->SetBranchAddress("ic",&ic);
   std::vector<float> vex,vtc;
   for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
      if(chi2ndf>chi2max) continue; if(ic<icMin||ic>icMax) continue;
      auto[ex,tcm]=e_kine(Ebeam,theta*TMath::DegToRad(),ke); if(std::isnan(ex)) continue;
      vex.push_back(ex+exShift); vtc.push_back(tcm); }
   f->Close();
   printf("excorr_pd: %zu selected deuterons\n",vex.size());

   // g.s. centroid per theta_cm bin -> measured systematic drift (the genfit drift is
   // inverted-U with a steep back-angle drop, too sharp for a quadratic, so we subtract
   // the MEASURED drift interpolated through the points instead of a pol fit).
   TGraph *gpk=new TGraph(); double bw=(tcHi-tcLo)/nProf;
   for(int b=0;b<nProf;++b){ double t0=tcLo+b*bw,t1=t0+bw,tc=0.5*(t0+t1);
      TH1F hb("hb","",120,-3,5); hb.SetDirectory(nullptr); long n=0;
      for(size_t i=0;i<vex.size();++i) if(vtc[i]>=t0&&vtc[i]<t1){ hb.Fill(vex[i]); ++n; }
      if(n>200){ double gs,fw; gsFit(&hb,gs,fw); gpk->SetPoint(gpk->GetN(),tc,gs); } }
   double xfirst,xlast,yd; gpk->GetPoint(0,xfirst,yd); gpk->GetPoint(gpk->GetN()-1,xlast,yd);
   auto drift=[&](double tc){ return gpk->Eval(std::min(xlast,std::max(xfirst,tc))); }; // measured drift, clamped
   TSpline3 *spl=new TSpline3("spl",gpk); // smooth version for the overlay
   printf("measured g.s. drift: %d points over theta_cm %.0f-%.0f\n",gpk->GetN(),xfirst,xlast);

   // before / after spectra (subtract the measured per-angle g.s. position)
   TH1F *hb=new TH1F("hbef","",200,-4,12); hb->SetDirectory(nullptr);
   TH1F *ha=new TH1F("haft","",200,-4,12); ha->SetDirectory(nullptr);
   for(size_t i=0;i<vex.size();++i){ hb->Fill(vex[i]); ha->Fill(vex[i]-drift(vtc[i])); }
   double gb,fbw,ga,faw; gsFit(hb,gb,fbw); gsFit(ha,ga,faw);
   printf("\n            g.s.pos   g.s.FWHM\n  BEFORE    %+6.3f    %6.3f\n  AFTER     %+6.3f    %6.3f\n  => FWHM %.0f%% better\n",
          gb,fbw,ga,faw,100*(fbw-faw)/fbw);

   auto*c=new TCanvas("c_ex","excorr",1400,560); c->Divide(2,1);
   c->cd(1); gpk->SetTitle("g.s. peak position vs #theta_{cm} (subtracted = red);#theta_{cm} [deg];E_{x}^{g.s.} [MeV]");
   gpk->SetMarkerStyle(20); gpk->SetMarkerColor(kBlue+1); gpk->Draw("AP"); spl->SetLineColor(kRed+1); spl->Draw("same");
   c->cd(2); hb->SetLineColor(kGray+2); hb->SetLineWidth(2); ha->SetLineColor(kBlue+1); ha->SetLineWidth(2);
   hb->SetTitle("^{15}C E_{x}: before (grey) vs #theta_{cm}-corrected (blue);E_{x}(^{15}C) [MeV];counts");
   hb->Draw("hist"); ha->Draw("hist same");
   auto*tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.035);
   tx->DrawLatex(0.5,0.84,Form("g.s. FWHM %.3f #rightarrow %.3f MeV",fbw,faw));
   c->SaveAs(out);
   printf("wrote %s\n",out.Data());
}
