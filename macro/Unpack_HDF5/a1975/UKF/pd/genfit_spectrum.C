/// @file genfit_spectrum.C
/// @brief Final full-statistics 15C(p,d) Ex spectrum from the 84-run GenFit fits,
/// with the theta_cm correction re-derived ON THE GENFIT DATA (genfit's systematic
/// differs from the UKF's). Reads reco_gf/ (genfit kinematics) + reco/ (PRA tracks
/// for the deuteron PID gate + theta for kinematics), event-matched.
///
///   root -b -q 'pd/genfit_spectrum.C'

#include <map>
static double omega2(double x, double y, double z)
{ return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
// returns {Ex, theta_cm[deg]} for inverse kinematics (recoil angle)
static std::pair<double,double> kine2b(double m1,double m2,double m3,double m4,double Kp,double thLab,double Ke)
{
   double Et1=Kp+m1, Et3=Ke+m3, Et4=Et1+m2-Et3;
   double s=m1*m1+m2*m2+2*m2*Et1, uu=m2*m2+m3*m3-2*m2*Et3;
   double m4e=std::sqrt((std::cos(thLab)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2);
   double ex=m4e-m4;
   double t=m2*m2+m4e*m4e-2*m2*Et4;
   double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4e*m4e)+(m1*m1-m2*m2)*(m3*m3-m4e*m4e))/(omega2(s,m1*m1,m2*m2)*omega2(s,m3*m3,m4e*m4e)));
   return {ex, tcm*TMath::RadToDeg()};
}
static double gsPos(TH1F *h)
{
   double mx=h->GetMaximum();
   TF1 dg("dgp","[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x",-1.5,2.3);
   dg.SetParameters(0.45,0.30,0.6*mx,0.6*mx,0.1*mx,0.0);
   dg.SetParLimits(0,-0.4,1.0); dg.SetParLimits(1,0.12,0.55);
   dg.SetParLimits(2,0,2*mx); dg.SetParLimits(3,0,2*mx);
   h->Fit(&dg,"QRN"); return dg.GetParameter(0);
}

void genfit_spectrum(double Ebeam=192.0, double chi2Cut=1e9, double exShift=0.0)
{
   gStyle->SetOptStat(0);
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtTools.so");
   const double u=931.49401, mC16=16.0147*u, mp=1.007825*u, md=2.01410178*u, mC15=15.0105993*u;
   auto pid=AtTools::AtParticleID::LoadJSON("pid/deuteron_band.json");
   AtTools::AtSpyralPID spy; spy.SetBField(2.85);

   std::vector<float> ex0, tc;
   long nD=0;
   void *dir=gSystem->OpenDirectory("/mnt/f/a1975/reco_gf");
   std::vector<TString> runs;
   const char*ent; while((ent=gSystem->GetDirEntry(dir))){ TString s(ent); if(s.BeginsWith("run_")&&s.EndsWith("_genfit.root")) runs.push_back(s(0,8)); }
   gSystem->FreeDirectory(dir);
   std::sort(runs.begin(),runs.end());
   printf("genfit runs found: %zu\n",(size_t)runs.size());

   for(auto&run:runs){
      TString gf=TString("/mnt/f/a1975/reco_gf/")+run+"_genfit.root";
      TString pf=TString("/mnt/f/a1975/reco_gf/")+run+"_pid.root"; // cached PID, no spy recompute
      if(gSystem->AccessPathName(gf)||gSystem->AccessPathName(pf)) continue;
      TFile*fg=TFile::Open(gf); TTree*tg=(TTree*)fg->Get("cbmsim");
      TFile*fr=TFile::Open(pf); TTree*tr=(TTree*)fr->Get("cbmsim");
      if(!tg||!tr){ if(fg)fg->Close(); if(fr)fr->Close(); continue; }
      TClonesArray*te=nullptr,*pe=nullptr;
      tg->SetBranchAddress("AtTrackingEvent",&te); tr->SetBranchAddress("AtPIDEvent",&pe);
      Long64_t N=std::min(tg->GetEntries(),tr->GetEntries());
      for(Long64_t i=0;i<N;++i){
         tg->GetEntry(i); if(te->GetEntries()==0)continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev)continue;
         tr->GetEntry(i); if(pe->GetEntries()==0)continue; auto*pidev=(AtPIDEvent*)pe->At(0); if(!pidev)continue;
         std::map<int,const AtTools::AtSpyralResult*> byID; for(auto&sr:pidev->GetSpyral())byID[sr.trackID]=&sr;
         for(auto&ft:ev->GetFittedTracks()){
            if(!ft)continue; auto&k=ft->GetKinematics();
            double ndf=ft->GetTrackMetadata()->GetNdf(),chi2=ft->GetTrackMetadata()->GetChi2(),c2n=ndf>0?chi2/ndf:1e9;
            double ke=k.kineticEnergy, thR=k.theta*TMath::DegToRad(); // genfit theta in DEGREES
            if(ke<=0||c2n>chi2Cut)continue; // physical KE + PID gate only (chi2 off by default)
            auto it=byID.find(ft->GetTrackID()); if(it==byID.end())continue;
            const auto&r=*it->second; if(!r.valid||!pid.IsInside(r.sqrtdEdx,r.brho))continue;
            auto [ex,thcm]=kine2b(mC16,mp,md,mC15,Ebeam,thR,ke);
            if(std::isnan(ex)||std::isnan(thcm))continue;
            if(thcm<10||thcm>160)continue;
            ex0.push_back(ex+exShift); tc.push_back(thcm); ++nD;
         }
      }
      fg->Close(); fr->Close();
   }
   printf("deuterons in spectrum: %ld\n",nD);

   // --- raw self-calibrated doublet (before correction) ---
   auto dres=[&](std::vector<float>&v,double&fw,double&pv,double&gs){
      TH1F*h=new TH1F("ht","",240,-6,14); h->SetDirectory(nullptr); for(float e:v)h->Fill(e);
      gs=gsPos(h); double mx=h->GetMaximum();
      TF1 dg("dg2","[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x",gs-1.5,gs+1.6);
      dg.SetParameters(gs,0.28,0.6*mx,0.6*mx,0.1*mx,0.0); dg.SetParLimits(0,gs-0.4,gs+0.4); dg.SetParLimits(1,0.12,0.55);
      h->Fit(&dg,"QRN"); fw=2.3548*dg.GetParameter(1);
      TH1F*hs=new TH1F("hs2","",240,-6-gs,14-gs); hs->SetDirectory(nullptr); for(float e:v)hs->Fill(e-gs); hs->Smooth(1);
      auto mm=[&](double lo,double hi,bool mn){double m=mn?1e18:0; for(int b=hs->FindBin(lo);b<=hs->FindBin(hi);++b)m=mn?std::min(m,hs->GetBinContent(b)):std::max(m,hs->GetBinContent(b)); return m;};
      pv=std::min(mm(-0.25,0.25,false),mm(0.5,1.0,false))/std::max(1.0,mm(0.25,0.55,true)); delete h; delete hs;
   };
   double f0,p0,g0; dres(ex0,f0,p0,g0);

   // --- theta_cm correction: per-bin g.s. position (re-derived on genfit) -> pol2 ---
   TF1 sysT("sysT","pol2",20,150);
   TGraph*gT=new TGraph(); int nb=13; double vlo=20,vhi=150,bw=(vhi-vlo)/nb;
   for(int b=0;b<nb;++b){ double c0=vlo+b*bw,c1=c0+bw; TH1F*hb=new TH1F("hbb","",120,-3,5); hb->SetDirectory(nullptr); long n=0;
      for(size_t i=0;i<ex0.size();++i) if(tc[i]>=c0&&tc[i]<c1){hb->Fill(ex0[i]-g0);++n;}
      if(n>250) gT->SetPoint(gT->GetN(),0.5*(c0+c1),gsPos(hb)); delete hb; }
   gT->Fit(&sysT,"QRN");
   std::vector<float> exc(ex0.size());
   for(size_t i=0;i<ex0.size();++i){ double a=std::min(150.0,std::max(20.0,(double)tc[i])); exc[i]=ex0[i]-g0-sysT.Eval(a); }
   double f1,p1,g1; dres(exc,f1,p1,g1);

   printf("\n                  doublet_FWHM   peak/valley   nD\n");
   printf("  raw (calib)       %6.3f        %.3f        %ld\n",f0,p0,nD);
   printf("  + theta_cm corr   %6.3f        %.3f\n",f1,p1);
   // NOTE: for GENFIT the theta_cm correction HURTS (the UKF needed it to fix a
   // theta_cm-dependent fit bias; genfit's fit is unbiased so the correction only
   // adds spread). The RAW self-calibrated spectrum is the best — plot that.
   printf("  => theta_cm correction %s for genfit (raw is best)\n", f1<f0?"helps":"HURTS");

   // cache kinematics for instant re-plots
   TFile *fc=new TFile("genfit_kin.root","RECREATE");
   TNtuple *nt=new TNtuple("gk","genfit (p,d) kin","ex:thcm");
   for(size_t i=0;i<ex0.size();++i) nt->Fill(ex0[i],tc[i]); nt->Write(); fc->Close();

   // --- draw RAW self-calibrated spectrum (the best genfit result) ---
   TH1F*hf=new TH1F("hf","",200,-3,10); hf->SetDirectory(nullptr);
   for(size_t i=0;i<ex0.size();++i) hf->Fill(ex0[i]-g0);
   TCanvas*c=new TCanvas("c","genfit_spectrum",1000,650);
   hf->SetLineColor(kBlue+1); hf->SetLineWidth(2); hf->SetFillColorAlpha(kAzure+1,0.15);
   hf->SetTitle(Form("^{15}C(p,d) GenFit full-stats, 84 runs (n=%ld);E_{x}(^{15}C) [MeV];counts",nD));
   hf->Draw("hist");
   TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035);
   tl.DrawLatex(0.50,0.83,Form("doublet FWHM = %.3f MeV  (no #theta_{cm} corr needed)",f0));
   tl.DrawLatex(0.50,0.78,Form("peak/valley = %.2f",p0));
   tl.DrawLatex(0.50,0.73,"g.s. self-calibrated (no energy offset)");
   c->SaveAs("pd/plots/genfit_spectrum_full.png");
   printf("saved pd/plots/genfit_spectrum_full.png (RAW, %.3f MeV)\n",f0);
}
