/// @file prod_ex_pp.C
/// @brief Production 16C(p,p) result over many runs: combined Ex spectrum and elastic
/// KE-vs-theta from the genfit proton fits (run_XXXX_genfitter_pp.root), IC-gated 16C
/// beam. Nominal beam energy (vertex E-loss correction left for later calibration).
///
///   root -b -q 'pp/prod_ex_pp.C("run_0106,...,run_0115","/mnt/f/a1975/reco/")'

#include <tuple>
static double omega2(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
static double Ex2b(double m1,double m2,double m3,double m4,double Kp,double thl,double Ke){
   double Et1=Kp+m1,Et3=Ke+m3;double s=m1*m1+m2*m2+2*m2*Et1,uu=m2*m2+m3*m3-2*m2*Et3;
   double a=(cos(thl)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   return a<0?NAN:sqrt(a)-m4; }

void prod_ex_pp(TString runsCSV="run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,run_0114,run_0115",
                TString inDir="/mnt/f/a1975/reco/", TString suffix="_genfitter_pp", double Ebeam=192.0,
                double icMin=950, double icMax=1350, int icTbLo=1000, int icTbHi=1350, double chi2Cut=5.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   const double u=931.49401, m_C16=16.0147*u, m_p=1.007825*u;

   TH1F *hEx=new TH1F("hEx","16C(p,p') E_{x} (genfit, all runs);E_{x} [MeV];protons",200,-5,20);
   TH2F *hKT=new TH2F("hKT","elastic-gated proton KE vs #theta;#theta_{lab} [deg];KE [MeV]",90,0,95,110,0,45);
   TH1F *hVz=new TH1F("hVz","vertex z;z [mm];protons",110,-50,1000);

   TObjArray *runs=runsCSV.Tokenize(",");
   long nIC=0,nP=0; int nRun=0;
   for(int ri=0;ri<runs->GetEntries();++ri){
      TString run=((TObjString*)runs->At(ri))->GetString();
      TString gf=inDir+run+suffix+".root", ff=inDir+run+"_FRIB.root";
      if(gSystem->AccessPathName(gf)||gSystem->AccessPathName(ff)){ printf("skip %s\n",run.Data()); continue; }
      TFile *fg=TFile::Open(gf); TTree *tg=(TTree*)fg->Get("cbmsim");
      TFile *fc=TFile::Open(ff); TTree *tc=(TTree*)fc->Get("cbmsim");
      TClonesArray *te=nullptr,*re=nullptr;
      tg->SetBranchAddress("AtTrackingEvent",&te); tc->SetBranchAddress("AtRawEvent",&re);
      Long64_t N=std::min(tg->GetEntries(),tc->GetEntries());
      for(Long64_t i=0;i<N;++i){
         tc->GetEntry(i); double ic=-1;
         if(re->GetEntries()>0){ auto*raw=(AtRawEvent*)re->At(0);
            if(raw&&!raw->GetGenTraces().empty()){ auto&adc=raw->GetGenTraces()[0]->GetADC();
               double mx=-1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();++b) mx=std::max(mx,adc[b]); ic=mx; } }
         if(ic<icMin||ic>icMax) continue; nIC++;
         tg->GetEntry(i); if(te->GetEntries()==0) continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
         for(auto&ft:ev->GetFittedTracks()){ // already proton-gated + theta-windowed by the fitter
            if(!ft) continue; auto&k=ft->GetKinematics(); auto&m=ft->GetTrackMetadata(); if(!m) continue;
            double ndf=m->GetNdf(), c2=ndf>0?m->GetChi2()/ndf:1e9, ke=k.kineticEnergy;
            if(ke<=0||ke>1000||c2>chi2Cut) continue;
            double ex=Ex2b(m_C16,m_p,m_p,m_C16,Ebeam,k.theta,ke); if(std::isnan(ex)) continue;
            hEx->Fill(ex); hKT->Fill(k.theta*TMath::RadToDeg(),ke); hVz->Fill(ft->GetVertex(0).Z()); nP++;
         }
      }
      fg->Close(); fc->Close(); nRun++;
      printf("processed %s (running total protons=%ld)\n",run.Data(),nP);
   }
   printf("\nruns=%d  IC-gated events=%ld  protons->Ex=%ld\n",nRun,nIC,nP);
   hEx->Fit("gaus","Q0","",-2,2);
   if(hEx->GetFunction("gaus"))
      printf("elastic peak: mean=%.3f MeV  sigma=%.3f MeV\n",
             hEx->GetFunction("gaus")->GetParameter(1), hEx->GetFunction("gaus")->GetParameter(2));

   TCanvas *c=new TCanvas("c","prod pp",1500,520); c->Divide(3,1);
   c->cd(1); hEx->Draw("hist"); if(hEx->GetFunction("gaus")){hEx->GetFunction("gaus")->SetNpx(500);hEx->GetFunction("gaus")->Draw("same");}
   c->cd(2); hKT->Draw("colz");
   c->cd(3); hVz->Draw();
   c->SaveAs("/tmp/prod_ex_pp.png");
   printf("wrote /tmp/prod_ex_pp.png\n");
}
