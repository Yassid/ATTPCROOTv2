/// @file ex_genfit_vs_ukf.C
/// @brief 16C excitation-energy spectrum from the genfit protons vs the UKF protons,
/// on the SAME events with IDENTICAL selection (IC gate + proton PID gate + chi2 +
/// theta window) and the SAME two-body kinematics as pp/ex_a1975.C. Isolates the
/// fitter (genfit vs UKF); elastic -> Ex~0, 16C excited states appear as peaks.
///
///   root -b -q 'pp/ex_genfit_vs_ukf.C("run_0106","/tmp/run_0106_genfitter_pp.root")'

#include <map>
#include <tuple>

static double omega2(double x, double y, double z) { return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }

static std::tuple<double,double> kine_2b(double m1,double m2,double m3,double m4,double K_proj,double thetalab,double K_eject)
{
   double Et1=K_proj+m1, Et2=m2, Et3=K_eject+m3, Et4=Et1+Et2-Et3;
   double s=m1*m1+m2*m2+2*m2*Et1, uu=m2*m2+m3*m3-2*m2*Et3;
   double m4_ex=std::sqrt((std::cos(thetalab)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)
                          -(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2);
   double Ex=m4_ex-m4;
   double t=m2*m2+m4_ex*m4_ex-2*m2*Et4;
   double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4_ex*m4_ex)+(m1*m1-m2*m2)*(m3*m3-m4_ex*m4_ex))
                                    /(omega2(s,m1*m1,m2*m2)*omega2(s,m3*m3,m4_ex*m4_ex)));
   return {Ex, tcm*TMath::RadToDeg()};
}

void ex_genfit_vs_ukf(TString run="run_0106", TString genfitFile="/tmp/run_0106_genfitter_pp.root",
                      TString inDir="/mnt/f/a1975/reco/", TString gateFile="pid/proton_band.json",
                      double Ebeam=192.0, double icMin=950, double icMax=1350, int icTbLo=1000, int icTbHi=1350,
                      double chi2Cut=5.0, double bField=2.85, double thMin=10, double thMax=170)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   const double u=931.49401, m_C16=16.0147*u, m_p=1.007825*u;
   auto pid=AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy; spy.SetBField(bField);

   TFile *fg=TFile::Open(genfitFile);                 TTree *tg=(TTree*)fg->Get("cbmsim");
   TFile *fu=TFile::Open(inDir+run+"_ukf.root");       TTree *tu=(TTree*)fu->Get("cbmsim");
   TFile *fc=TFile::Open(inDir+run+"_FRIB.root");      TTree *tc=(TTree*)fc->Get("cbmsim");
   if(!tg||!tu||!tc){printf("missing input\n");return;}
   TClonesArray *teg=nullptr,*teu=nullptr,*re=nullptr;
   tg->SetBranchAddress("AtTrackingEvent",&teg);
   tu->SetBranchAddress("AtTrackingEvent",&teu);
   tc->SetBranchAddress("AtRawEvent",&re);

   TH1F *hG=new TH1F("hG","16C(p,p') E_{x}: genfit vs UKF;E_{x} [MeV];protons",120,-5,25);
   TH1F *hU=new TH1F("hU","UKF",120,-5,25);
   hG->SetLineColor(kRed+1); hG->SetLineWidth(2);
   hU->SetLineColor(kAzure+2); hU->SetLineWidth(2);

   // process one fit-event into a histogram (shared selection)
   auto fill=[&](TClonesArray*te, TH1F*h)->int{
      if(te->GetEntries()==0) return 0;
      auto *ev=(AtTrackingEvent*)te->At(0); if(!ev) return 0;
      std::vector<AtTrack> orig=ev->GetTrackArray();
      std::map<int,AtTrack*> byID; for(auto&tr:orig) byID[tr.GetTrackID()]=&tr;
      int nf=0;
      for(auto&ft:ev->GetFittedTracks()){
         if(!ft) continue;
         auto&k=ft->GetKinematics();
         auto&m=ft->GetTrackMetadata(); if(!m) continue;
         double ndf=m->GetNdf(), c2n=ndf>0?m->GetChi2()/ndf:1e9;
         double ke=k.kineticEnergy, thDeg=k.theta*TMath::RadToDeg();
         if(ke<=0||ke>1000||c2n>chi2Cut||thDeg<thMin||thDeg>thMax) continue;
         auto it=byID.find(ft->GetTrackID()); if(it==byID.end()) continue;
         auto r=spy.Estimate(*it->second); if(!r.valid||!pid.IsInside(r.sqrtdEdx,r.brho)) continue;
         auto [ex,thcm]=kine_2b(m_C16,m_p,m_p,m_C16,Ebeam,k.theta,ke);
         if(std::isnan(ex)) continue;
         h->Fill(ex); nf++;
      }
      return nf;
   };

   Long64_t N=std::min(tg->GetEntries(),std::min(tu->GetEntries(),tc->GetEntries()));
   long nIC=0,nG=0,nU=0;
   for(Long64_t i=0;i<N;++i){
      tc->GetEntry(i);
      double ic=-1;
      if(re->GetEntries()>0){ auto*raw=(AtRawEvent*)re->At(0);
         if(raw&&!raw->GetGenTraces().empty()){ auto&adc=raw->GetGenTraces()[0]->GetADC();
            double mx=-1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();++b) mx=std::max(mx,adc[b]); ic=mx; } }
      if(ic<icMin||ic>icMax) continue; // IC gate (16C beam)
      nIC++;
      tg->GetEntry(i); nG+=fill(teg,hG);
      tu->GetEntry(i); nU+=fill(teu,hU);
   }
   printf("events=%lld  IC-gated=%ld  genfit protons=%ld  UKF protons=%ld\n",N,nIC,nG,nU);

   TCanvas *c=new TCanvas("c","ex cmp",900,650);
   hU->Draw("hist"); hG->Draw("hist same");
   double ymax=std::max(hG->GetMaximum(),hU->GetMaximum()); hU->SetMaximum(1.2*ymax);
   TLegend *lg=new TLegend(0.6,0.75,0.88,0.88);
   lg->AddEntry(hG,Form("genfit (%ld)",nG),"l");
   lg->AddEntry(hU,Form("UKF (%ld)",nU),"l");
   lg->Draw();
   TLine *l0=new TLine(0,0,0,1.2*ymax); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
   c->SaveAs("/tmp/ex_genfit_vs_ukf.png");
   printf("genfit Ex: mean=%.2f rms=%.2f  |  UKF Ex: mean=%.2f rms=%.2f  | wrote /tmp/ex_genfit_vs_ukf.png\n",
          hG->GetMean(),hG->GetRMS(),hU->GetMean(),hU->GetRMS());
}
