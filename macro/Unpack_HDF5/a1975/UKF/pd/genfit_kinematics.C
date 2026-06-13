/// @file genfit_kinematics.C
/// @brief GenFit (p,d) kinematics + EFFICIENCY cut-flow over the 84-run set.
/// Quantifies where deuterons are lost (fit -> KE -> chi2 -> PRA match -> spy -> PID
/// gate), plots KE-vs-theta_lab with theory curves, and the PID plane vs the gate.
/// Caches a rich ntuple (genfit_kin_rich.root) for instant re-analysis.
///
///   root -b -q 'pd/genfit_kinematics.C'

#include <map>
static double omega2(double x,double y,double z){ return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
static std::pair<double,double> kine2b(double m1,double m2,double m3,double m4,double Kp,double thLab,double Ke){
   double Et1=Kp+m1,Et3=Ke+m3,Et4=Et1+m2-Et3;
   double s=m1*m1+m2*m2+2*m2*Et1, uu=m2*m2+m3*m3-2*m2*Et3;
   double m4e=std::sqrt((std::cos(thLab)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2);
   double t=m2*m2+m4e*m4e-2*m2*Et4;
   double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4e*m4e)+(m1*m1-m2*m2)*(m3*m3-m4e*m4e))/(omega2(s,m1*m1,m2*m2)*omega2(s,m3*m3,m4e*m4e)));
   return {m4e-m4, tcm*TMath::RadToDeg()};
}

void genfit_kinematics(double Ebeam=192.0, double chi2Cut=5.0)
{
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtTools.so");
   const double u=931.49401, m1=16.0147*u, m2=1.007825*u, m3=2.01410178*u, m4_0=15.0105993*u;
   auto pid=AtTools::AtParticleID::LoadJSON("pid/deuteron_band.json");
   AtTools::AtSpyralPID spy; spy.SetBField(2.85);

   // cut-flow counters
   long nFit=0,nKE=0,nChi2=0,nMatch=0,nSpy=0,nGate=0;
   TH2F *hk=new TH2F("hk","^{16}C(p,d)^{15}C GenFit kinematics + theory;#theta_{lab} [deg];KE_{d} [MeV]",250,0,60,250,0,50);
   hk->SetDirectory(nullptr);
   TH2F *hpid=new TH2F("hpid","GenFit PID plane (all valid);#sqrt{dEdx};B#rho [Tm]",300,0,30,300,0,1.2); hpid->SetDirectory(nullptr);
   TH2F *hpidg=new TH2F("hpidg","",300,0,30,300,0,1.2); hpidg->SetDirectory(nullptr);
   TFile *fc=new TFile("genfit_kin_rich.root","RECREATE");
   TNtuple *nt=new TNtuple("gk","genfit pd rich","ke:thlab:thcm:ex:brho:sqrtdedx:chi2ndf:gated");

   void *dir=gSystem->OpenDirectory("/mnt/f/a1975/reco_gf"); std::vector<TString> runs;
   const char*ent; while((ent=gSystem->GetDirEntry(dir))){ TString s(ent); if(s.BeginsWith("run_")&&s.EndsWith("_genfit.root")) runs.push_back(s(0,8)); }
   gSystem->FreeDirectory(dir); std::sort(runs.begin(),runs.end());

   for(auto&run:runs){
      TString gf=TString("/mnt/f/a1975/reco_gf/")+run+"_genfit.root", rc=TString("/mnt/f/a1975/reco/")+run+"_reco.root";
      if(gSystem->AccessPathName(gf)||gSystem->AccessPathName(rc))continue;
      TFile*fg=TFile::Open(gf); TTree*tg=(TTree*)fg->Get("cbmsim");
      TFile*fr=TFile::Open(rc); TTree*tr=(TTree*)fr->Get("cbmsim");
      if(!tg||!tr){ if(fg)fg->Close(); if(fr)fr->Close(); continue; }
      TClonesArray*te=nullptr,*pe=nullptr; tg->SetBranchAddress("AtTrackingEvent",&te); tr->SetBranchAddress("AtPatternEvent",&pe);
      Long64_t N=std::min(tg->GetEntries(),tr->GetEntries());
      for(Long64_t i=0;i<N;++i){
         tg->GetEntry(i); if(te->GetEntries()==0)continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev)continue;
         tr->GetEntry(i); if(pe->GetEntries()==0)continue; auto*pat=(AtPatternEvent*)pe->At(0); if(!pat)continue;
         std::vector<AtTrack>&cand=pat->GetTrackCand(); std::map<int,AtTrack*> byID; for(auto&t2:cand)byID[t2.GetTrackID()]=&t2;
         for(auto&ft:ev->GetFittedTracks()){
            if(!ft)continue; ++nFit; auto&k=ft->GetKinematics();
            double ndf=ft->GetTrackMetadata()->GetNdf(),chi2=ft->GetTrackMetadata()->GetChi2(),c2n=ndf>0?chi2/ndf:1e9;
            double ke=k.kineticEnergy, thd=k.theta; // genfit theta already in DEGREES
            if(ke<=0||ke>1000)continue; ++nKE;
            if(c2n>chi2Cut)continue; ++nChi2;
            auto it=byID.find(ft->GetTrackID()); if(it==byID.end())continue; ++nMatch;
            auto r=spy.Estimate(*it->second); if(!r.valid)continue; ++nSpy;
            bool g=pid.IsInside(r.sqrtdEdx,r.brho); if(g)++nGate;
            hpid->Fill(r.sqrtdEdx,r.brho); if(g)hpidg->Fill(r.sqrtdEdx,r.brho);
            double thR=thd*TMath::DegToRad(); auto[ex,thcm]=kine2b(m1,m2,m3,m4_0,Ebeam,thR,ke);
            if(g){ hk->Fill(thd,ke); }
            nt->Fill(ke,thd,thcm,ex,r.brho,r.sqrtdEdx,c2n,g?1:0);
         }
      }
      fg->Close(); fr->Close();
   }
   nt->Write(); fc->Close();

   printf("\n=== GenFit (p,d) cut-flow efficiency (84 runs) ===\n");
   printf("  fitted tracks      : %ld\n",nFit);
   printf("  KE in (0,1000)     : %ld  (%.1f%%)\n",nKE,100.0*nKE/nFit);
   printf("  chi2/ndf < %.0f      : %ld  (%.1f%% of fit, %.1f%% of KE-ok)\n",chi2Cut,nChi2,100.0*nChi2/nFit,100.0*nChi2/nKE);
   printf("  matched to PRA     : %ld  (%.1f%%)\n",nMatch,100.0*nMatch/nChi2);
   printf("  spy valid          : %ld  (%.1f%%)\n",nSpy,100.0*nSpy/nMatch);
   printf("  deuteron PID gated : %ld  (%.1f%% of spy-valid)\n",nGate,100.0*nGate/nSpy);

   // KE vs theta_lab + theory curves
   auto curve=[&](double Ex,int col)->TGraph*{
      double E1=Ebeam+m1,p1=std::sqrt(E1*E1-m1*m1),Etot=E1+m2,m4=m4_0+Ex; std::vector<double>ut,uk,lt,lk;
      for(double thd=0;thd<=55;thd+=0.1){ double c=std::cos(thd*TMath::DegToRad());
         double A=Etot*Etot-m4*m4-p1*p1+m3*m3, a=4*Etot*Etot-4*p1*p1*c*c, b=-4*Etot*A, cc=A*A+4*p1*p1*c*c*m3*m3, disc=b*b-4*a*cc;
         if(disc<0||a==0)continue; for(int s=-1;s<=1;s+=2){ double E3=(-b+s*std::sqrt(disc))/(2*a); if(E3<=m3||(2*Etot*E3-A)*c<=0)continue;
            if(s>0){ut.push_back(thd);uk.push_back(E3-m3);}else{lt.push_back(thd);lk.push_back(E3-m3);} } }
      TGraph*g=new TGraph(); for(size_t i=0;i<ut.size();++i)g->SetPoint(g->GetN(),ut[i],uk[i]); for(size_t i=lt.size();i-->0;)g->SetPoint(g->GetN(),lt[i],lk[i]);
      g->SetLineColor(col); g->SetLineWidth(2); return g; };
   TCanvas*c1=new TCanvas("c1","kin",950,700); gPad->SetLogz(); hk->Draw("colz");
   struct{double ex;int col;const char*l;}st[]={{0,kRed,"g.s."},{0.74,kMagenta+2,"0.74"},{3.10,kGreen+2,"3.10"},{4.66,kOrange+7,"4.66"}};
   TLegend*lg=new TLegend(0.6,0.66,0.88,0.88); lg->SetHeader("^{15}C E_{x}");
   for(auto&s:st){TGraph*g=curve(s.ex,s.col); g->Draw("L same"); lg->AddEntry(g,s.l,"l");} lg->Draw();
   c1->SaveAs("pd/plots/genfit_kinematics.png");

   // PID plane: all valid (grey) + gated (color)
   TCanvas*c2=new TCanvas("c2","pid",950,700); gPad->SetLogz();
   hpid->SetTitle(Form("GenFit (p,d) PID: %ld valid, %ld gated (%.0f%%);#sqrt{dEdx};B#rho [Tm]",nSpy,nGate,100.0*nGate/nSpy));
   hpid->Draw("colz"); hpidg->SetMarkerColor(kRed); hpidg->SetMarkerStyle(1); hpidg->Draw("same");
   c2->SaveAs("pd/plots/genfit_pidplane.png");
   printf("saved pd/plots/genfit_kinematics.png + genfit_pidplane.png (cache genfit_kin_rich.root)\n");
}
