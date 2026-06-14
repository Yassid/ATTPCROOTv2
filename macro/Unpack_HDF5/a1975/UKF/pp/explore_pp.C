/// @file explore_pp.C
/// @brief Interactive explorer for the 16C(p,p) genfit protons. Reads the cached
/// kinematics ntuple (/tmp/pp_kin.root from pp/cache_pp_run.C) into memory and lets
/// you change beam energy / cuts / binning with GUI controls and redraw instantly.
/// Shows: Ex spectrum (Gaussian-fit elastic peak), KE-vs-theta_lab, Ex-vs-theta_cm.
/// 2D TCanvas only (works over X11; no Eve/OpenGL).
///
/// Build the cache once (parallel), then run interactively (NOT -b):
///   ls .../run_*_reco.root | grep -oE 'run_[0-9]+'|sed 's/run_//'|xargs -P4 -I{} \
///     root -b -q 'pp/cache_pp_run.C("run_{}","/tmp/ppkin_run{}.root")'
///   hadd -f /tmp/pp_kin.root /tmp/ppkin_run*.root
///   root -l 'pp/explore_pp.C'
///
/// Controls: Ebeam, chi2/ndf max, theta min/max, IC min/max, Ex bins/lo/hi, [Redraw] [Save].

#include <vector>

static double omega2(double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
static std::pair<double,double> kine2b(double m1,double m2,double m3,double m4,double Kp,double thl,double Ke){
   double Et1=Kp+m1,Et3=Ke+m3,Et4=Et1+m2-Et3; double s=m1*m1+m2*m2+2*m2*Et1,uu=m2*m2+m3*m3-2*m2*Et3;
   double a=(cos(thl)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   if(a<0) return {NAN,NAN}; double m4x=std::sqrt(a); double ex=m4x-m4;
   double t=m2*m2+m4x*m4x-2*m2*Et4;
   double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4x*m4x)+(m1*m1-m2*m2)*(m3*m3-m4x*m4x))
                                    /(omega2(s,m1*m1,m2*m2)*omega2(s,m3*m3,m4x*m4x)));
   return {ex, tcm*TMath::RadToDeg()};
}

// proton KE vs lab angle locus for p(16C,16C*)p at excitation Ex (relativistic 2-body,
// CM-angle scan): beam 16C (m1,KE=Eb), target p at rest, ejectile p (m3), residual 16C* (m4=mC+Ex).
static TGraph *kinLine(double Eb,double Ex,double m1,double m2,double m3,double m4_0,Color_t col,int style){
   double m4=m4_0+Ex;
   double E1=Eb+m1, P=std::sqrt(E1*E1-m1*m1), W=E1+m2, s=W*W-P*P;
   auto*g=new TGraph(); g->SetLineColor(col); g->SetLineWidth(2); g->SetLineStyle(style);
   if(s<=(m3+m4)*(m3+m4)) return g; // below threshold
   double rs=std::sqrt(s), E3cm=(s+m3*m3-m4*m4)/(2*rs), p3cm=std::sqrt(std::max(0.,E3cm*E3cm-m3*m3));
   double beta=P/W, gamma=W/rs;
   for(double tc=0;tc<=180;tc+=0.5){ double c=std::cos(tc*TMath::DegToRad()),sn=std::sin(tc*TMath::DegToRad());
      double Elab=gamma*(E3cm+beta*p3cm*c), pz=gamma*(p3cm*c+beta*E3cm), pperp=p3cm*sn;
      double th=std::atan2(pperp,pz)*TMath::RadToDeg(), ke=Elab-m3;
      if(ke>0&&th>=0) g->SetPoint(g->GetN(),th,ke); }
   return g;
}

class PPExplorer : public TObject {
public:
   PPExplorer(TString cache, double mEjectAmu, double mResidAmu, TString tag)
   {
      const double u=931.49401;
      fMbeam=16.0147*u; fMtarg=1.007825*u; fMeject=mEjectAmu*u; fMresid=mResidAmu*u; fTag=tag;
      TFile *f=TFile::Open(cache);
      if(!f||f->IsZombie()){ std::cerr<<"cannot open "<<cache<<"\n"; return; }
      TTree *t=(TTree*)f->Get("pk");
      if(!t){ std::cerr<<"no tree pk\n"; return; }
      float ke,theta,vz,chi2ndf,ic; int run;
      t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta); t->SetBranchAddress("vz",&vz);
      t->SetBranchAddress("chi2ndf",&chi2ndf); t->SetBranchAddress("ic",&ic);
      Long64_t N=t->GetEntries();
      fKe.reserve(N); fTh.reserve(N); fVz.reserve(N); fC2.reserve(N); fIc.reserve(N);
      for(Long64_t i=0;i<N;++i){ t->GetEntry(i);
         fKe.push_back(ke); fTh.push_back(theta); fVz.push_back(vz); fC2.push_back(chi2ndf); fIc.push_back(ic); }
      std::cout<<"PPExplorer: loaded "<<fKe.size()<<" protons from "<<cache<<"\n";
      MakeGui();
      Redraw();
   }

   void Redraw()
   {
      double Eb=fEbeam->GetNumber(), c2max=fChi2->GetNumber();
      double thLo=fThLo->GetNumber(), thHi=fThHi->GetNumber();
      double icLo=fIcLo->GetNumber(), icHi=fIcHi->GetNumber();
      int exb=(int)fExBins->GetNumber(); double exLo=fExLo->GetNumber(), exHi=fExHi->GetNumber();

      int thLabB=(int)fThLabBins->GetNumber(), keB=(int)fKEBins->GetNumber(), thCmB=(int)fThCMBins->GetNumber();
      double keMax=fKEMax->GetNumber(), keCutLo=fKECutLo->GetNumber(), keCutHi=fKECutHi->GetNumber();
      delete gROOT->FindObject("hEx"); delete gROOT->FindObject("hKT");
      delete gROOT->FindObject("hEt"); delete gROOT->FindObject("hZE");
      TH1F *hEx=new TH1F("hEx",Form("%s E_{x} (E_{beam}=%.0f);E_{x} [MeV];ejectiles",fTag.Data(),Eb),exb,exLo,exHi);
      TH2F *hKT=new TH2F("hKT","KE vs #theta_{lab};#theta_{lab} [deg];KE [MeV]",thLabB,0,95,keB,0,keMax);
      TH2F *hEt=new TH2F("hEt","E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]",thCmB,0,180,exb,exLo,exHi);
      TH2F *hZE=new TH2F("hZE","vertex z vs E_{x};E_{x} [MeV];vertex z [mm]",exb,exLo,exHi,110,-50,1000);

      long n=0;
      for(size_t i=0;i<fKe.size();++i){
         if(fC2[i]>c2max) continue;
         if(fTh[i]<thLo||fTh[i]>thHi) continue;
         if(fIc[i]<icLo||fIc[i]>icHi) continue;
         if(fKe[i]<keCutLo||fKe[i]>keCutHi) continue; // KE physics cut (affects Ex too)
         double thr=fTh[i]*TMath::DegToRad();
         auto [ex,thcm]=kine2b(fMbeam,fMtarg,fMeject,fMresid,Eb,thr,fKe[i]);
         hKT->Fill(fTh[i],fKe[i]);
         if(!std::isnan(ex)){ hEx->Fill(ex); hEt->Fill(thcm,ex); hZE->Fill(ex,fVz[i]); }
         ++n;
      }
      double mean=0,sig=0;
      if(hEx->GetEntries()>50){ hEx->Fit("gaus","Q0","",exLo,std::min(exHi,3.0));
         if(hEx->GetFunction("gaus")){ mean=hEx->GetFunction("gaus")->GetParameter(1); sig=hEx->GetFunction("gaus")->GetParameter(2); } }
      fLabel->SetText(Form("N=%ld   elastic peak: mean=%.3f  sigma=%.3f MeV",n,mean,sig));

      fCanvas->cd(1); gPad->SetLogy(0); hEx->Draw("hist");
      if(hEx->GetFunction("gaus")){ hEx->GetFunction("gaus")->SetNpx(500); hEx->GetFunction("gaus")->SetLineColor(kRed+1); hEx->GetFunction("gaus")->Draw("same"); }
      fCanvas->cd(2); hKT->Draw("colz");
      for(auto*g:fLines) delete g; fLines.clear();
      if(fKinLines->IsDown()){ // elastic (red) + excited-state reference lines at current Ebeam
         struct L{double ex;Color_t col;int st;};
         for(L l:{L{0,kRed+1,1},L{2,kGray+2,2},L{4,kGray+2,2},L{6,kGray+2,2}}){
            TGraph*g=kinLine(Eb,l.ex,fMbeam,fMtarg,fMeject,fMresid,l.col,l.st);
            if(g->GetN()>0){ g->Draw("L same"); fLines.push_back(g); } else delete g;
         }
         auto*tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.04);
         tx->SetTextColor(kRed+1);  tx->DrawLatex(0.5,0.85,"E_{x}=0 (g.s.)");
         tx->SetTextColor(kGray+2); tx->DrawLatex(0.5,0.80,"E_{x}=2,4,6 MeV");
      }
      fCanvas->cd(3); hEt->Draw("colz");
      fCanvas->cd(4); hZE->Draw("colz");
      fCanvas->Modified(); fCanvas->Update(); gSystem->ProcessEvents();
   }
   void Save(){ TString p="/tmp/explore_pp.png"; fCanvas->SaveAs(p); std::cout<<"saved "<<p<<"\n"; }

private:
   TGNumberEntry* mkNum(TGCompositeFrame*bar,const char*lab,double val,double lo,double hi,int dig=0){
      bar->AddFrame(new TGLabel(bar,lab),new TGLayoutHints(kLHintsLeft|kLHintsCenterY,6,2,3,3));
      auto*ne=new TGNumberEntry(bar,val,6,-1, dig?TGNumberFormat::kNESRealTwo:TGNumberFormat::kNESInteger,
                                TGNumberFormat::kNEAAnyNumber, TGNumberFormat::kNELLimitMinMax,lo,hi);
      bar->AddFrame(ne,new TGLayoutHints(kLHintsLeft,2,4,3,3)); return ne;
   }
   void MakeGui(){
      auto*main=new TGMainFrame(gClient->GetRoot(),1700,820); main->SetWindowName("16C(p,p) explorer");
      auto*bar=new TGHorizontalFrame(main);
      fEbeam=mkNum(bar,"Ebeam",192,50,400,1); fChi2=mkNum(bar,"chi2/ndf<",5,0,1000,1);
      fThLo=mkNum(bar,"thetaLo",10,0,180,1); fThHi=mkNum(bar,"thetaHi",90,0,180,1);
      fIcLo=mkNum(bar,"ICmin",950,0,5000,1); fIcHi=mkNum(bar,"ICmax",1350,0,5000,1);
      main->AddFrame(bar,new TGLayoutHints(kLHintsTop|kLHintsExpandX));
      auto*bar2=new TGHorizontalFrame(main);
      fExBins=mkNum(bar2,"ExBins",200,10,2000); fExLo=mkNum(bar2,"ExLo",-5,-20,0,1); fExHi=mkNum(bar2,"ExHi",20,1,50,1);
      auto*bR=new TGTextButton(bar2,"  Redraw  "); bR->Connect("Clicked()","PPExplorer",this,"Redraw()");
      bar2->AddFrame(bR,new TGLayoutHints(kLHintsLeft,12,4,3,3));
      auto*bS=new TGTextButton(bar2,"  Save PNG  "); bS->Connect("Clicked()","PPExplorer",this,"Save()");
      bar2->AddFrame(bS,new TGLayoutHints(kLHintsLeft,4,4,3,3));
      fKinLines=new TGCheckButton(bar2,"kin lines"); fKinLines->SetState(kButtonDown);
      fKinLines->Connect("Clicked()","PPExplorer",this,"Redraw()");
      bar2->AddFrame(fKinLines,new TGLayoutHints(kLHintsLeft|kLHintsCenterY,10,4,3,3));
      fLabel=new TGLabel(bar2,"                                                            ");
      bar2->AddFrame(fLabel,new TGLayoutHints(kLHintsLeft|kLHintsCenterY,16,4,3,3));
      main->AddFrame(bar2,new TGLayoutHints(kLHintsTop|kLHintsExpandX));
      auto*bar3=new TGHorizontalFrame(main);  // KE physics cut + 2D-plot binning
      fKECutLo=mkNum(bar3,"KEcutLo",0,0,500,1); fKECutHi=mkNum(bar3,"KEcutHi",1000,0,2000,1);
      fThLabBins=mkNum(bar3,"#theta_{lab} bins",100,10,1000); fKEBins=mkNum(bar3,"KE bins",120,10,1000);
      fKEMax=mkNum(bar3,"KE max",45,5,500,1); fThCMBins=mkNum(bar3,"#theta_{cm} bins",90,10,1000);
      main->AddFrame(bar3,new TGLayoutHints(kLHintsTop|kLHintsExpandX));
      auto*ec=new TRootEmbeddedCanvas("ec_pp",main,1680,740);
      main->AddFrame(ec,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));
      fCanvas=ec->GetCanvas(); fCanvas->Divide(2,2);
      main->MapSubwindows(); main->Resize(main->GetDefaultSize()); main->MapWindow();
   }
   std::vector<float> fKe,fTh,fVz,fC2,fIc;
   double fMbeam,fMtarg,fMeject,fMresid; TString fTag;
   TCanvas*fCanvas{nullptr};
   TGNumberEntry *fEbeam{nullptr},*fChi2{nullptr},*fThLo{nullptr},*fThHi{nullptr},*fIcLo{nullptr},*fIcHi{nullptr},
                 *fExBins{nullptr},*fExLo{nullptr},*fExHi{nullptr},
                 *fThLabBins{nullptr},*fKEBins{nullptr},*fKEMax{nullptr},*fThCMBins{nullptr},
                 *fKECutLo{nullptr},*fKECutHi{nullptr};
   TGCheckButton*fKinLines{nullptr};
   std::vector<TGraph*> fLines;
   TGLabel*fLabel{nullptr};
   ClassDef(PPExplorer,0);
};

/// (p,p): defaults. (p,d) 16C(p,d)15C: ejectile=deuteron(2.014102), residual=15C(15.010599):
///   root -l 'pp/explore_pp.C("/tmp/pd_kin.root",2.014102,15.010599,"16C(p,d)15C")'
void explore_pp(TString cache="/tmp/pp_kin.root", double mEjectAmu=1.007825, double mResidAmu=16.0147,
                TString tag="16C(p,p)")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   new PPExplorer(cache, mEjectAmu, mResidAmu, tag);
}
