/// @file ex_vtxcorr.C
/// @brief 16C(p,p) Ex with a VERTEX ENERGY-LOSS correction: the 16C beam loses energy
/// penetrating the gas (~0.027 MeV/mm), so the beam energy at the reaction vertex is
/// vertex-z dependent. Using a constant nominal Ebeam smears/shifts Ex with vertex z.
/// Here Ebeam(z) = CATIMA energy of 16C after traveling from the entrance to the vertex
/// (AtELossCATIMA::GetEnergy), anchored so the elastic peak sits at Ex~0 on average.
/// Recomputes Ex per event and compares to the uncorrected spectrum.
///
///   root -b -q 'pp/ex_vtxcorr.C("run_0106","/tmp/run_0106_genfitter_pp.root")'

#include <map>
#include <tuple>
#include <vector>

static double omega2(double x, double y, double z) { return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
static double Ex_2b(double m1,double m2,double m3,double m4,double Kp,double thlab,double Ke){
   double Et1=Kp+m1,Et3=Ke+m3; double s=m1*m1+m2*m2+2*m2*Et1, uu=m2*m2+m3*m3-2*m2*Et3;
   double arg=(std::cos(thlab)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)
              -(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
   if(arg<0) return NAN; return std::sqrt(arg)-m4;
}

void ex_vtxcorr(TString run="run_0106", TString fitFile="/tmp/run_0106_genfitter_pp.root",
                TString inDir="/mnt/f/a1975/reco/", TString gateFile="pid/proton_band.json",
                double Ebeam=192.0, double icMin=950, double icMax=1350, int icTbLo=1000, int icTbHi=1350,
                double chi2Cut=5.0, double bField=2.85, double thMin=10, double thMax=170, double densGcm3=8.3147e-5)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   const double u=931.49401, m_C16=16.0147*u, m_p=1.007825*u;
   auto pid=AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   auto beamEloss=new AtTools::AtELossCATIMA(densGcm3);
   beamEloss->SetProjectile(16,6,16.0147);
   std::vector<std::tuple<int,int,int>> mat; mat.push_back(std::make_tuple(1,1,1)); beamEloss->SetMaterial(mat);

   TFile *ff=TFile::Open(fitFile);              TTree *tf=(TTree*)ff->Get("cbmsim");
   TFile *fc=TFile::Open(inDir+run+"_FRIB.root"); TTree *tc=(TTree*)fc->Get("cbmsim");
   TClonesArray *te=nullptr,*re=nullptr;
   tf->SetBranchAddress("AtTrackingEvent",&te);
   tc->SetBranchAddress("AtRawEvent",&re);

   // pass 1: collect selected protons (ke, theta, vertex z) + nominal Ex
   std::vector<double> vKe,vTh,vZ,vExNom;
   Long64_t N=std::min(tf->GetEntries(),tc->GetEntries());
   for(Long64_t i=0;i<N;++i){
      tc->GetEntry(i); double ic=-1;
      if(re->GetEntries()>0){ auto*raw=(AtRawEvent*)re->At(0);
         if(raw&&!raw->GetGenTraces().empty()){ auto&adc=raw->GetGenTraces()[0]->GetADC();
            double mx=-1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();++b) mx=std::max(mx,adc[b]); ic=mx; } }
      if(ic<icMin||ic>icMax) continue;
      tf->GetEntry(i); if(te->GetEntries()==0) continue;
      auto*ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      std::vector<AtTrack> orig=ev->GetTrackArray(); std::map<int,AtTrack*> byID;
      for(auto&tr:orig) byID[tr.GetTrackID()]=&tr;
      for(auto&ft:ev->GetFittedTracks()){
         if(!ft) continue; auto&k=ft->GetKinematics(); auto&m=ft->GetTrackMetadata(); if(!m) continue;
         double ndf=m->GetNdf(), c2n=ndf>0?m->GetChi2()/ndf:1e9, ke=k.kineticEnergy, thDeg=k.theta*TMath::RadToDeg();
         if(ke<=0||ke>1000||c2n>chi2Cut||thDeg<thMin||thDeg>thMax) continue;
         auto it=byID.find(ft->GetTrackID()); if(it==byID.end()) continue;
         auto r=spy.Estimate(*it->second); if(!r.valid||!pid.IsInside(r.sqrtdEdx,r.brho)) continue;
         double ex=Ex_2b(m_C16,m_p,m_p,m_C16,Ebeam,k.theta,ke); if(std::isnan(ex)) continue;
         vKe.push_back(ke); vTh.push_back(k.theta); vZ.push_back(ft->GetVertex(0).Z()); vExNom.push_back(ex);
      }
   }
   printf("selected protons: %zu\n",vKe.size());
   if(vKe.empty()) return;

   double zmin=*std::min_element(vZ.begin(),vZ.end()), zmax=*std::max_element(vZ.begin(),vZ.end());
   printf("vertex z range [%.0f, %.0f] mm\n",zmin,zmax);

   // The beam-eloss systematic only shows in the beam-energy-SENSITIVE forward protons
   // (large angles are kinematically insensitive to Ebeam). Fit Ex_nom(z) for forward
   // elastic events -> z0 where the nominal beam energy is correct (Ex=0). The beam
   // ENTERS at low z and loses energy toward high z (Ex rises with z, confirmed).
   double zref=0; for(double z:vZ) zref+=z; zref/=vZ.size(); // reference (mean) vertex z
   // forward (beam-energy-sensitive) elastic events carry the systematic
   std::vector<size_t> fwd;
   for(size_t j=0;j<vKe.size();++j)
      if(vTh[j]*TMath::RadToDeg()<45 && std::fabs(vExNom[j])<4) fwd.push_back(j);

   // linear beam-energy profile E(z) = Ebeam - k*(z - zref); self-calibrate k by
   // flattening the forward-elastic Ex slope (avoids over-correcting; k should land
   // near the CATIMA 16C dE/dx if vertex z tracks beam penetration).
   auto fwdSlope=[&](double k)->double{
      if(fwd.size()<5) return 0;
      TGraph g; for(size_t j:fwd){ double eb=Ebeam-k*(vZ[j]-zref);
         double ex=Ex_2b(m_C16,m_p,m_p,m_C16,eb,vTh[j],vKe[j]); if(!std::isnan(ex)) g.SetPoint(g.GetN(),vZ[j],ex); }
      g.Fit("pol1","Q"); return g.GetFunction("pol1")->GetParameter(1);
   };
   double slopeBefore=fwdSlope(0.0), kBest=0, best=1e9;
   for(double k=-0.02;k<=0.06;k+=0.001){ double s=std::fabs(fwdSlope(k)); if(s<best){best=s;kBest=k;} }
   printf("forward-elastic slope nominal=%.5f MeV/mm; self-calibrated beam gradient k=%.4f MeV/mm "
          "(CATIMA 16C dE/dx@192=%.4f); residual slope=%.5f\n",
          slopeBefore, kBest, beamEloss->GetdEdx(Ebeam), fwdSlope(kBest));
   auto ebeamAt=[&](double z)->double{ return Ebeam - kBest*(z-zref); };

   TH1F *hNom=new TH1F("hNom","16C(p,p') Ex: vertex E-loss correction;E_{x} [MeV];protons",120,-5,15);
   TH1F *hCor=new TH1F("hCor","corrected",120,-5,15);
   hNom->SetLineColor(kAzure+2); hNom->SetLineWidth(2);
   hCor->SetLineColor(kRed+1);   hCor->SetLineWidth(2);
   TH2F *h2n=new TH2F("h2n","Ex_{nom} vs vertex z;z [mm];E_{x}",100,zmin,zmax,100,-5,10);
   TH2F *h2c=new TH2F("h2c","Ex_{corr} vs vertex z;z [mm];E_{x}",100,zmin,zmax,100,-5,10);
   TGraph gFwdC;
   for(size_t j=0;j<vKe.size();++j){
      hNom->Fill(vExNom[j]); h2n->Fill(vZ[j],vExNom[j]);
      double ex=Ex_2b(m_C16,m_p,m_p,m_C16,ebeamAt(vZ[j]),vTh[j],vKe[j]);
      if(!std::isnan(ex)){ hCor->Fill(ex); h2c->Fill(vZ[j],ex);
         if(vTh[j]*TMath::RadToDeg()<45 && std::fabs(ex)<4) gFwdC.SetPoint(gFwdC.GetN(),vZ[j],ex); }
   }
   if(gFwdC.GetN()>=5){ gFwdC.Fit("pol1","Q"); auto*fn=gFwdC.GetFunction("pol1");
      printf("forward-elastic Ex slope: before=%.5f -> after=%.5f MeV/mm (chamber tilt %.2f -> %.2f MeV)\n",
             slopeBefore, fn->GetParameter(1), slopeBefore*(zmax-zmin), fn->GetParameter(1)*(zmax-zmin)); }
   auto sig=[&](TH1F*h){ h->Fit("gaus","Q0","",-3,3); return h->GetFunction("gaus")?h->GetFunction("gaus")->GetParameter(2):-1; };
   double sNom=sig(hNom), sCor=sig(hCor);
   printf("elastic peak sigma: nominal=%.3f MeV -> corrected=%.3f MeV\n",sNom,sCor);

   TCanvas *c=new TCanvas("c","vtxcorr",1400,520); c->Divide(3,1);
   c->cd(1); hNom->Draw("hist"); hCor->Draw("hist same");
   TLegend*lg=new TLegend(0.55,0.7,0.88,0.88);
   lg->AddEntry(hNom,Form("nominal #sigma=%.2f",sNom),"l");
   lg->AddEntry(hCor,Form("vtx-corr #sigma=%.2f",sCor),"l"); lg->Draw();
   c->cd(2); h2n->Draw("colz"); h2n->ProfileX()->Draw("same");
   c->cd(3); h2c->Draw("colz"); h2c->ProfileX()->Draw("same");
   c->SaveAs("/tmp/ex_vtxcorr.png");
   printf("wrote /tmp/ex_vtxcorr.png\n");
}
