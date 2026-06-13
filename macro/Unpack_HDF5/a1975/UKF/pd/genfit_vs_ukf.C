/// @file genfit_vs_ukf.C
/// @brief Overlay the 15C(p,d) g.s./0.74 doublet: GENFIT vs UKF, same runs, same
/// PID gate. GenFit kinematics from reco_gf + PRA from reco; UKF from reco_pd.
/// Each spectrum self-calibrated (g.s.->0) and peak-normalized. Saves pd/plots/genfit_vs_ukf.png
///   root -b -q 'pd/genfit_vs_ukf.C'

#include <map>
static double omega2(double x, double y, double z)
{ return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
static double exOf(double m1,double m2,double m3,double m4,double Kp,double th,double Ke)
{
   double Et1=Kp+m1,Et3=Ke+m3; double s=m1*m1+m2*m2+2*m2*Et1, uu=m2*m2+m3*m3-2*m2*Et3;
   double m4e=std::sqrt((std::cos(th)*omega2(s,m1*m1,m2*m2)*omega2(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2);
   return m4e-m4;
}
static double gsPosOf(TH1F *h)
{
   double mx=h->GetMaximum();
   TF1 dg("dgp","[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x",-1.5,2.3);
   dg.SetParameters(0.45,0.32,0.6*mx,0.6*mx,0.1*mx,0.0);
   dg.SetParLimits(0,-0.3,1.0); dg.SetParLimits(1,0.10,0.55);
   h->Fit(&dg,"QRN"); return dg.GetParameter(0);
}

void genfit_vs_ukf(double Ebeam=192.0, double chi2Cut=5.0, int nbins=70)
{
   gStyle->SetOptStat(0);
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtTools.so");
   const double u=931.49401, mC16=16.0147*u, mp=1.007825*u, md=2.01410178*u, mC15=15.0105993*u;
   auto pid=AtTools::AtParticleID::LoadJSON("pid/deuteron_band.json");
   AtTools::AtSpyralPID spy; spy.SetBField(2.85);
   const char* runs[]={"run_0106","run_0107","run_0108","run_0109","run_0110","run_0111","run_0112","run_0113"};

   TH1F *hg=new TH1F("hg","",240,-6,14); hg->SetDirectory(nullptr);
   TH1F *hu=new TH1F("hu","",240,-6,14); hu->SetDirectory(nullptr);

   // --- GENFIT: reco_gf (kinematics) + reco (PRA for gate) ---
   for (auto run : runs) {
      TString gf=TString("/mnt/f/a1975/reco_gf/")+run+"_genfit.root";
      TString rc=TString("/mnt/f/a1975/reco/")+run+"_reco.root";
      if (gSystem->AccessPathName(gf)||gSystem->AccessPathName(rc)) continue;
      TFile *fg=TFile::Open(gf); TTree *tg=(TTree*)fg->Get("cbmsim");
      TFile *fr=TFile::Open(rc); TTree *tr=(TTree*)fr->Get("cbmsim");
      TClonesArray *te=nullptr,*pe=nullptr;
      tg->SetBranchAddress("AtTrackingEvent",&te); tr->SetBranchAddress("AtPatternEvent",&pe);
      Long64_t N=std::min(tg->GetEntries(),tr->GetEntries());
      for (Long64_t i=0;i<N;++i){ tg->GetEntry(i); if(te->GetEntries()==0)continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev)continue;
         tr->GetEntry(i); if(pe->GetEntries()==0)continue; auto*pat=(AtPatternEvent*)pe->At(0); if(!pat)continue;
         std::vector<AtTrack>&cand=pat->GetTrackCand(); std::map<int,AtTrack*> byID; for(auto&t2:cand)byID[t2.GetTrackID()]=&t2;
         for(auto&ft:ev->GetFittedTracks()){ if(!ft)continue; auto&k=ft->GetKinematics();
            double ndf=ft->GetTrackMetadata()->GetNdf(),chi2=ft->GetTrackMetadata()->GetChi2(),c2n=ndf>0?chi2/ndf:1e9;
            if(k.kineticEnergy<=0||k.kineticEnergy>1000||c2n>chi2Cut)continue;
            auto it=byID.find(ft->GetTrackID()); if(it==byID.end())continue;
            auto r=spy.Estimate(*it->second); if(!r.valid||!pid.IsInside(r.sqrtdEdx,r.brho))continue;
            // GenFit theta is in DEGREES (UKF is radians) -> convert for exOf.
            double ex=exOf(mC16,mp,md,mC15,Ebeam,k.theta*TMath::DegToRad(),k.kineticEnergy); if(!std::isnan(ex))hg->Fill(ex);
         } }
      fg->Close(); fr->Close();
   }

   // --- UKF: reco_pd _ukf_d (kinematics + PRA both in tracking event) ---
   for (auto run : runs) {
      TString uf=TString("/mnt/f/a1975/reco_pd/")+run+"_ukf_d.root";
      if (gSystem->AccessPathName(uf)) continue;
      TFile *fu=TFile::Open(uf); TTree *tu=(TTree*)fu->Get("cbmsim");
      TClonesArray *te=nullptr; tu->SetBranchAddress("AtTrackingEvent",&te);
      for (Long64_t i=0;i<tu->GetEntries();++i){ tu->GetEntry(i); if(te->GetEntries()==0)continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev)continue;
         std::vector<AtTrack> orig=ev->GetTrackArray(); std::map<int,AtTrack*> byID; for(auto&t2:orig)byID[t2.GetTrackID()]=&t2;
         for(auto&ft:ev->GetFittedTracks()){ if(!ft)continue; auto&k=ft->GetKinematics();
            double ndf=ft->GetTrackMetadata()->GetNdf(),chi2=ft->GetTrackMetadata()->GetChi2(),c2n=ndf>0?chi2/ndf:1e9;
            if(k.kineticEnergy<=0||k.kineticEnergy>1000||c2n>chi2Cut)continue;
            auto it=byID.find(ft->GetTrackID()); if(it==byID.end())continue;
            auto r=spy.Estimate(*it->second); if(!r.valid||!pid.IsInside(r.sqrtdEdx,r.brho))continue;
            double ex=exOf(mC16,mp,md,mC15,Ebeam,k.theta,k.kineticEnergy); if(!std::isnan(ex))hu->Fill(ex);
         } }
      fu->Close();
   }

   // self-calibrate each to g.s.->0, fill display hists peak-normalized
   double gsg=gsPosOf(hg), gsu=gsPosOf(hu);
   TH1F *dg=new TH1F("dg2","",nbins,-2,5); dg->SetDirectory(nullptr);
   TH1F *du=new TH1F("du2","",nbins,-2,5); du->SetDirectory(nullptr);
   for(int b=1;b<=hg->GetNbinsX();++b) dg->Fill(hg->GetBinCenter(b)-gsg, hg->GetBinContent(b));
   for(int b=1;b<=hu->GetNbinsX();++b) du->Fill(hu->GetBinCenter(b)-gsu, hu->GetBinContent(b));
   // integral-normalize each to unit area -> fair shape comparison despite the yield gap
   dg->Scale(1.0/std::max(1.0,dg->Integral()));
   du->Scale(1.0/std::max(1.0,du->Integral()));

   // fit each self-calibrated doublet (common-sigma, fixed-0.74 sep, linear bg) for FWHM + overlay curve
   auto fitDoublet=[&](TH1F*h,int col)->TF1*{
      double mx=h->GetMaximum();
      TF1*f=new TF1(Form("ft_%d",col),"[2]*exp(-0.5*((x-[0])/[1])^2)+[3]*exp(-0.5*((x-[0]-0.74)/[1])^2)+[4]+[5]*x",-1.6,2.4);
      f->SetParameters(0.0,0.28,0.6*mx,0.6*mx,0.05*mx,0.0);
      f->SetParLimits(0,-0.4,0.5); f->SetParLimits(1,0.12,0.55);
      h->Fit(f,"QRN"); f->SetLineColor(col); f->SetLineWidth(3); f->SetNpx(500); return f;
   };
   TF1*fg=fitDoublet(dg,kRed+1); TF1*fu=fitDoublet(du,kAzure+2);
   double fwg=2.3548*fg->GetParameter(1), fwu=2.3548*fu->GetParameter(1);

   TCanvas *c=new TCanvas("c","gfvsukf",900,650);
   double ymax=std::max(dg->GetMaximum(),du->GetMaximum())*1.2;
   du->SetLineColor(kGray+1); du->SetLineWidth(1); du->SetFillColorAlpha(kAzure+1,0.12);
   du->SetTitle("^{15}C(p,d) g.s./0.74 doublet: GenFit vs UKF (same runs, same PID gate);E_{x}-E_{x}^{gs} [MeV];fraction / bin (area-normalized)");
   du->SetMaximum(ymax); du->GetXaxis()->SetRangeUser(-1.5,3.0);
   du->Draw("hist"); dg->SetLineColor(kRed+1); dg->SetLineWidth(1); dg->SetFillColorAlpha(kRed,0.10); dg->Draw("hist same");
   fu->Draw("same"); fg->Draw("same"); // bold fit curves on top
   TLegend *lg=new TLegend(0.50,0.70,0.88,0.88);
   lg->AddEntry(fg,Form("GenFit: FWHM %.2f MeV  (n=%.0f)",fwg,hg->Integral()),"l");
   lg->AddEntry(fu,Form("UKF: FWHM %.2f MeV  (n=%.0f)",fwu,hu->Integral()),"l");
   lg->Draw();
   c->SaveAs("pd/plots/genfit_vs_ukf.png");
   printf("genfit FWHM=%.3f (n=%.0f)  ukf FWHM=%.3f (n=%.0f)  -> saved pd/plots/genfit_vs_ukf.png\n",
          fwg,hg->Integral(),fwu,hu->Integral());
}
