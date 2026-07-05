/// @file diag_momentum_shape.C
/// @brief Is the reconstructed momentum distribution GAUSSIAN for a monoenergetic beam?
///        A wide-but-Gaussian dp/p is honest resolution; a skewed / heavy-tailed /
///        bimodal one signals a reconstruction pathology that an IQR "resolution"
///        number hides. Plots, for the monoenergetic branch-8 pi sample (truth
///        |p|=374.9 MeV/c), the distribution SHAPE at three stages:
///          1) PRA rigidity p = 0.2998*B*R (raw geometry, R = circle radius)
///          2) sagitta s = L^2/(8R) over the hit arc  (should be ~Gaussian: linear in
///             the measurement) -> tests whether non-Gaussianity enters at the R=c/s step
///          3) UKF and genfit fitted |p|
///        Each with a Gaussian fit to the core + mean/median/RMS/skew/kurtosis.
/// Run: root -b -q 'diag_momentum_shape.C("/mnt/f/puma_sweep/output_digi_pi_pid.root")'
double skewness(const std::vector<double>&v){ if(v.size()<3)return 0; double m=0; for(double x:v)m+=x; m/=v.size();
   double s2=0,s3=0; for(double x:v){double d=x-m; s2+=d*d; s3+=d*d*d;} s2/=v.size(); s3/=v.size();
   return s2>0? s3/std::pow(s2,1.5):0; }
double kurtosis(const std::vector<double>&v){ if(v.size()<4)return 0; double m=0; for(double x:v)m+=x; m/=v.size();
   double s2=0,s4=0; for(double x:v){double d=x-m; s2+=d*d; s4+=d*d*d*d;} s2/=v.size(); s4/=v.size();
   return s2>0? s4/(s2*s2)-3.0:0; } // excess kurtosis (0 = Gaussian)
double medv(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void fillHist(TH1F*h,const std::vector<double>&v){ for(double x:v) h->Fill(x); }

void diag_momentum_shape(TString digiFile="/mnt/f/puma_sweep/output_digi_pi_pid.root",
                         double B=4.0, double p0=374.9)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetOptFit(0);

   TFile f(digiFile); TTree*t=(TTree*)f.Get("cbmsim");
   if(!t){ printf("missing %s\n",digiFile.Data()); return; }
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent",&pat);
   TClonesArray*ukf=new TClonesArray("AtTrackingEvent"); t->SetBranchAddress("AtTrackingEventUKF",&ukf);
   TClonesArray*gf =new TClonesArray("AtTrackingEvent"); t->SetBranchAddress("AtTrackingEventGenfit",&gf);
   const double mPi=139.57039;

   std::vector<double> praP, sag, ukfP, gfP; // rigidity, sagitta[mm], fitted |p|
   for(Long64_t e=0;e<t->GetEntries();++e){
      t->GetEntry(e);
      if(pat->GetEntries()) for(auto&tr:((AtPatternEvent*)pat->At(0))->GetTrackCand()){
         double R=tr.GetGeoRadius();
         if(R>0&&R<1e5){ praP.push_back(0.299792458*B*R);
            // sagitta from the chord spanned by the hits: s = R - sqrt(R^2 - (c/2)^2)
            auto&hits=tr.GetHitArray(); if(hits.size()>=2){
               auto&a=hits.front()->GetPosition(); auto&b=hits.back()->GetPosition();
               double c=std::hypot(a.X()-b.X(),a.Y()-b.Y()); double arg=R*R-0.25*c*c;
               if(arg>0) sag.push_back(R-std::sqrt(arg)); } }
      }
      TClonesArray*arr[2]={ukf,gf}; std::vector<double>*dst[2]={&ukfP,&gfP};
      for(int fi=0;fi<2;++fi) if(arr[fi]->GetEntries())
         for(const auto&ft:((AtTrackingEvent*)arr[fi]->At(0))->GetFittedTracks()){
            double KE=ft->GetKinematics(0).kineticEnergy; if(KE>0) dst[fi]->push_back(std::sqrt(KE*KE+2*KE*mPi)); }
   }

   auto report=[&](const char*nm,std::vector<double>&v,double truth){
      if(v.empty()){printf("  %-10s: (empty)\n",nm); return;}
      std::vector<double> dp; for(double x:v) dp.push_back(100*(x-truth)/truth);
      double mean=0; for(double x:dp)mean+=x; mean/=dp.size();
      double rms=0; for(double x:dp)rms+=(x-mean)*(x-mean); rms=std::sqrt(rms/dp.size());
      printf("  %-10s: n=%zu  median=%+.1f%%  mean=%+.1f%%  RMS=%.1f%%  skew=%+.2f  exKurt=%+.2f\n",
             nm,v.size(),medv(dp),mean,rms,skewness(dp),kurtosis(dp));
   };
   printf("\n===== momentum-distribution SHAPE (monoenergetic pi, |p|=%.1f) =====\n",p0);
   printf("  (Gaussian => skew~0, exKurt~0; positive skew/kurt => right/heavy tail)\n");
   report("PRA rigid",praP,p0); report("UKF |p|",ukfP,p0); report("genfit |p|",gfP,p0);
   double sTrue=medv(sag); // just for reference scale
   printf("  sagitta   : median=%.2f mm (truth arc; near-Gaussian expected)  skew=%+.2f\n",sTrue,skewness(sag));

   // draw: 4 panels
   auto*c=new TCanvas("ms","momentum shape",1300,900); c->Divide(2,2);
   auto mkDP=[&](std::vector<double>&v,const char*ti,int col){
      auto*h=new TH1F(Form("h_%s",ti),Form("%s;#Deltap/p [%%];tracks",ti),80,-60,60);
      for(double x:v) h->Fill(100*(x-p0)/p0); h->SetLineColor(col); h->SetLineWidth(2); return h; };
   c->cd(1); { auto*h=mkDP(praP,"PRA rigidity p=0.2998BR",kBlue+1); h->Draw();
      h->Fit("gaus","Q0","",-25,25); auto*fn=h->GetFunction("gaus");
      if(fn){fn->SetLineColor(kRed);fn->SetNpx(500);fn->Draw("same");
         TLatex L;L.SetNDC();L.SetTextSize(0.035);
         L.DrawLatex(0.14,0.85,Form("Gaus core #sigma=%.1f%%",fn->GetParameter(2)));
         L.DrawLatex(0.14,0.80,Form("skew=%+.2f exKurt=%+.2f",skewness(praP),kurtosis(praP)));} }
   c->cd(2); { auto*h=mkDP(ukfP,"UKF |p|",kGreen+2); h->Draw(); h->Fit("gaus","Q0","",-25,25);
      auto*fn=h->GetFunction("gaus"); if(fn){fn->SetLineColor(kRed);fn->SetNpx(500);fn->Draw("same");} }
   c->cd(3); { auto*h=mkDP(gfP,"genfit |p|",kMagenta+1); h->Draw(); h->Fit("gaus","Q0","",-25,25);
      auto*fn=h->GetFunction("gaus"); if(fn){fn->SetLineColor(kRed);fn->SetNpx(500);fn->Draw("same");} }
   c->cd(4); { auto*h=new TH1F("hR",";PRA radius R [mm];tracks",100,0,1200);
      for(double p:praP) h->Fill(p/(0.299792458*B)); h->SetLineColor(kBlack); h->SetLineWidth(2); h->Draw();
      TLatex L;L.SetNDC();L.SetTextSize(0.032);L.DrawLatex(0.4,0.85,Form("truth R=%.0f mm",p0/(0.299792458*B)));
      L.DrawLatex(0.4,0.80,"(1/s tail -> high R = high p)"); }
   c->SaveAs("./data/diag_momentum_shape.png");
   printf("wrote ./data/diag_momentum_shape.png\n\n");
}
