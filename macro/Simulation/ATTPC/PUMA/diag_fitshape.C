/// @file diag_fitshape.C
/// @brief Diagnose the AtPSAMultiFit fit: fit the GET response A*exp(-3u)sin(u)u^3
///        (u=(t-t0)/tau, tau=25 TB fixed, as the macro uses) to real pad traces,
///        and report the FIT UNCERTAINTY on t0 (impulse time) vs the peak. Shows
///        why z(t0) is noisy: the u^3 onset is nearly flat -> t0 weakly constrained.
/// Run: root -b -q 'diag_fitshape.C("data/output_digi_pi1000_maxring.root")'
double Resp(double t, double t0, double tau){ double u=(t-t0)/tau; if(u<=0)return 0; return exp(-3*u)*sin(u)*u*u*u; }

void diag_fitshape(TString file="data/output_digi_pi1000_maxring.root",
                   TString out="/Users/quantumlab/fair_install/puma_slides/figs/psa_fitshape.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62,"xyz"); gStyle->SetTitleFont(62,"xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const double tau=25.0, uPeak=1.166, peakLag=uPeak*tau; // TB
   const double mmPerTB=0.3; // 20ns * 1.5cm/us

   TFile f(file); auto*t=(TTree*)f.Get("cbmsim"); TClonesArray*raw=nullptr; t->SetBranchAddress("AtRawEvent",&raw);

   std::vector<double> errT0_TB, errPeak_TB, redChi2;
   auto*c=new TCanvas("c","",1200,850); c->Divide(2,2); int panel=0;
   for(Long64_t e=0;e<200 && panel<4;e++){ t->GetEntry(e);
      auto*ev=(AtRawEvent*)(raw?raw->At(0):nullptr); if(!ev)continue;
      // rank pads by peak amplitude
      std::vector<std::pair<double,const AtPad*>> pads;
      for(const auto&p:ev->GetPads()){double mx=0;for(double a:p->GetADC())mx=std::max(mx,a);pads.push_back({mx,p.get()});}
      std::sort(pads.begin(),pads.end(),[](auto&a,auto&b){return a.first>b.first;});
      for(auto&pr:pads){ if(pr.first<50)continue; const auto&adc=pr.second->GetADC();
         int pk=0; double pkv=0; for(int i=0;i<512;i++)if(adc[i]>pkv){pkv=adc[i];pk=i;}
         if(pk<40||pk>480)continue;
         int lo=std::max(1,pk-50),hi=std::min(510,pk+60);
         // noise from a quiet region
         double noise=1; { std::vector<double>q; for(int i=20;i<80;i++)q.push_back(adc[i]);
            std::sort(q.begin(),q.end()); double med=q[q.size()/2],s=0; for(double x:q)s+=fabs(x-med); noise=std::max(1.0,1.4826*s/q.size()); }
         auto*g=new TGraphErrors();
         for(int i=lo;i<=hi;i++){int n=g->GetN(); g->SetPoint(n,i,adc[i]); g->SetPointError(n,0,noise);}
         // fixed-tau single-component fit: params [A, t0]
         auto*fit=new TF1(Form("f%lld",e),[tau](double*x,double*p){return p[0]*Resp(x[0],p[1],tau);},lo,hi,2);
         fit->SetParameters(pkv/0.0441, pk-peakLag); fit->SetParLimits(0,0,10*pkv/0.0441+10);
         fit->SetParLimits(1,pk-peakLag-tau,pk-peakLag+tau);
         auto r=g->Fit(fit,"SQNR");
         double et0=fit->GetParError(1); double chi=r->Ndf()>0?r->Chi2()/r->Ndf():0;
         errT0_TB.push_back(et0); redChi2.push_back(chi);
         // draw the first 4
         if(panel<4){ c->cd(panel+1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
            g->SetTitle(Form("pad %d;time bucket;ADC",pr.second->GetPadNum()));
            g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetMarkerColor(kAzure+2); g->SetLineColor(kAzure+2); g->Draw("AP");
            fit->SetLineColor(kRed+1); fit->SetLineWidth(3); fit->Draw("same");
            double t0f=fit->GetParameter(1);
            auto*l0=new TLine(t0f,0,t0f,pkv*0.4); l0->SetLineColor(kGreen+2); l0->SetLineStyle(2); l0->SetLineWidth(2); l0->Draw();
            auto*lp=new TLine(t0f+peakLag,0,t0f+peakLag,pkv); lp->SetLineColor(kOrange+8); lp->SetLineStyle(2); lp->SetLineWidth(2); lp->Draw();
            if(panel==0){auto*lg=new TLegend(0.5,0.68,0.88,0.9); lg->SetTextFont(62); lg->SetTextSize(0.045);
               lg->AddEntry(g,"ADC trace","p"); lg->AddEntry(fit,"GET-response fit","l");
               lg->AddEntry(l0,Form("fitted t_{0} (#pm%.0f TB)",et0),"l"); lg->AddEntry(lp,"peak (t_{0}+29 TB)","l"); lg->Draw();}
            panel++;
         }
         break; // one pad per event for the plot
      }
   }
   c->SaveAs(out);
   auto md=[](std::vector<double>v){std::sort(v.begin(),v.end());return v.empty()?0:v[v.size()/2];};
   printf("\nFIT DIAGNOSTIC (single-component, tau fixed):\n");
   printf("  median fit error on t0   : %.1f TB  = %.2f mm\n", md(errT0_TB), md(errT0_TB)*mmPerTB);
   printf("  median reduced chi2      : %.2f  (>>1 = poor shape match / multi-pulse)\n", md(redChi2));
   printf("  (peak position is the same fit's t0+29TB but robust; err on peak << err on t0)\n");
   printf("FITSHAPE_DONE\n");
}
