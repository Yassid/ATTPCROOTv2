/// @file dlc_ukf_genfit.C
/// @brief UKF vs GENFIT momentum resolution & bias vs |p| with DLC ON, winning
///        config (HDBSCAN + MultiFit + ring). Reads the high-stat reco files
///        directly and draws with statistical error bars
///        (sigma err ~ sigma/sqrt(2N); median err ~ 1.253 sigma/sqrt(N)).
/// Run: root -b -q 'dlc_ukf_genfit.C'
static double medOf(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }
static double iqrS(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }

// fill dp (%) for a fitter branch; returns N
static int grab(TString file, const char *branch, double p0MeV, std::vector<double> &dp)
{
   const double m=139.57039; TFile f(file); auto*t=(TTree*)f.Get("cbmsim"); if(!t)return 0;
   auto*arr=new TClonesArray("AtTrackingEvent"); t->SetBranchAddress(branch,&arr);
   for(Long64_t e=0;e<t->GetEntries();e++){ t->GetEntry(e); if(arr->GetEntries()==0)continue;
      for(const auto&ft:((AtTrackingEvent*)arr->At(0))->GetFittedTracks()){ const auto&k=ft->GetKinematics(0);
         double KE=k.kineticEnergy; if(!(KE>0))continue; double p=std::sqrt(KE*KE+2*KE*m);
         dp.push_back(100.*(p-p0MeV)/p0MeV); } }
   return dp.size();
}

void dlc_ukf_genfit(TString out = "/Users/quantumlab/fair_install/puma_slides/figs/dlc_ukf_genfit.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62,"xyz"); gStyle->SetTitleFont(62,"xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   struct Pt{ double p; TString f; }; // p [MeV], reco file
   std::vector<Pt> pts={{149.6,"data/hs_reco150.root"},{375.0,"data/hs_reco375.root"},{600.0,"data/hs_reco600.root"}};
   const int N=pts.size();
   std::vector<double> pv(N),sU(N),sG(N),bU(N),bG(N),esU(N),esG(N),ebU(N),ebG(N),ex(N,0);
   printf("  p[MeV]  UKF bias/sig   GENFIT bias/sig   (N)\n");
   for(int i=0;i<N;i++){ std::vector<double> du,dg; int nu=grab(pts[i].f,"AtTrackingEventUKF",pts[i].p,du);
      int ng=grab(pts[i].f,"AtTrackingEventGenfit",pts[i].p,dg);
      pv[i]=pts[i].p; sU[i]=iqrS(du); sG[i]=iqrS(dg); bU[i]=medOf(du); bG[i]=medOf(dg);
      esU[i]=nu>1?sU[i]/std::sqrt(2.0*nu):0; esG[i]=ng>1?sG[i]/std::sqrt(2.0*ng):0;
      ebU[i]=nu>1?1.253*sU[i]/std::sqrt((double)nu):0; ebG[i]=ng>1?1.253*sG[i]/std::sqrt((double)ng):0;
      printf("  %5.0f   %+5.1f/%4.1f     %+5.1f/%4.1f     (%d)\n",pv[i],bU[i],sU[i],bG[i],sG[i],nu); }

   auto *c = new TCanvas("c","",1050,500); c->Divide(2,1);
   c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGrid();
   auto *gU=new TGraphErrors(N,pv.data(),sU.data(),ex.data(),esU.data());
   auto *gG=new TGraphErrors(N,pv.data(),sG.data(),ex.data(),esG.data());
   gU->SetTitle("momentum resolution (DLC on);truth |p| [MeV/c];#sigma_{IQR}(p)/p  [%]");
   gU->SetMarkerStyle(20); gU->SetMarkerSize(1.5); gU->SetMarkerColor(kAzure+2); gU->SetLineColor(kAzure+2); gU->SetLineWidth(3);
   gU->GetYaxis()->SetRangeUser(0,48); gU->GetXaxis()->SetLimits(100,650); gU->Draw("ALP");
   gG->SetMarkerStyle(21); gG->SetMarkerSize(1.5); gG->SetMarkerColor(kRed+1); gG->SetLineColor(kRed+1); gG->SetLineWidth(3); gG->SetLineStyle(2); gG->Draw("LP");
   { auto*l=new TLine(100,15,650,15); l->SetLineColor(kGray+2); l->SetLineStyle(3); l->Draw(); }
   { auto*lg=new TLegend(0.42,0.72,0.86,0.88); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
     lg->AddEntry(gU,"UKF","lp"); lg->AddEntry(gG,"GENFIT","lp"); lg->Draw(); }
   { auto*t=new TLatex(); t->SetTextFont(62); t->SetTextSize(0.033); t->SetTextColor(kGray+3); t->DrawLatex(380,16.5,"no-DLC #sigma"); }
   c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGrid();
   auto *bUg=new TGraphErrors(N,pv.data(),bU.data(),ex.data(),ebU.data());
   auto *bGg=new TGraphErrors(N,pv.data(),bG.data(),ex.data(),ebG.data());
   bUg->SetTitle("momentum bias (DLC on);truth |p| [MeV/c];median bias [%]");
   bUg->SetMarkerStyle(20); bUg->SetMarkerSize(1.5); bUg->SetMarkerColor(kAzure+2); bUg->SetLineColor(kAzure+2); bUg->SetLineWidth(3);
   bUg->GetYaxis()->SetRangeUser(-19,6); bUg->GetXaxis()->SetLimits(100,650); bUg->Draw("ALP");
   bGg->SetMarkerStyle(21); bGg->SetMarkerSize(1.5); bGg->SetMarkerColor(kRed+1); bGg->SetLineColor(kRed+1); bGg->SetLineWidth(3); bGg->SetLineStyle(2); bGg->Draw("LP");
   { auto*l=new TLine(100,0,650,0); l->SetLineColor(kGray+2); l->SetLineStyle(3); l->Draw(); }
   { auto*lg=new TLegend(0.40,0.20,0.86,0.36); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
     lg->AddEntry(bUg,"UKF","lp"); lg->AddEntry(bGg,"GENFIT (unbiased)","lp"); lg->Draw(); }
   c->SaveAs(out);
   printf("DLC_UKFGF_DONE\n");
}
