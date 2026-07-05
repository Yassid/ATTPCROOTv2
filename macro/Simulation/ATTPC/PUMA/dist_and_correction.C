/// @file dist_and_correction.C
/// @brief (1) Delta p/p DISTRIBUTIONS per pion momentum for UKF and genfit, and
///        (2) the deterministic material correction before/after (genfit median bias
///        flattened toward 0). The correction adds back Delta_p(p) = 18.2 * S(p)/S(375),
///        S(p)=BetheBloch(p)/beta — the same formula AtGenfitter::SetVertexMaterial-
///        Correction applies in-pipeline, here shown offline on the fitted momenta.
/// Run: root -b -q dist_and_correction.C
double medv(std::vector<double>v){if(v.empty())return 0;std::sort(v.begin(),v.end());return v[v.size()/2];}
double iqrv(std::vector<double>v){if(v.size()<4)return 0;std::sort(v.begin(),v.end());return(v[3*v.size()/4]-v[v.size()/4])/1.349;}
const double MPI=139.57039;
double bbShape(double p){ double me=0.510999,I=40e-6,m=MPI; double bg=p/m,b2=bg*bg/(1+bg*bg),g2=1+bg*bg;
  double Tmax=2*me*bg*bg/(1+2*sqrt(g2)*me/m+(me/m)*(me/m)); return (1/b2)*(0.5*log(2*me*b2*g2*Tmax/(I*I))-b2); }
double sfun(double p){ double b=p/sqrt(p*p+MPI*MPI); return b>0? bbShape(p)/b : 0; }
double dpCorr(double p){ const double refP=375,refDp=18.2; double s0=sfun(refP); if(s0<=0)return 0;
  double d1=refDp*sfun(p)/s0; return refDp*sfun(p+d1)/s0; } // one iteration: eval at corrected p

// fill matched-fit momenta for one fitter; returns raw fitted |p| (MeV) per track
void loadP(TString digi,TString sim,int fitter,std::vector<double>&outP){
  TFile fD(digi); TTree*tD=(TTree*)fD.Get("cbmsim");
  TFile fS(sim);  TTree*tS=(TTree*)fS.Get("cbmsim"); if(!tD||!tS)return;
  const char*br=fitter==0?"AtTrackingEventUKF":"AtTrackingEventGenfit";
  TClonesArray*fe=new TClonesArray("AtTrackingEvent"); tD->SetBranchAddress(br,&fe);
  TClonesArray*pat=new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent",&pat);
  TClonesArray*mcP=new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint",&mcP);
  Long64_t n=std::min(tD->GetEntries(),tS->GetEntries());
  for(Long64_t e=0;e<n;++e){tD->GetEntry(e);tS->GetEntry(e);
    if(!pat->GetEntries()||!fe->GetEntries())continue;
    int nMC=mcP->GetEntries(); std::vector<double>mx(nMC),my(nMC); std::vector<int>mid(nMC);
    for(int k=0;k<nMC;++k){auto*p=(AtMCPoint*)mcP->At(k);mx[k]=p->GetX()*10;my[k]=p->GetY()*10;mid[k]=p->GetTrackID();}
    auto&tc=((AtPatternEvent*)pat->At(0))->GetTrackCand();
    for(const auto&ft:((AtTrackingEvent*)fe->At(0))->GetFittedTracks()){
      double KE=ft->GetKinematics(0).kineticEnergy; if(!(KE>0))continue;
      int tid=ft->GetTrackID(); if(tid<0||tid>=(int)tc.size())continue;
      std::map<int,int>votes; for(auto&h:tc[tid].GetHitArray()){auto&p3=h->GetPosition();double best=9;int bid=-99;
        for(int k=0;k<nMC;++k){double d2=(p3.X()-mx[k])*(p3.X()-mx[k])+(p3.Y()-my[k])*(p3.Y()-my[k]);if(d2<best){best=d2;bid=mid[k];}}
        if(bid!=-99)votes[bid]++;}
      int bv=0,bt=-99;for(auto&kv:votes)if(kv.second>bv){bv=kv.second;bt=kv.first;}
      if(!(bt==0||bt==1))continue;
      outP.push_back(sqrt(KE*KE+2*KE*MPI));
    }
  }
}

void dist_and_correction(){
  gSystem->Load("libAtReconstruction.so"); gStyle->SetOptStat(0);
  struct P{const char*tag;double p0;int col;};
  std::vector<P> pts={{"p150",150,kViolet+1},{"p200",200,kAzure+1},{"p300",300,kGreen+2},{"p375",375,kOrange+7},{"p500",500,kRed+1}};
  TString D="/mnt/f/puma_sweep/",S="./data/";

  // ---- (1) distributions: 2x3 grid, one panel per energy, UKF vs genfit ----
  auto*c1=new TCanvas("dist","distributions",1400,820); c1->Divide(3,2);
  std::vector<double> vp, medG_before, medG_after, sigG, deficitPct;
  for(size_t i=0;i<pts.size();++i){
    TString dg=D+"output_digi_es_"+pts[i].tag+".root", sm=S+"attpcsim_"+pts[i].tag+".root";
    std::vector<double> pu,pg; loadP(dg,sm,0,pu); loadP(dg,sm,1,pg);
    double p0=pts[i].p0;
    c1->cd(i+1); gPad->SetGrid();
    auto*hu=new TH1F(Form("hu%zu",i),Form("p = %.0f MeV/c;#Deltap/p [%%];tracks (norm)",p0),60,-60,40);
    auto*hg=new TH1F(Form("hg%zu",i),"",60,-60,40);
    std::vector<double> dpu,dpg;
    for(double p:pu){double d=100*(p-p0)/p0; hu->Fill(d); dpu.push_back(d);}
    for(double p:pg){double d=100*(p-p0)/p0; hg->Fill(d); dpg.push_back(d);}
    if(hu->Integral()>0)hu->Scale(1./hu->Integral()); if(hg->Integral()>0)hg->Scale(1./hg->Integral());
    hu->SetLineColor(kBlue+1);hu->SetLineWidth(2); hg->SetLineColor(kRed+1);hg->SetLineWidth(2);
    hu->SetMaximum(1.3*std::max(hu->GetMaximum(),hg->GetMaximum())); hu->Draw("hist");hg->Draw("hist same");
    TLine*l0=new TLine(0,0,0,hu->GetMaximum());l0->SetLineStyle(2);l0->SetLineColor(kGray+2);l0->Draw();
    TLatex L;L.SetNDC();L.SetTextSize(0.05);
    L.SetTextColor(kBlue+1);L.DrawLatex(0.15,0.84,Form("UKF #sigma=%.0f%%",iqrv(dpu)));
    L.SetTextColor(kRed+1);L.DrawLatex(0.15,0.77,Form("genfit #sigma=%.0f%%",iqrv(dpg)));
    // correction bookkeeping (genfit)
    std::vector<double> dpB,dpA; for(double p:pg){dpB.push_back(100*(p-p0)/p0); double pc=p+dpCorr(p); dpA.push_back(100*(pc-p0)/p0);}
    vp.push_back(p0); medG_before.push_back(medv(dpB)); medG_after.push_back(medv(dpA)); sigG.push_back(iqrv(dpA));
  }
  c1->cd(6); // legend panel
  auto*lg=new TLegend(0.1,0.4,0.9,0.7); auto*a=new TH1F("la","",1,0,1),*b=new TH1F("lb","",1,0,1);
  a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kRed+1);b->SetLineWidth(2);
  lg->AddEntry(a,"UKF","l"); lg->AddEntry(b,"genfit","l"); lg->SetHeader("#Deltap/p per momentum");
  lg->SetTextSize(0.06); lg->Draw();
  TLatex T;T.SetNDC();T.SetTextSize(0.045);T.SetTextColor(kGray+3);
  T.DrawLatex(0.1,0.32,"dashed = truth (0%)"); T.DrawLatex(0.1,0.24,"left shift grows at low p"); T.DrawLatex(0.1,0.18,"= Cu/Al material loss");
  c1->SaveAs("./data/dist_by_energy.png");

  // ---- (2) material correction: median bias vs p, before vs after ----
  auto*c2=new TCanvas("corr","material correction",760,600); c2->SetGrid();
  auto mk=[&](std::vector<double>&y,int col,int mst){auto*g=new TGraph(vp.size(),vp.data(),y.data());
    g->SetLineColor(col);g->SetMarkerColor(col);g->SetMarkerStyle(mst);g->SetLineWidth(2);g->SetMarkerSize(1.5);return g;};
  auto*gb=mk(medG_before,kRed+1,21),*ga=mk(medG_after,kGreen+2,20);
  auto*mg=new TMultiGraph();mg->Add(gb);mg->Add(ga);
  mg->SetTitle("genfit vertex-momentum bias: deterministic material correction;pion momentum [MeV/c];median #Deltap/p [%]");
  mg->Draw("ALP"); mg->GetYaxis()->SetRangeUser(-22,8);
  TLine*z=new TLine(150,0,500,0);z->SetLineStyle(2);z->SetLineColor(kGray+2);z->Draw();
  auto*lg2=new TLegend(0.5,0.18,0.88,0.34);lg2->AddEntry(gb,"before (in-gas momentum)","lp");
  lg2->AddEntry(ga,"after material correction","lp");lg2->Draw();
  c2->SaveAs("./data/material_correction.png");

  printf("\n  p     genfit bias before -> after (material correction)\n");
  for(size_t i=0;i<vp.size();++i) printf("  %-4.0f   %+6.1f%%  ->  %+6.1f%%   (res %.1f%%)\n",vp[i],medG_before[i],medG_after[i],sigG[i]);
  printf("  wrote dist_by_energy.png + material_correction.png\n\n");
}
