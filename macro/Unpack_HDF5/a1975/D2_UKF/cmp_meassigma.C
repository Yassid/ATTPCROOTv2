/// Compare the KE-vs-theta locus across measSigma values against the (d,t) g.s. kinematics.
/// The test: does <KE_meas> at FORWARD angles come down toward the 0.85-1.3 MeV the
/// kinematics demand, or stay pinned at the ~2.2-2.9 MeV floor? Binned in theta, never in
/// KE -- binning on the noisy variable manufactures a slope from the scatter alone.
void cmp_meassigma(const char* runs="run_0016_multifit,run_0019_multifit,run_0022_multifit"){
  const double u=931.49401,m1=16.0147013*u,m2=2.0135532*u,m3=3.01550072*u,m4=15.0105993*u,Eb=184.17;
  double E1=Eb+m1,pb=std::sqrt(E1*E1-m1*m1);
  const int NT=240; double thg[NT],kel[NT]; int n=0;
  for(double th=1;th<=56;th+=0.25){ double t=th*TMath::DegToRad(); double kl=-1;
    for(double ke=0.05;ke<140;ke+=0.005){ double E3=ke+m3,p3=std::sqrt(E3*E3-m3*m3),E4=E1+m2-E3;
      double px=p3*std::sin(t),pz=p3*std::cos(t);
      double m4x2=E4*E4-(px*px+(pb-pz)*(pb-pz)); if(m4x2<0) continue;
      if(std::fabs(std::sqrt(m4x2)-m4)<0.04){ kl=ke; break; } }
    if(kl>0&&n<NT){ thg[n]=th; kel[n]=kl; ++n; } }
  auto kin=[&](double th)->double{ if(th<thg[0]||th>thg[n-1]) return -1;
    for(int i=1;i<n;i++) if(thg[i]>=th){ double f=(th-thg[i-1])/(thg[i]-thg[i-1]);
      return kel[i-1]+f*(kel[i]-kel[i-1]); } return -1; };
  const char* dirs[3]={"/mnt/f/a1975/gf_dt_ms40/","/mnt/f/a1975/gf_dt_ms20/","/mnt/f/a1975/gf_dt_ms10/"};
  const char* labs[3]={"measSigma 4.0 (production)","measSigma 2.0","measSigma 1.0"};
  const int NB=9; double lo=12,w=5;
  printf("\n  KINEMATICS require: ");
  for(int b=0;b<NB;b++) printf("%6.2f",kin(lo+(b+0.5)*w)); printf("   <- target\n");
  auto*cv=new TCanvas("cv","",1100,700); auto*mg=new TMultiGraph();
  int col[3]={kBlack,kBlue+1,kRed+1};
  for(int d=0;d<3;d++){
    double sum[NB]={0}; long cnt[NB]={0}; long tot=0;
    TObjArray*a=TString(runs).Tokenize(","); 
    for(int ir=0;ir<a->GetEntries();++ir){
      TString fn=TString(dirs[d])+((TObjString*)a->At(ir))->GetString()+"_genfitter_t.root";
      if(gSystem->AccessPathName(fn)) continue;
      TFile*f=TFile::Open(fn); if(!f||f->IsZombie()) continue;
      TTree*t=(TTree*)f->Get("cbmsim"); if(!t){f->Close();continue;}
      TClonesArray*te=nullptr; t->SetBranchAddress("AtTrackingEvent",&te);
      for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
        if(te->GetEntries()==0) continue; auto*ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
        for(auto&ft:ev->GetFittedTracks()){ if(!ft) continue;
          auto&k=ft->GetKinematicsXtr(); auto&m=ft->GetTrackMetadata(); if(!m) continue;
          double ndf=m->GetNdf(); if(ndf<=0||m->GetChi2()/ndf>5) continue;
          double ke=k.kineticEnergy, th=k.theta*TMath::RadToDeg();
          double kk=kin(th); if(kk<0||std::fabs(ke-kk)>4) continue;
          int b=(int)((th-lo)/w); if(b<0||b>=NB) continue;
          sum[b]+=ke; ++cnt[b]; ++tot; } }
      f->Close(); }
    printf("  %-26s ",labs[d]);
    auto*g=new TGraph();
    for(int b=0;b<NB;b++){ if(cnt[b]<20){printf("   ---");continue;}
      double mean=sum[b]/cnt[b]; printf("%6.2f",mean);
      g->SetPoint(g->GetN(),lo+(b+0.5)*w,mean); }
    printf("   N=%ld\n",tot);
    g->SetLineColor(col[d]); g->SetMarkerColor(col[d]); g->SetMarkerStyle(20); g->SetLineWidth(2);
    g->SetTitle(labs[d]); mg->Add(g,"LP"); }
  auto*gk=new TGraph(); for(int b=0;b<NB;b++) gk->SetPoint(b,lo+(b+0.5)*w,kin(lo+(b+0.5)*w));
  gk->SetLineColor(kGreen+2); gk->SetLineWidth(3); gk->SetTitle("(d,t) g.s. kinematics"); mg->Add(gk,"L");
  mg->SetTitle("<KE> vs #theta_{lab} for each measSigma;#theta_{lab} [deg];<KE> [MeV]");
  mg->Draw("A"); cv->BuildLegend(0.55,0.68,0.88,0.88);
  cv->SaveAs("/tmp/claude-1000/-home-yassid/b025789e-3ded-4e0d-8d25-35b205d047eb/scratchpad/meassigma.png");
}
