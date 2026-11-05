// Does the real integral buy separation at the cost of WIDTH?  Resolving power is what matters:
// R = separation / (sigma1 + sigma2). Fit the deuteron and triton peaks in the ln-mass spectrum
// for peak-charge vs real-integral on the SAME tracks.
void tauwidth(){
  gSystem->Load("libAtReconstruction.so");
  const char* runs[2]={"run_0016_reco.root","run_0019_reco.root"};
  for(int mode=0;mode<2;mode++){
    AtTools::AtSpyralPID e; e.SetBField(2.85); e.SetMinPoints(30); e.SetZTieTolerance(0.0);
    e.SetUseTraceIntegral(mode==1);
    std::vector<double> v;
    for(int k=0;k<2;k++){
      TFile* F=TFile::Open(Form("/home/yassid/Data/C15d_tautest/flo2/%s",runs[k]));
      if(!F||F->IsZombie()) continue;
      TTree* t=(TTree*)F->Get("cbmsim");
      TClonesArray* pa=nullptr; t->SetBranchAddress("AtPatternEvent",&pa);
      for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
        if(!pa||pa->GetEntriesFast()==0) continue;
        auto* ev=dynamic_cast<AtPatternEvent*>(pa->At(0)); if(!ev) continue;
        for(auto& tr : ev->GetTrackCand()){
          auto r=e.Estimate(const_cast<AtTrack&>(tr));
          if(!r.valid||r.brho<=0||r.sqrtdEdx<=0) continue;
          if(std::hypot(r.vertex.X(),r.vertex.Y())>=25) continue;
          v.push_back(r.sqrtdEdx*r.brho); } }
      F->Close();
    }
    std::sort(v.begin(),v.end()); double med=v[v.size()/2];
    auto* h=new TH1D(Form("tw%d",mode),"",240,-1.3,1.3); h->SetDirectory(nullptr);
    for(double q:v) h->Fill(log(q/med));
    // fit d at 0 and the t bump at its observed place
    double tpos = (mode? 0.53 : 0.42);
    TF1 fd("fd","gaus(0)+pol1(3)",-0.26,0.22);
    fd.SetParameters(h->GetMaximum(),0,0.10,10,0);
    fd.SetParLimits(2,0.02,0.30); h->Fit(&fd,"RQ0N");
    TF1 ft("ft","gaus(0)+pol1(3)",tpos-0.17,tpos+0.17);
    ft.SetParameters(h->GetBinContent(h->FindBin(tpos)),tpos,0.08,10,0);
    ft.SetParLimits(1,tpos-0.10,tpos+0.10); ft.SetParLimits(2,0.02,0.25);
    h->Fit(&ft,"RQ0N");
    double sd=fabs(fd.GetParameter(2)), st=fabs(ft.GetParameter(2));
    double sep=ft.GetParameter(1)-fd.GetParameter(1);
    printf("%-26s  d: mu %+.3f sigma %.4f | t: mu %+.3f sigma %.4f | separation %.3f | R = %.2f\n",
           mode? "REAL INTEGRAL":"PEAK (production)",
           fd.GetParameter(1), sd, ft.GetParameter(1), st, sep, sep/(sd+st));
    delete h;
  }
}
