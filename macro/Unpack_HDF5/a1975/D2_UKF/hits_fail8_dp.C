/// @file hits_fail8_dp.C -- hits of census_arc_dp.C bound candidates (vz 200-400, thcm lo-hi) that AtSpyralPID
/// rejects with failCode fc, data only. Per candidate: x-y and rho-z, candidate red over the event's other tracks.
void hits_fail8_dp(int fc=8, double lo=15, double hi=30, int nmax=8, TString png="hits_fail8_dp.png", int skip=0){
  gSystem->Load("libAtReconstruction.so"); FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel=kError; gStyle->SetOptStat(0);
  TChain d("c"); d.Add("/mnt/f/a1975/census_arc/data_*.root");
  TString cut=Form("thl>89&&thl<161&&ke>1.2&&ex>-1.5&&ex<1.5&&vz>=200&&vz<400&&tcm>=%g&&tcm<%g&&%s",lo,hi,fc>=0?Form("!sv&&sf==%d",fc):"sv");
  d.Draw("run:ev:tid:d0:nh:ke",cut,"goff"); int n=d.GetSelectedRows(); printf("%d candidates\n",n);
  std::vector<std::array<double,6>> L; for(int i=0;i<n;i++) L.push_back({d.GetVal(0)[i],d.GetVal(1)[i],d.GetVal(2)[i],d.GetVal(3)[i],d.GetVal(4)[i],d.GetVal(5)[i]});
  int np=std::min(nmax,n-skip); auto c=new TCanvas("c","c",1600,400*((np+1)/2)); c->Divide(4,(np+1)/2);
  for(int k=0;k<np;k++){ auto&e=L[k*std::max(1,(n-skip)/np)+skip]; int run=e[0]; Long64_t ev=e[1]; int tid=e[2];
    TFile f(Form("/mnt/f/a1975/reco_d2_dv1104/run_%04d_multifit_reco.root",run)); auto t=(TTree*)f.Get("cbmsim"); TClonesArray*pa=nullptr;
    t->SetBranchStatus("*",0); t->SetBranchStatus("AtPatternEvent*",1); t->SetBranchAddress("AtPatternEvent",&pa); t->GetEntry(ev); auto pe=(AtPatternEvent*)pa->At(0);
    auto gxy=new TGraph, grz=new TGraph, oxy=new TGraph, orz=new TGraph;
    for(auto&tr:pe->GetTrackCand()) for(auto&h:tr.GetHitArray()){ auto p=h->GetPosition(); bool me=tr.GetTrackID()==tid;
      (me?gxy:oxy)->AddPoint(p.X(),p.Y()); (me?grz:orz)->AddPoint(p.Z(),std::hypot(p.X(),p.Y())); }
    c->cd(2*k+1); auto fr=gPad->DrawFrame(-250,-250,250,250,Form("run %d ev %lld tid %d  d0 %.0f nh %.0f KE %.2f;x;y",run,ev,tid,e[3],e[4],e[5]));
    oxy->SetMarkerStyle(7); oxy->SetMarkerColor(kGray+1); oxy->Draw("P"); gxy->SetMarkerStyle(7); gxy->SetMarkerColor(kRed); gxy->Draw("P");
    c->cd(2*k+2); gPad->DrawFrame(0,0,1000,260,";z_{hit};#rho"); orz->SetMarkerStyle(7); orz->SetMarkerColor(kGray+1); orz->Draw("P"); grz->SetMarkerStyle(7); grz->SetMarkerColor(kRed); grz->Draw("P");
    f.Close(); }
  c->SaveAs(png);
}
