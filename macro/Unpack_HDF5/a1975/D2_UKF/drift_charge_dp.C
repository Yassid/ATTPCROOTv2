/// @file drift_charge_dp.C -- hit charge vs drift (z_hit) for backward protons, data vs g.s. sim (2026-10-12).
/// Tests whether long-drift signal loss (attachment / diffusion vs threshold) is data-only: the forward-cm bound
/// protons that data lose at vz 200-400 live entirely at z_hit 600-864, the longest drift.
/// Tracks: census_arc_dp.C rows, Spyral valid, thl 89-161, KE 1.2-4, (data) inside the proton gate.
/// Per track the hit charges are binned by z_hit and divided by the track's own median hit charge, so the gain
/// difference cancels and only the z dependence is left. Also hits per mm of z travelled, same normalisation.
void drift_charge_dp(int runMax = 25){
  gSystem->Load("libAtReconstruction.so"); FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel=kError; gStyle->SetOptStat(0);
  TGraph G; { std::ifstream in("pid/proton_d2_dv1104.json"); std::string s((std::istreambuf_iterator<char>(in)),{});
    auto p=s.find("vertices"); p=s.find('[',p)+1; while(true){ auto a=s.find('[',p); if(a==std::string::npos) break; auto b=s.find(']',a); double x,y; sscanf(s.substr(a+1,b-a-1).c_str(),"%lf , %lf",&x,&y); G.AddPoint(x,y); p=b+1; if(s[s.find_first_not_of(" \n\t,",p)]==']') break; } }
  TProfile *pq[2], *pk[2]; TH1D* hq[2][3];
  for(int k=0;k<2;k++){ pq[k]=new TProfile(Form("pq%d",k),";z_{hit} (mm);hit charge / track median",18,0,900); pk[k]=new TProfile(Form("pk%d",k),";z_{hit};hits per mm of z / track mean",18,0,900);
    for(int j=0;j<3;j++) hq[k][j]=new TH1D(Form("hq%d%d",k,j),";hit charge / track median",60,0,3); }
  auto doFile=[&](int k,TString cen,TString reco){
    TFile fc(cen); auto c=(TTree*)fc.Get("c"); int ev,tid,sv; double sth,ske,ssq,sbr;
    c->SetBranchAddress("ev",&ev); c->SetBranchAddress("tid",&tid); c->SetBranchAddress("sv",&sv); c->SetBranchAddress("sth",&sth);
    c->SetBranchAddress("ske",&ske); c->SetBranchAddress("ssq",&ssq); c->SetBranchAddress("sbr",&sbr);
    std::map<int,std::set<int>> want;
    for(Long64_t i=0;i<c->GetEntries();i++){ c->GetEntry(i); if(!sv||sth<89||sth>161||ske<1.2||ske>4) continue; if(k==0&&!G.IsInside(ssq,sbr)) continue; want[ev].insert(tid); }
    TFile fr(reco); auto t=(TTree*)fr.Get("cbmsim"); TClonesArray*pa=nullptr; t->SetBranchStatus("*",0); t->SetBranchStatus("AtPatternEvent*",1); t->SetBranchAddress("AtPatternEvent",&pa);
    for(auto&[e,ids]:want){ t->GetEntry(e); auto pe=(AtPatternEvent*)pa->At(0); if(!pe) continue;
      for(auto&tr:pe->GetTrackCand()){ if(!ids.count(tr.GetTrackID())) continue; auto&hs=tr.GetHitArray(); if(hs.size()<20) continue;
        std::vector<double> q; double z0=1e9,z1=-1e9; for(auto&h:hs){ q.push_back(h->GetCharge()); z0=std::min(z0,h->GetPosition().Z()); z1=std::max(z1,h->GetPosition().Z()); }
        auto qs=q; std::nth_element(qs.begin(),qs.begin()+qs.size()/2,qs.end()); double med=qs[qs.size()/2]; if(med<=0||z1-z0<50) continue;
        std::map<int,int> nb; for(auto&h:hs){ double z=h->GetPosition().Z(); pq[k]->Fill(z,h->GetCharge()/med); nb[int(z/50)]++;
          int j=z<300?0:(z<600?1:2); hq[k][j]->Fill(h->GetCharge()/med); }
        double mean=double(hs.size())/(z1-z0);
        for(int b=int(z0/50)+1;b<int(z1/50);b++) pk[k]->Fill(b*50+25,(nb[b]/50.)/mean); } } };
  for(int r=16;r<=runMax;r++){ TString cen=Form("/mnt/f/a1975/census_arc/data_%04d.root",r); if(gSystem->AccessPathName(cen)) continue;
    doFile(0,cen,Form("/mnt/f/a1975/reco_d2_dv1104/run_%04d_multifit_reco.root",r)); }
  doFile(1,"/mnt/f/a1975/census_arc/sim_gs_s4001.root","/mnt/f/a1975_C16_dp_sim/gs_s4001_reco_g75000.root");
  printf(" z_hit   | charge/median  data   sim  | hits-per-mm/mean  data   sim\n");
  for(int b=1;b<=18;b++) printf(" %3.0f-%3.0f |              %.3f  %.3f |                  %.3f  %.3f\n",(b-1)*50.,b*50.,pq[0]->GetBinContent(b),pq[1]->GetBinContent(b),pk[0]->GetBinContent(b),pk[1]->GetBinContent(b));
  auto cv=new TCanvas("cv","cv",1500,450); cv->Divide(3,1);
  cv->cd(1); pq[0]->SetLineColor(kBlack); pq[1]->SetLineColor(kRed); pq[0]->SetMinimum(0); pq[0]->Draw(); pq[1]->Draw("same");
  auto lg=new TLegend(0.15,0.75,0.5,0.88); lg->AddEntry(pq[0],"data (gated p)","l"); lg->AddEntry(pq[1],"sim g.s.","l"); lg->Draw();
  cv->cd(2); pk[0]->SetLineColor(kBlack); pk[1]->SetLineColor(kRed); pk[0]->SetMinimum(0); pk[0]->Draw(); pk[1]->Draw("same");
  cv->cd(3); const char* nm[3]={"z<300","300-600","600+"}; int st[3]={1,2,3};
  for(int k=0;k<2;k++) for(int j=0;j<3;j++){ auto h=hq[k][j]; if(h->Integral()>0) h->Scale(1/h->Integral()); h->SetLineColor(k?kRed:kBlack); h->SetLineStyle(st[j]); h->SetTitle(nm[j]); h->Draw(k||j?"hist same":"hist"); }
  cv->SaveAs("drift_charge_dp.png");
}
