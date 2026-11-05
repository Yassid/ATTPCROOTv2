// Is Spyral's dE/dx a LINEAR rescale of ATTPCROOT's, or a non-linear stretch?  A power law with
// exponent > 1 would spread the bands apart and explain better separation without any change in
// intrinsic resolution.  Match tracks on (run, event) -- same physical track, both codes.
void scalecmp(){
  gStyle->SetOptStat(0);
  std::map<std::pair<int,int>,std::pair<float,float>> A;   // (run,ev) -> (sqrtdedx, brho)
  {
    TFile f("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root");
    TTree* t=(TTree*)f.Get("pts");
    float x,y,ic; int np,run,ev;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    t->SetBranchAddress("run",&run);t->SetBranchAddress("event",&ev);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      if(ic<1045||ic>1225||np!=1||x<=0||y<=0) continue;
      auto k=std::make_pair(run,ev);
      auto it=A.find(k); if(it==A.end()||y>it->second.second) A[k]={x,y}; }
    printf("ATTPCROOT tracks (one per event, highest Brho): %zu\n", A.size());
  }
  auto* h2=new TH2F("h2","Spyral vs ATTPCROOT sqrt(dE/dx), matched tracks;ATTPCROOT #sqrt{dE/dx};Spyral #sqrt{dE/dx}",
                    180,0,45, 180,0,80);
  h2->SetDirectory(nullptr);
  long m=0;
  std::vector<double> lx, ly;
  {
    TFile f("/home/yassid/Data/spyral/C15d_a2091/points_spyral.root");
    TTree* t=(TTree*)f.Get("pts");
    float x,y,ic; int np,run;
    t->SetBranchAddress("sqrtdedx",&x);t->SetBranchAddress("brho",&y);
    t->SetBranchAddress("ic",&ic);t->SetBranchAddress("npulse",&np);
    t->SetBranchAddress("run",&run);
    TBranch* be=t->GetBranch("event");
    int ev=0; if(be) t->SetBranchAddress("event",&ev);
    if(!be){ printf("NO event branch on the Spyral plane -- cannot match\n"); return; }
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      if(ic<1000||ic>1250||np!=1||x<=0||y<=0) continue;
      auto it=A.find(std::make_pair(run,ev)); if(it==A.end()) continue;
      // require the two to agree on Brho, else they are different tracks of the same event
      if(fabs(y-it->second.second)/it->second.second > 0.10) continue;
      m++; h2->Fill(it->second.first, x);
      if(it->second.first>1&&x>1){ lx.push_back(log(it->second.first)); ly.push_back(log(x)); } }
  }
  printf("matched tracks (same event AND Brho within 10%%): %ld\n", m);
  if(lx.size()>1000){
    double sx=0,sy=0,sxx=0,sxy=0; int n=lx.size();
    for(int i=0;i<n;i++){ sx+=lx[i]; sy+=ly[i]; sxx+=lx[i]*lx[i]; sxy+=lx[i]*ly[i]; }
    double a=(n*sxy-sx*sy)/(n*sxx-sx*sx), b=(sy-a*sx)/n;
    printf("\n  log-log fit:  ln(Spyral) = %.4f * ln(ATTPCROOT) + %.4f\n", a, b);
    printf("  -> Spyral = %.3f * ATTPCROOT^%.4f\n", exp(b), a);
    printf("  exponent %.4f : %s\n", a, fabs(a-1)<0.02? "LINEAR rescale (no stretch)" :
           (a>1? "NON-LINEAR, Spyral STRETCHES the scale -> bands spread apart"
               : "NON-LINEAR, Spyral COMPRESSES the scale"));
  }
  auto* c=new TCanvas("c","",800,750); gPad->SetLogz();
  h2->Draw("colz");
  c->SaveAs("/home/yassid/C15d_dedx_scale.png");
  printf("wrote /home/yassid/C15d_dedx_scale.png\n");
}
