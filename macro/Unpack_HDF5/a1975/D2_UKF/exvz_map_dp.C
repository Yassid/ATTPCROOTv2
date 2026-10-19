/// exvz_map_dp.C (2026-10-19): Ex vs vz tables per thcm bin, data (proton gate) vs far sim. Run from D2_UKF.
void exvz_map_dp(){
 TChain d("c"),s("c"); d.Add("/mnt/f/a1975/census_arc/data_*.root"); s.Add("/mnt/f/a1975_C16_dp_sim_far/far_s41*_census.root");
 std::ifstream in("pid/proton_d2_dv1104.json"); std::string str((std::istreambuf_iterator<char>(in)),{}); TGraph G;
 auto p=str.find("vertices"); p=str.find('[',p)+1;
 while(true){auto a=str.find('[',p); if(a==std::string::npos)break; auto b=str.find(']',a); double x,y; sscanf(str.substr(a+1,b-a-1).c_str(),"%lf , %lf",&x,&y); G.AddPoint(x,y); p=b+1; if(str[str.find_first_not_of(" \n\t,",p)]==']')break;}
 auto cg=new TCutG("gp",G.GetN(),G.GetX(),G.GetY()); cg->SetVarX("ssq"); cg->SetVarY("sbr");
 TString W="sv&&sth>89&&sth<161&&ske>1.2";
 double T[4]={10,20,30,40};
 for(int t=0;t<3;t++) for(int k=0;k<2;k++){
  TH2D h(Form("h%d%d",t,k),"",9,100,1000,20,-2,3);
  TString c=W+Form("&&stcm>=%g&&stcm<%g",T[t],T[t+1]); if(k==0)c+="&&gp";
  (k?s:d).Draw(Form("sex:svz>>h%d%d",t,k),c,"goff");
  printf("\n%s thcm %g-%g : rows Ex (0.25 MeV), columns vz 100..1000 by 100\n  Ex   |",k?"SIM":"DATA",T[t],T[t+1]);
  for(int x=1;x<=9;x++)printf(" %4.0f",h.GetXaxis()->GetBinLowEdge(x)); printf("\n");
  for(int y=20;y>=1;y--){printf(" %5.2f|",h.GetYaxis()->GetBinLowEdge(y)); for(int x=1;x<=9;x++)printf(" %4.0f",h.GetBinContent(x,y)); printf("\n");}
 }
}
