/// far_ex_dp.C (2026-10-18): Ex by vz region, thcm 10-20, data (proton gate) vs far sim, Spyral kinematics. Run from D2_UKF.
void far_ex_dp(){
 TChain d("c"),s("c"); d.Add("/mnt/f/a1975/census_arc/data_*.root"); s.Add("/mnt/f/a1975_C16_dp_sim_far/far_s41*_census.root");
 std::ifstream in("pid/proton_d2_dv1104.json"); std::string str((std::istreambuf_iterator<char>(in)),{}); TGraph G;
 auto p=str.find("vertices"); p=str.find('[',p)+1;
 while(true){auto a=str.find('[',p); if(a==std::string::npos)break; auto b=str.find(']',a); double x,y; sscanf(str.substr(a+1,b-a-1).c_str(),"%lf , %lf",&x,&y); G.AddPoint(x,y); p=b+1; if(str[str.find_first_not_of(" \n\t,",p)]==']')break;}
 auto cg=new TCutG("gp",G.GetN(),G.GetX(),G.GetY()); cg->SetVarX("ssq"); cg->SetVarY("sbr");
 TString W="sv&&sth>89&&sth<161&&ske>1.2&&stcm>=10&&stcm<20";
 const char* R[4]={"svz>=200&&svz<400","svz>=900&&svz<980","svz>=980&&svz<1010","svz>=1010&&svz<1050"};
 TH1D* h[8]; for(int i=0;i<4;i++){ h[i]=new TH1D(Form("d%d",i),"",16,-3,5); h[i+4]=new TH1D(Form("s%d",i),"",16,-3,5);
  d.Draw(Form("sex>>d%d",i),W+"&&gp&&"+R[i],"goff"); s.Draw(Form("sex>>s%d",i),W+"&&"+R[i],"goff");}
 printf("Ex bin    | DATA 200-400  900-980  980-1010  1010-1050 | SIM 200-400 900-980 980-1010 1010-1050\n");
 for(int b=1;b<=16;b++){ printf("%5.1f-%4.1f |",h[0]->GetBinLowEdge(b),h[0]->GetBinLowEdge(b)+0.5);
  for(int i=0;i<4;i++)printf(" %8.0f",h[i]->GetBinContent(b)); printf("   |"); for(int i=4;i<8;i++)printf(" %8.0f",h[i]->GetBinContent(b)); printf("\n");}
}
