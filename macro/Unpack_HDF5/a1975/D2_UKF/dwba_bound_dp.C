/// @file dwba_bound_dp.C -- compare dp_bound_angdist_pertheta.csv with the FRESCO ADWA bound-state curves
/// (~/fair_install/16C_dp/fresco/inputs/dp17C_bound_adwa*.nin; unit cfp amplitudes, so fitted scales are S).
/// state1 = g.s. 3/2+ n(1d3/2) | state2 = 0.217 1/2+ n(2s1/2) | state3 = 0.331 5/2+ n(1d5/2)
/// Theory is AVERAGED over each 5-deg bin with sin(theta) weight before comparing.
TGraph* rd(const char*f){ auto g=new TGraph; std::ifstream in(f); double a,b; while(in>>a>>b) g->AddPoint(a,b); return g; }
double avg(TGraph*g,double a,double b){ double s=0,w=0; for(double t=a;t<=b;t+=0.05){ double ww=sin(t*TMath::DegToRad()); s+=ww*g->Eval(t); w+=ww;} return s/w; }
void dwba_bound_dp(const char*csv="dp_bound_angdist_pertheta.csv"){
  std::vector<double> lo,hi,y,e;
  { std::ifstream in(csv); std::string l; std::map<std::string,int> c;
    while(std::getline(in,l)){ if(l.empty()||l[0]=='#') continue; std::vector<std::string> t; std::stringstream ss(l); std::string k; while(std::getline(ss,k,',')) t.push_back(k);
      if(l[0]=='t'){ for(size_t i=0;i<t.size();++i) c[t[i]]=i; continue; }
      if(std::stod(t[c["use"]])<0.5) continue;
      lo.push_back(std::stod(t[c["thcm_lo"]])); hi.push_back(std::stod(t[c["thcm_hi"]])); y.push_back(std::stod(t[c["dsdo_mb_sr"]])); e.push_back(std::stod(t[c["err_mb_sr"]])); } }
  const int N=y.size();
  const char*TAG[2]={"dp17C_bound_adwa","dp17C_bound_adwa_23p10"}; const char*EL[2]={"E_d 23.55","E_d 23.10"};
  TString base="/home/yassid/fair_install/16C_dp/fresco/outputs/";
  TGraph*G[2][3];
  for(int k=0;k<2;++k) for(int s=0;s<3;++s) G[k][s]=rd(base+TAG[k]+Form("_state%d.dat",s+1));
  for(int k=0;k<2;++k){
    std::vector<double> gs(N),s0(N),d5(N);
    for(int i=0;i<N;++i){ gs[i]=avg(G[k][0],lo[i],hi[i]); s0[i]=avg(G[k][1],lo[i],hi[i]); d5[i]=avg(G[k][2],lo[i],hi[i]); }
    printf("\n=== %s  (bin-averaged theory, mb/sr, unit S) ===\n  thcm    data          gs(l=2 d3/2)  0.217(l=0)  0.331(l=2 d5/2)   sum(S=1)\n",EL[k]);
    for(int i=0;i<N;++i) printf("  %2.0f-%2.0f  %5.2f +- %4.2f   %7.3f      %7.3f     %7.3f          %7.3f\n",lo[i],hi[i],y[i],e[i],gs[i],s0[i],d5[i],gs[i]+s0[i]+d5[i]);
    // linear least squares fits
    auto fit1=[&](std::vector<double>&m,const char*lab){ double sxy=0,sxx=0; for(int i=0;i<N;++i){ sxy+=y[i]*m[i]/(e[i]*e[i]); sxx+=m[i]*m[i]/(e[i]*e[i]); }
      double S=sxy/sxx, eS=1/sqrt(sxx), c2=0; for(int i=0;i<N;++i) c2+=pow((y[i]-S*m[i])/e[i],2);
      printf("  %-34s S = %6.3f +- %.3f   chi2/ndf = %5.1f/%d\n",lab,S,eS,c2,N-1); };
    std::vector<double> dd(N),all(N); for(int i=0;i<N;++i){ dd[i]=gs[i]+d5[i]; all[i]=gs[i]+s0[i]+d5[i]; }
    fit1(s0,"pure l=0 (0.217 shape)"); fit1(dd,"pure l=2 (gs+0.331, common S)"); fit1(all,"all three, common S");
    // 2-param: a*s0 + b*dd
    TMatrixD A(2,2); TVectorD B(2); for(int i=0;i<N;++i){ double w=1/(e[i]*e[i]); A(0,0)+=s0[i]*s0[i]*w; A(0,1)+=s0[i]*dd[i]*w; A(1,1)+=dd[i]*dd[i]*w; B(0)+=y[i]*s0[i]*w; B(1)+=y[i]*dd[i]*w; }
    A(1,0)=A(0,1); TMatrixD C=A; C.Invert(); TVectorD P=C*B; double c2=0; for(int i=0;i<N;++i) c2+=pow((y[i]-P(0)*s0[i]-P(1)*dd[i])/e[i],2);
    printf("  %-34s S(l=0) = %6.3f +- %.3f  S(l=2) = %6.3f +- %.3f  corr %+.2f  chi2/ndf = %5.1f/%d\n","l=0 + l=2 free",P(0),sqrt(C(0,0)),P(1),sqrt(C(1,1)),C(0,1)/sqrt(C(0,0)*C(1,1)),c2,N-2);
    if(k==0){ // plot: SHAPE comparison -- each curve scaled by its own single-S fit (the 2-param fit goes unphysical, S(l=0)<0)
      auto S1=[&](std::vector<double>&m){ double a=0,b=0; for(int i=0;i<N;++i){ a+=y[i]*m[i]/(e[i]*e[i]); b+=m[i]*m[i]/(e[i]*e[i]); } return a/b; };
      double Ss=S1(s0), Sd=S1(dd);
      gStyle->SetOptStat(0); TCanvas c("c","",1500,650); c.Divide(2,1);
      c.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetGrid();
      TH1F*fr=gPad->DrawFrame(0,0.3,60,40); fr->SetTitle("^{16}C(d,p)^{17}C bound states vs ADWA (E_{d}=23.55 MeV): SHAPES;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]");
      TGraphErrors gd; for(int i=0;i<N;++i){ gd.SetPoint(i,0.5*(lo[i]+hi[i]),y[i]); gd.SetPointError(i,0.5*(hi[i]-lo[i]),e[i]); }
      auto gl0=new TGraph, gl2=new TGraph;
      for(int i=0;i<G[0][0]->GetN();++i){ double x=G[0][0]->GetX()[i]; if(x>60) break; gl0->AddPoint(x,Ss*G[0][1]->GetY()[i]); gl2->AddPoint(x,Sd*(G[0][0]->GetY()[i]+G[0][2]->GetY()[i])); }
      gl0->SetLineColor(kRed+1); gl0->SetLineWidth(2); gl0->Draw("L"); gl2->SetLineColor(kGreen+2); gl2->SetLineWidth(2); gl2->Draw("L");
      gd.SetMarkerStyle(20); gd.SetMarkerSize(1.2); gd.SetMarkerColor(kBlue+1); gd.SetLineColor(kBlue+1); gd.SetLineWidth(2); gd.Draw("P");
      TLegend L(0.15,0.13,0.62,0.30); L.SetTextSize(0.032);
      L.AddEntry(&gd,"data, per-#theta v_{z} ranges (stat.)","pe");
      L.AddEntry(gl0,Form("l=0 (0.217 s_{1/2}) #times %.2f",Ss),"l"); L.AddEntry(gl2,Form("l=2 (g.s. d_{3/2} + 0.331 d_{5/2}) #times %.3f",Sd),"l"); L.Draw();
      c.cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.13); gPad->SetGrid();
      TH1F*f2=gPad->DrawFrame(10,0.002,45,50); f2->SetTitle("data / theory (unit S), bin-averaged;#theta_{cm} [deg];data / #sigma_{DWBA}(S=1)");
      TGraphErrors r0,r2; for(int i=0;i<N;++i){ double x=0.5*(lo[i]+hi[i]); r0.SetPoint(i,x,y[i]/s0[i]); r0.SetPointError(i,2.5,e[i]/s0[i]); r2.SetPoint(i,x+0.4,y[i]/dd[i]); r2.SetPointError(i,2.5,e[i]/dd[i]); }
      r0.SetMarkerStyle(21); r0.SetMarkerColor(kRed+1); r0.SetLineColor(kRed+1); r0.Draw("P"); r2.SetMarkerStyle(20); r2.SetMarkerColor(kGreen+2); r2.SetLineColor(kGreen+2); r2.Draw("P");
      TLegend L2(0.15,0.75,0.6,0.88); L2.SetTextSize(0.032); L2.AddEntry(&r0,"vs l=0 (0.217)","pe"); L2.AddEntry(&r2,"vs l=2 (g.s.+0.331)","pe"); L2.Draw();
      c.SaveAs("dp_bound_vs_dwba.png"); }
  }
}
