#include "ridge.C"
// Draw the saved deuteron gate on top of the identified bands.
void gateplot(){
  gStyle->SetOptStat(0); gStyle->SetPalette(kBird);
  std::vector<double> gx,gy;
  std::ifstream in("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/deuteron_C15d.json");
  std::string all((std::istreambuf_iterator<char>(in)),std::istreambuf_iterator<char>());
  size_t p=all.find("\"vertices\"");
  while(true){ size_t a=all.find('[',p+1); if(a==std::string::npos)break;
    size_t b=all.find(']',a); if(b==std::string::npos)break;
    double u,v; if(sscanf(all.substr(a+1,b-a-1).c_str()," %lf , %lf",&u,&v)==2){gx.push_back(u);gy.push_back(v);} p=b; }
  TH2F* h=loadPlane("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mf_C15d.root",1045,1225,1.0,"gp");
  auto* c=new TCanvas("c","",900,760);
  h->SetTitle("ATTPCROOT plane: saved deuteron gate vs the identified bands");
  h->GetXaxis()->SetTitle("sqrt(dE/dx)"); h->GetYaxis()->SetTitle("B#rho [T m]");
  h->GetXaxis()->SetRangeUser(0,40); h->SetMinimum(0.5); h->SetMaximum(200);
  gPad->SetLogz(); gPad->SetRightMargin(0.13); h->Draw("colz");
  double K0=0.4242;
  struct S{const char*n;double A,z;int col;};
  S sp[4]={{"p",1,1,kWhite},{"d",2,1,kRed},{"t",3,1,kMagenta},{"#alpha",4,2,kGreen+1}};
  auto* l=new TLatex(); l->SetNDC(); l->SetTextSize(0.03);
  for(int s=0;s<4;s++){
    auto* g=new TGraph();
    for(double br=0.15;br<1.34;br+=0.006){ double v=K0*sqrt(betheRel(br,sp[s].A,sp[s].z));
      if(v>0.5&&v<40) g->SetPoint(g->GetN(),v,br); }
    g->SetLineColor(sp[s].col); g->SetLineWidth(3); g->Draw("L same");
    l->SetTextColor(sp[s].col); l->DrawLatex(0.70,0.86-0.045*s,sp[s].n);
  }
  auto* gg=new TGraph(gx.size()+1);
  for(size_t i=0;i<gx.size();++i) gg->SetPoint(i,gx[i],gy[i]);
  gg->SetPoint(gx.size(),gx[0],gy[0]);
  gg->SetLineColor(kBlack); gg->SetLineWidth(4); gg->Draw("L same");
  l->SetTextColor(kBlack); l->DrawLatex(0.60,0.92,"black = saved deuteron gate");
  c->SaveAs("/home/yassid/C15d_gate_vs_bands.png");
  printf("wrote /home/yassid/C15d_gate_vs_bands.png\n");
}
