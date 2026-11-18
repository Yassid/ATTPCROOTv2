#include "binacc.h"
#include "binbuild.C"
// IS THE BACKGROUND PROTON OR DEUTERON?  The PID plane is ln(sqrt(dEdx)*Brho).  At equal Brho a
// proton and a deuteron carry equal momentum, so the deuteron is slower and loses MORE energy:
// deuterons sit ABOVE protons.  Left: the ungated plane, showing both loci.  Right: where each
// Ex window sits, against the elastic (certain deuterons) and the proton locus.
void pidfig(){
  gStyle->SetOptStat(0);
  TFile* fp=TFile::Open("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/pid/points_mfi_C15d.root");
  TTree* pt=(TTree*)fp->Get("pts");
  float sq,br; int rn,ev,tid;
  pt->SetBranchAddress("sqrtdedx",&sq); pt->SetBranchAddress("brho",&br);
  pt->SetBranchAddress("run",&rn); pt->SetBranchAddress("event",&ev); pt->SetBranchAddress("trackID",&tid);
  std::unordered_map<long long,float> M; M.reserve(pt->GetEntries()*1.4);
  TH1D* hall=new TH1D("hall",";ln(#sqrt{dE/dx} #upoint B#rho);tracks",300,1.3,3.9);
  hall->SetDirectory(nullptr);
  for(Long64_t i=0;i<pt->GetEntries();i++){ pt->GetEntry(i);
    if(sq<=0||br<=0) continue; double L=std::log(sq*br); hall->Fill(L);
    M[((long long)rn<<40)|((long long)ev<<10)|(tid&0x3FF)]=(float)L; }
  fp->Close();
  Cfg k=CFG[0]; k.vzHi=900; k.chi2=5.0;
  Pars A=parsOf(k); A.vzLo=300;
  TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd=0; float icf=0; int r2,e2,t2;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("run",&r2);t->SetBranchAddress("event",&e2);t->SetBranchAddress("trackID",&t2);
  const int NW=7; double wlo[NW]={-0.5,2,3,4,5,7,9}, whi[NW]={0.5,3,4,5,7,9,12};
  const char* wn[NW]={"elastic","2-3","3-4","4-5","5-7","7-9","9-12"};
  long N[NW]={0}; double S[NW]={0},S2[NW]={0};
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    if(c2>k.chi2||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double eb=ebAt(A,vz,true), thc=thCorr(A,thf,ke);
    double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
    double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>60) continue;
    auto it=M.find(((long long)r2<<40)|((long long)e2<<10)|(t2&0x3FF));
    if(it==M.end()) continue;
    for(int wq=0;wq<NW;wq++) if(e>=wlo[wq]&&e<whi[wq]){ N[wq]++; S[wq]+=it->second; S2[wq]+=it->second*it->second; } }
  f->Close();
  const double DEU=2.4888, PRO=1.9712;
  TCanvas* c=new TCanvas("c","c",1250,520); c->Divide(2,1);
  c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetLogy();
  hall->SetLineColor(kBlack); hall->SetLineWidth(2); hall->Draw("HIST");
  { auto* l1=new TLine(DEU,0.5,DEU,hall->GetBinContent(hall->GetMaximumBin())*1.3);
    l1->SetLineColor(kAzure+2); l1->SetLineWidth(3); l1->SetLineStyle(2); l1->Draw();
    auto* l2=new TLine(PRO,0.5,PRO,hall->GetBinContent(hall->GetMaximumBin())*1.3);
    l2->SetLineColor(kGreen+3); l2->SetLineWidth(3); l2->SetLineStyle(2); l2->Draw();
    auto* tx=new TLatex(); tx->SetTextSize(0.034);
    tx->SetTextColor(kAzure+2); tx->DrawLatex(DEU+0.05,hall->GetBinContent(hall->GetMaximumBin())*0.9,"deuteron");
    tx->SetTextColor(kGreen+3); tx->DrawLatex(PRO-0.55,hall->GetBinContent(hall->GetMaximumBin())*0.9,"proton"); }
  hall->SetTitle("the PID plane, ungated");
  c->cd(2); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetGridy();
  auto* fr=gPad->DrawFrame(-0.6,PRO-0.12,6.6,DEU+0.30);
  fr->SetTitle("where each E_{x} window sits;;ln(#sqrt{dE/dx} #upoint B#rho)");
  fr->GetXaxis()->SetNdivisions(0);
  { auto* b=new TBox(-0.6,PRO-0.12,6.6,PRO+0.12); b->SetFillColorAlpha(kGreen+3,0.15);
    b->SetLineColor(0); b->Draw();
    auto* l=new TLine(-0.6,PRO,6.6,PRO); l->SetLineColor(kGreen+3); l->SetLineWidth(3);
    l->SetLineStyle(2); l->Draw(); }
  auto* g=new TGraphErrors();
  for(int wq=0;wq<NW;wq++){ if(!N[wq]) continue;
    double m=S[wq]/N[wq], r=std::sqrt(std::max(0.0,S2[wq]/N[wq]-m*m));
    int p=g->GetN(); g->SetPoint(p,wq,m); g->SetPointError(p,0,r/std::sqrt((double)N[wq]));
    auto* tx=new TLatex(wq,PRO-0.06,wn[wq]); tx->SetTextAlign(22); tx->SetTextSize(0.030); tx->Draw(); }
  g->SetMarkerStyle(20); g->SetMarkerSize(1.5); g->SetMarkerColor(kAzure+2);
  g->SetLineColor(kAzure+2); g->SetLineWidth(3); g->Draw("PL");
  { auto* l=new TLine(-0.6,DEU,6.6,DEU); l->SetLineColor(kAzure+2); l->SetLineWidth(2);
    l->SetLineStyle(3); l->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.032);
    tx->SetTextColor(kGreen+3); tx->DrawLatex(0.22,0.17,"proton locus  (0.52 BELOW)");
    tx->SetTextColor(kGray+3); tx->SetTextSize(0.030);
    tx->DrawLatex(0.30,0.80,"#splitline{the background drifts UP, away from protons;}{rms unchanged from the elastic -- one species}"); }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/background_is_deuteron.png");
  printf("   wrote background_is_deuteron.png\n");
}
