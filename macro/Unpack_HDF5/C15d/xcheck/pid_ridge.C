// RIDGE FOLLOWER. Every earlier measurement was limited by per-slice peak finding, which mislabels
// bands (it called 12.93 a triton when the deuteron was at 12.56, and split a proton band in two).
// Follow each band CONTINUOUSLY in Brho instead: seed at a Brho where the band is unambiguous via
// Bethe, then step and search only near the previous centroid. Continuity is what makes it robust.
//
// Outputs per band, per Brho bin:  mu, sigma, sigma/mu, amplitude, local background, S/B.
#include <vector>
#include <cstdio>
double betheRel(double brho,double A,double z){
  double m=A*931.494, pm=299.792458*brho*z;
  double E=sqrt(pm*pm+m*m),b=pm/E,g=E/m,b2=b*b; if(b2<1e-9)b2=1e-9;
  const double I=19.2e-6,me=0.511; double a=2*me*b2*g*g/I; if(a<1.)a=1.;
  return z*z/b2*(log(a)-b2);
}
struct Pt { double brho, mu, sg, amp, bkg; bool ok; };

// fit one band in one Brho bin, searching near `guess`
Pt fitBin(TH2F* h2, int ib, double guess, double halfwin){
  Pt r{0,0,0,0,0,false};
  r.brho = h2->GetYaxis()->GetBinCenter(ib);
  TH1D* px = h2->ProjectionX("px_tmp", ib, ib);
  px->SetDirectory(nullptr);
  if(px->Integral() < 60){ delete px; return r; }
  double lo=guess*(1-halfwin), hi=guess*(1+halfwin);
  if(lo < px->GetXaxis()->GetXmin()) lo = px->GetXaxis()->GetXmin();
  if(hi > px->GetXaxis()->GetXmax()) hi = px->GetXaxis()->GetXmax();
  if(px->Integral(px->FindBin(lo), px->FindBin(hi)) < 40){ delete px; return r; }
  // local maximum inside the window = the ridge point
  int bbest=0; double vbest=-1;
  for(int b=px->FindBin(lo); b<=px->FindBin(hi); ++b)
    if(px->GetBinContent(b) > vbest){ vbest=px->GetBinContent(b); bbest=b; }
  double seed = px->GetBinCenter(bbest);
  TF1 f("rf","gaus(0)+pol1(3)", seed*(1-halfwin*0.9), seed*(1+halfwin*0.9));
  f.SetParameters(vbest, seed, 0.09*seed, 1, 0);
  f.SetParLimits(0, 0, 1e9);
  f.SetParLimits(1, seed*(1-halfwin*0.8), seed*(1+halfwin*0.8));
  f.SetParLimits(2, 0.015*seed, 0.30*seed);
  int st = px->Fit(&f, "RQ0N");
  if(st==0){
    r.mu=f.GetParameter(1); r.sg=fabs(f.GetParameter(2)); r.amp=f.GetParameter(0);
    r.bkg=std::max(0.0, f.GetParameter(3)+f.GetParameter(4)*r.mu);
    r.ok = r.sg>0 && r.amp>0 && r.mu>lo && r.mu<hi;
  }
  delete px; return r;
}

// RAIL-GUIDED tracking. Free-running continuity DRIFTS: seeded on the triton at 16.08 it slid
// down onto the deuteron at 9.5-11.2 and reported 24 happy bins of the wrong band. So do not
// search near the PREVIOUS fit -- search near the band's own Bethe locus, whose normalisation K is
// fixed once at the seed. The locus is then a rail the search cannot leave, and a lost bin is
// simply dropped instead of dragging the chain onto a neighbour.
std::vector<Pt> followRail(TH2F* h2, double brhoSeed, double muSeed, double A, double halfwin=0.13){
  std::vector<Pt> out;
  double K = muSeed / sqrt(betheRel(brhoSeed, A, 1));     // fixed normalisation
  for(int ib=1; ib<=h2->GetNbinsY(); ++ib){
    double br = h2->GetYaxis()->GetBinCenter(ib);
    if(br<=0.02) continue;
    double pred = K * sqrt(betheRel(br, A, 1));
    if(pred<1.0 || pred>44.0) continue;
    Pt p = fitBin(h2, ib, pred, halfwin);
    // keep only if the fit stayed on the rail
    if(p.ok && fabs(p.mu-pred)/pred < halfwin*0.85) out.push_back(p);
  }
  return out;
}

TH2F* loadPlane(const char* file, double icLo, double icHi, double sc, const char* nm){
  TFile* f = TFile::Open(file);
  if(!f || f->IsZombie()){ printf("CANNOT OPEN %s\n", file); return nullptr; }
  TTree* t=(TTree*)f->Get("pts"); if(!t){ printf("no pts in %s\n", file); return nullptr; }
  auto* h=new TH2F(nm,nm, 450,0,45, 120,0.15,1.35);
  h->SetDirectory(nullptr);
  float x,y,ic; int np;
  t->SetBranchAddress("sqrtdedx",&x); t->SetBranchAddress("brho",&y);
  t->SetBranchAddress("ic",&ic); t->SetBranchAddress("npulse",&np);
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    if(ic<icLo||ic>icHi||np!=1||x<=0||y<=0) continue;
    h->Fill(x*sc, y); }
  f->Close();
  return h;
}
