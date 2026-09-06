/// @file arcz_vertex_dp.C
/// @brief Yassid's proposal: use the (arc length, z) LINE FIT to give both the direction and the
///        vertex, from where the extrapolated line reaches the track's closest approach to the
///        beam axis.
///
/// This is the same machinery AtPRA::SetTrackInitialParameters already runs, but with two
/// differences that matter:
///   (a) the azimuth is unwrapped along the CLUSTER array, which OrderClustersAlongTrack has
///       already put in nearest-neighbour order. AtPRA unwraps along the HIT array, which nothing
///       orders -- that is the bug that makes its theta unreliable here.
///   (b) the fitted line is then EXTRAPOLATED to arc*, the arc position where the circle passes
///       closest to the beam axis, giving a predicted vertex (x*, y*, z*). My earlier
///       vertex_end_dp.C used only the circle angles and scored 79%; the z brings a genuine extra
///       constraint, because a vertex predicted OUTSIDE the active volume excludes that direction.
///
/// Scored against the trusted label (fit both ways, keep the vertex nearest the beam axis).
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

static double azWrap(double d){ while(d>M_PI) d-=2*M_PI; while(d<=-M_PI) d+=2*M_PI; return d; }

static bool azCircle(const std::vector<ROOT::Math::XYZPoint> &P, double &cx, double &cy, double &R)
{
   const int m=P.size(); if(m<5) return false;
   double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
   for(auto&p:P){double x=p.X(),y=p.Y();
      Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y;}
   double C=m*Sxx-Sx*Sx,D=m*Sxy-Sx*Sy,E=m*Sxxx+m*Sxyy-(Sxx+Syy)*Sx;
   double G=m*Syy-Sy*Sy,H=m*Sxxy+m*Syyy-(Sxx+Syy)*Sy;
   double den=2*(C*G-D*D); if(std::fabs(den)<1e-9) return false;
   cx=(E*G-D*H)/den; cy=(C*H-D*E)/den; R=0;
   for(auto&p:P) R+=std::hypot(p.X()-cx,p.Y()-cy);
   R/=m; return R>1.0;
}

void arcz_vertex_dp(TString listFile="/mnt/f/a1975/caches/.explorer/refit_list.txt",
                    TString csvOut="/tmp/arcz.csv",
                    TString gfDir="/mnt/f/a1975/gf_dp_cateloss/",
                    double zLo=0.0, double zHi=1000.0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while(in>>r>>e>>ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit",e,ti});
   std::ofstream o(csvOut.Data());
   o << "run,entry,tid,ncl,ok,reverse,slope,zvtxF,zvtxB,dArcF,dArcB,inVolF,inVolB,dstar,resid\n";
   TString open; TFile *f=nullptr; TTree *t=nullptr; TClonesArray *te=nullptr; long n=0;
   for(auto &j:jobs){
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if(rt!=open){ if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); te=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&te); open=rt; }
      t->GetEntry(en);
      auto *ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      auto trks=ev->GetTrackArray(); AtTrack *src=nullptr;
      for(auto &tr:trks) if(tr.GetTrackID()==id) src=&tr;
      if(!src) continue;
      auto *hc=src->GetHitClusterArray(); if(!hc||hc->size()<8) continue;
      std::vector<ROOT::Math::XYZPoint> P;
      for(auto &c:*hc) P.push_back(c.GetPosition());
      double cx,cy,R;
      if(!azCircle(P,cx,cy,R)){ ++n; continue; }
      const double cR=std::hypot(cx,cy);
      if(cR<1e-6){ ++n; continue; }
      const double dstar=std::fabs(cR-R);
      // unwrap along the CLUSTER order (already NN-ordered by OrderClustersAlongTrack)
      std::vector<double> u(P.size());
      u[0]=std::atan2(P[0].Y()-cy,P[0].X()-cx);
      for(size_t i=1;i<P.size();++i){
         double a=std::atan2(P[i].Y()-cy,P[i].X()-cx);
         u[i]=u[i-1]+azWrap(a-u[i-1]);
      }
      // (arc, z) least-squares line:  z = A + B*arc,   arc = R*u
      const int m=P.size();
      double Sx=0,Sy=0,Sxx=0,Sxy=0;
      for(int i=0;i<m;++i){ double x=R*u[i], y=P[i].Z(); Sx+=x;Sy+=y;Sxx+=x*x;Sxy+=x*y; }
      double den=m*Sxx-Sx*Sx;
      if(std::fabs(den)<1e-9){ ++n; continue; }
      const double B=(m*Sxy-Sx*Sy)/den, A=(Sy-B*Sx)/m;
      double resid=0;
      for(int i=0;i<m;++i){ double d=P[i].Z()-(A+B*R*u[i]); resid+=d*d; }
      resid=std::sqrt(resid/m);
      // arc* : azimuth of the circle's closest approach to the beam axis
      const double aStar=std::atan2(cy*(cR-R)/cR-cy, cx*(cR-R)/cR-cx);
      // nearest branch of aStar OUTWARD from each end
      auto outward=[&](double from,int dir){ double d=azWrap(aStar-from)*dir; if(d<0) d+=2*M_PI; return d; };
      const int sense=(u.back()>=u.front())?+1:-1;
      const double dF=outward(u.front(),-sense), dB=outward(u.back(),+sense);
      // z of the predicted vertex, from the LINE, at each candidate end
      const double zF=A+B*R*(u.front()-sense*dF);
      const double zB=A+B*R*(u.back() +sense*dB);
      const bool inF=(zF>zLo&&zF<zHi), inB=(zB>zLo&&zB<zHi);
      // direction: prefer the end whose extrapolated vertex is INSIDE the volume; if both or
      // neither, fall back to the shorter extrapolation
      bool rev;
      if(inF!=inB) rev=inB;
      else         rev=(dB<dF);
      o<<rt<<","<<en<<","<<id<<","<<P.size()<<",1,"<<(rev?1:0)<<","<<B<<","<<zF<<","<<zB<<","
       <<dF<<","<<dB<<","<<(inF?1:0)<<","<<(inB?1:0)<<","<<dstar<<","<<resid<<"\n";
      ++n;
   }
   if(f) f->Close(); o.close();
   printf("  wrote %s  (%ld tracks)\n", csvOut.Data(), n);
}
