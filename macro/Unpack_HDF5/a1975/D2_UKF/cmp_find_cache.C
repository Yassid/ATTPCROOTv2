/// Baseline (d,p) cache vs the FIND cache. Population level -- the pk tree carries no entry/tid,
/// so track-by-track matching has to be done at the genfit level with cmp_gf_dirs.C.
///
/// The referee is NOT the raw Ex IQR: that is dominated by the level structure plus continuum, so
/// it cannot see a resolution change. What CAN see a curvature-averaging bias is the drift of Ex
/// with theta_cm -- Ex must not depend on angle, and track length (hence spiral tightness)
/// correlates strongly with angle. A bias that grows with spiralling shows up as a slope.
#include <vector>
static void exVsTheta(TTree *t, const char *cut, double &slope, double &err, long &n)
{
   slope = err = 0; n = t->Draw("ex:thcm", cut, "goff");
   if (n < 20) return;
   double *y = t->GetV1(), *x = t->GetV2();
   double Sx=0,Sy=0,Sxx=0,Sxy=0;
   for (long i=0;i<n;++i){ Sx+=x[i]; Sy+=y[i]; Sxx+=x[i]*x[i]; Sxy+=x[i]*y[i]; }
   double den = n*Sxx - Sx*Sx; if (std::fabs(den)<1e-9) return;
   slope = (n*Sxy - Sx*Sy)/den;
   double a = (Sy - slope*Sx)/n, s2 = 0;
   for (long i=0;i<n;++i){ double r=y[i]-(a+slope*x[i]); s2+=r*r; }
   err = std::sqrt(s2/(n-2) * n/den);
}
static void stat(TTree *t, const char *var, const char *cut, double &med, double &iqr, long &n)
{
   med=iqr=0; n=t->Draw(var,cut,"goff"); if(n<5) return;
   std::vector<double> v(t->GetV1(), t->GetV1()+n); std::sort(v.begin(),v.end());
   med=v[v.size()/2]; iqr=v[(size_t)(0.75*v.size())]-v[(size_t)(0.25*v.size())];
}
void cmp_find_cache(TString a="/mnt/f/a1975/caches/dp_kin_dv1104.root",
                    TString b="/mnt/f/a1975/caches/dp_kin_find.root")
{
   TFile *fa=TFile::Open(a), *fb=TFile::Open(b);
   TTree *ta=(TTree*)fa->Get("pk"), *tb=(TTree*)fb->Get("pk");
   printf("\n  candidates      baseline %8lld    FIND %8lld    (%+.1f%%)\n",
          ta->GetEntries(), tb->GetEntries(),
          100.0*(tb->GetEntries()-ta->GetEntries())/ta->GetEntries());
   const char *cuts[3]={"", "theta>90", "chi2ndf<5&&theta>90"};
   const char *lbl[3]={"all", "backward", "backward, chi2<5"};
   for (int c=0;c<3;++c){
      double m1,q1,m2,q2; long n1,n2;
      stat(ta,"ex",cuts[c],m1,q1,n1); stat(tb,"ex",cuts[c],m2,q2,n2);
      printf("\n  --- %s ---\n", lbl[c]);
      printf("    n            %8ld   %8ld   (%+.1f%%)\n", n1,n2, n1?100.0*(n2-n1)/n1:0.0);
      printf("    median Ex    %8.3f   %8.3f   (%+.3f MeV)\n", m1,m2,m2-m1);
      printf("    IQR Ex       %8.3f   %8.3f   (%+.3f)   [weak metric -- level structure]\n", q1,q2,q2-q1);
      double s1,e1,s2,e2; long k1,k2;
      exVsTheta(ta,cuts[c],s1,e1,k1); exVsTheta(tb,cuts[c],s2,e2,k2);
      printf("    dEx/dtheta_cm %7.4f+-%.4f  %7.4f+-%.4f   MeV/deg   <== THE REFEREE\n", s1,e1,s2,e2);
   }
   double m1,q1,m2,q2; long n1,n2;
   stat(ta,"chi2ndf","",m1,q1,n1); stat(tb,"chi2ndf","",m2,q2,n2);
   printf("\n  median chi2/ndf  %8.3f   %8.3f\n", m1, m2);
   stat(ta,"chi2ndf","theta>90",m1,q1,n1); stat(tb,"chi2ndf","theta>90",m2,q2,n2);
   printf("  median chi2/ndf backward %.3f   %.3f\n", m1, m2);
   printf("  passing chi2<5   %.1f%%      %.1f%%\n",
          100.0*ta->GetEntries("chi2ndf<5")/ta->GetEntries(),
          100.0*tb->GetEntries("chi2ndf<5")/tb->GetEntries());
}
