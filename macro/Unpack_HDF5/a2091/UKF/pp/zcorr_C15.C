/// @file zcorr_C15.C
/// @brief Vertex-z tilt correction for the 15C(p,d)14C ground state, and its resolution gain.
///
/// Why z and not theta_cm, and why this angular region:
///   * theta_cm is computed FROM Ex inside kine_2b, so subtracting a function of it subtracts a
///     function of each event's own Ex -- that compresses the spectrum rather than sharpening it.
///   * vertex z is measured independently of Ex, so a trend in z is a real instrumental effect
///     (path length / energy loss) and removing it is a genuine resolution gain.
///   * the g.s. is only cleanly populated at theta_lab ~ 40-48 deg (at 16-32 deg the 6-7 MeV 14C
///     cluster dominates and the g.s. has <100 counts per bin; at 32-40 deg ~80% of events have
///     the unphysical Ex < 0). Restricting to 40-48 keeps ONE state, so no multi-anchor stitching
///     and no assumption about the 4-level cluster centroid is needed.
///
/// Reports sigma before and after so the gain is explicit, and the surviving theta_lab trend so it
/// is clear how much of the tilt was really a z effect.
///
///   root -b -q 'pp/zcorr_C15.C()'

static double om(double x,double y,double z){ return sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
static double exOf(double m1,double m2,double m3,double m4,double K,double th,double Ke)
{
   double E1=K+m1,E3=Ke+m3,s=m1*m1+m2*m2+2*m2*E1,u=m2*m2+m3*m3-2*m2*E3;
   double a=(cos(th)*om(s,m1*m1,m2*m2)*om(u,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)+s+u-m2*m2;
   return a<0?std::nan(""):sqrt(a)-m4;
}
static bool pk(TH1F*h,double c,double win,double&mu,double&er,double&sg)
{
   if(h->GetEntries()<80) return false;
   int bm=h->GetMaximumBin();
   TF1 g("g","gaus(0)+pol1(3)",c-win,c+win);
   g.SetParameters(h->GetBinContent(bm),c,0.35,0,0);
   g.SetParLimits(1,c-win,c+win); g.SetParLimits(2,0.03,1.5);
   if(h->Fit(&g,"QRN")!=0) return false;
   mu=g.GetParameter(1); er=g.GetParError(1); sg=g.GetParameter(2);
   return er>0 && er<0.5;
}

void zcorr_C15(TString cache="plots/proton_kin_pd_ukf.root", double Ebeam=170.0,
               double mEjectAmu=2.014102, double mResidAmu=14.003242,
               double thLo=40, double thHi=48, double exWin=2.0)
{
   gStyle->SetOptStat(0);
   const double u=931.49401,m1=15.0105993*u,m2=1.007825*u,m3=mEjectAmu*u,m4=mResidAmu*u;
   TString dir=gSystem->DirName(__FILE__);
   TFile*f=TFile::Open(cache.BeginsWith("/")?cache:dir+"/"+cache);
   TNtuple*n=(TNtuple*)f->Get("pk");
   float ke,th,vz,tc,ex,c2;
   n->SetBranchAddress("ke",&ke);n->SetBranchAddress("theta",&th);n->SetBranchAddress("vertexz",&vz);
   n->SetBranchAddress("thcm",&tc);n->SetBranchAddress("ex",&ex);n->SetBranchAddress("chi2ndf",&c2);

   struct Ev{double th,vz,ex;};
   std::vector<Ev> ev;
   for(Long64_t i=0;i<n->GetEntries();i++){
      n->GetEntry(i);
      if(c2>5||ke<=0||th<thLo||th>=thHi) continue;
      double e=exOf(m1,m2,m3,m4,Ebeam,th*TMath::DegToRad(),ke);
      if(std::isnan(e)||fabs(e)>exWin) continue;      // the g.s. band only
      ev.push_back({(double)th,(double)vz,e});
   }
   printf("\n=== g.s. z-correction, theta_lab %.0f-%.0f, |Ex|<%.1f ===\n",thLo,thHi,exWin);
   printf("g.s. events: %zu\n",ev.size());
   if(ev.size()<500){ printf("too few\n"); return; }

   // ---- 1. Ex vs vertex z ----
   printf("\n%-14s %8s %9s %8s %8s\n","z [mm]","N","mu","err","sigma");
   std::vector<double> zc,zm,ze;
   for(double z=0;z<1050;z+=100){
      TH1F h("h","",200,-exWin,exWin);
      for(auto&e:ev) if(e.vz>=z&&e.vz<z+100) h.Fill(e.ex);
      double mu,er,sg;
      if(pk(&h,0.0,1.4,mu,er,sg)){
         printf("%-14s %8.0f %9.3f %8.3f %8.3f\n",Form("%.0f-%.0f",z,z+100),h.GetEntries(),mu,er,sg);
         zc.push_back(z+50); zm.push_back(mu); ze.push_back(er);
      } else printf("%-14s %8.0f      (no fit)\n",Form("%.0f-%.0f",z,z+100),h.GetEntries());
   }
   if(zc.size()<4){ printf("too few z slices\n"); return; }

   TGraphErrors gz(zc.size());
   for(size_t i=0;i<zc.size();i++){ gz.SetPoint(i,zc[i],zm[i]); gz.SetPointError(i,0,sqrt(ze[i]*ze[i]+0.03*0.03)); }
   TF1 lz("lz","pol1",0,1050);
   gz.Fit(&lz,"QRN");
   printf("  Delta(z) = %+.4f %+.5f*z   MeV (z in mm)\n",lz.GetParameter(0),lz.GetParameter(1));
   printf("  slope %+.5f +- %.5f MeV/mm  -> %.3f MeV over 600 mm ; chi2/ndf %.2f\n",
          lz.GetParameter(1),lz.GetParError(1),600*lz.GetParameter(1),
          lz.GetNDF()>0?lz.GetChisquare()/lz.GetNDF():0.);

   // ---- 2. before / after ----
   printf("\n-- projected g.s. peak, before and after subtracting Delta(z) --\n");
   printf("%-26s %9s %9s %9s %9s\n","","mu","err","sigma","FWHM");
   TH1F h0("h0","",300,-exWin,exWin), h1("h1","",300,-exWin,exWin);
   for(auto&e:ev){ h0.Fill(e.ex); h1.Fill(e.ex-lz.Eval(e.vz)); }
   double m0,e0,s0,mm1,e1,s1;
   bool o0=pk(&h0,0.0,1.4,m0,e0,s0), o1=pk(&h1,0.0,1.4,mm1,e1,s1);
   if(o0) printf("%-26s %9.3f %9.3f %9.3f %9.3f\n","no correction",m0,e0,s0,2.355*s0);
   if(o1) printf("%-26s %9.3f %9.3f %9.3f %9.3f\n","minus Delta(z)",mm1,e1,s1,2.355*s1);
   if(o0&&o1) printf("  --> sigma %.3f -> %.3f MeV  (%+.1f%%)\n",s0,s1,100*(s1-s0)/s0);

   // ---- 3. what theta trend survives? ----
   printf("\n-- residual Ex vs theta_lab inside this region, after the z correction --\n");
   printf("%-12s %8s %9s %9s %9s\n","theta_lab","N","mu_raw","mu_zcorr","sigma_zc");
   for(double t=thLo;t<thHi;t+=2){
      TH1F a("a","",200,-exWin,exWin), b("b","",200,-exWin,exWin);
      for(auto&e:ev) if(e.th>=t&&e.th<t+2){ a.Fill(e.ex); b.Fill(e.ex-lz.Eval(e.vz)); }
      double ma,ea,sa,mb,eb,sb;
      bool oa=pk(&a,0.0,1.4,ma,ea,sa), ob=pk(&b,0.0,1.4,mb,eb,sb);
      printf("%-12s %8.0f %9s %9s %9s\n",Form("%.0f-%.0f",t,t+2),a.GetEntries(),
             oa?Form("%+.3f",ma):"-", ob?Form("%+.3f",mb):"-", ob?Form("%.3f",sb):"-");
   }

   TCanvas*c=new TCanvas("cz","zcorr",1500,520); c->Divide(3,1);
   c->cd(1); gz.SetMarkerStyle(20);
   gz.SetTitle("g.s. E_{x} vs vertex z;z [mm];E_{x} [MeV]"); gz.Draw("AP"); lz.SetLineColor(kRed); lz.Draw("same");
   c->cd(2); h0.SetTitle("before;E_{x} [MeV];counts"); h0.Draw();
   c->cd(3); h1.SetTitle("after #minus #Delta(z);E_{x} [MeV];counts"); h1.Draw();
   c->SaveAs(dir+"/plots/zcorr_pd.png");
   f->Close();
}
