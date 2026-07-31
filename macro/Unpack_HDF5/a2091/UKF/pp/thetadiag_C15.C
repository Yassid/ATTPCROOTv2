/// @file thetadiag_C15.C
/// @brief Ex-versus-theta_LAB diagnostic for a2091, and the KE bias it implies.
///
/// Use theta_lab, NOT theta_cm. theta_cm is computed FROM Ex inside kine_2b (via m4_ex), so the
/// two are algebraically coupled: an event with a mismeasured KE gets a wrong Ex AND a wrong
/// theta_cm, correlated along dEx/dtheta_cm (about -0.07 MeV/deg for (p,p') at 77 deg, -0.11 for
/// (p,d) at 30 deg). Binning in theta_cm therefore partly sorts events by their own Ex and slides
/// the peak, producing a tilt even for perfect data. theta_lab is measured independently of Ex,
/// so a trend there is real.
///
/// For each theta_lab slice this fits the reference peak (g.s. / elastic, expected at Ex = 0) in a
/// window ANCHORED near zero -- never at the maximum bin, which in (p,d) jumps onto the 6-7 MeV
/// 14C cluster and in (p,p') onto threshold junk.
///
/// It then converts the measured offset into the KE bias that would explain it,
///     dKE = -mu / (dEx/dKE)  evaluated at that slice's (theta, <KE>),
/// which is the quantity to compare against an energy-loss/path-length model. Slices whose <KE>
/// is within keFlag of the 0.79 MeV reconstruction threshold are marked (!) and excluded from the
/// straight-line fit, because there the measured ridge is the threshold, not the physics.
///
///   root -b -q 'pp/thetadiag_C15.C("plots/proton_kin_g_ukf.root",170,1.007825,15.0105993)'
///   root -b -q 'pp/thetadiag_C15.C("plots/proton_kin_pd_ukf.root",170,2.014102,14.003242,10,52,3)'

static double om(double x, double y, double z) { return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z); }

static double exOf(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3;
   double s = m1*m1 + m2*m2 + 2*m2*E1, u = m2*m2 + m3*m3 - 2*m2*E3;
   double a = (cos(th)*om(s,m1*m1,m2*m2)*om(u,m2*m2,m3*m3) - (s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)
              + s + u - m2*m2;
   return a < 0 ? std::nan("") : sqrt(a) - m4;
}

static bool anchoredPeak(TH1F *h, double &mu, double &er, double &sg, double lo, double hi)
{
   if (h->GetEntries() < 100) return false;
   int b1 = h->FindBin(lo), b2 = h->FindBin(hi), bm = b1; double best = -1;
   for (int b = b1; b <= b2; b++) if (h->GetBinContent(b) > best) { best = h->GetBinContent(b); bm = b; }
   double x0 = h->GetBinCenter(bm);
   TF1 g("g", "gaus(0)+pol1(3)", x0 - 0.9, x0 + 0.9);
   g.SetParameters(best, x0, 0.3, 0, 0);
   g.SetParLimits(2, 0.03, 1.2); g.SetParLimits(1, x0 - 0.9, x0 + 0.9);
   if (h->Fit(&g, "QRN") != 0) return false;
   mu = g.GetParameter(1); er = g.GetParError(1); sg = g.GetParameter(2);
   return er > 0 && er < 0.4;
}

void thetadiag_C15(TString cache = "plots/proton_kin_g_ukf.root", double Ebeam = 170.0,
                   double mEjectAmu = 1.007825, double mResidAmu = 15.0105993, double thLo = 60,
                   double thHi = 92, double dth = 4, double chi2Cut = 5.0, double keThresh = 0.79,
                   double keFlag = 2.0, TString outTag = "")
{
   gStyle->SetOptStat(0);
   const double u = 931.49401, m1 = 15.0105993*u, m2 = 1.007825*u, m3 = mEjectAmu*u, m4 = mResidAmu*u;
   TString dir = gSystem->DirName(__FILE__);
   TString path = cache.BeginsWith("/") ? cache : dir + "/" + cache;

   TFile *f = TFile::Open(path);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", path.Data()); return; }
   TNtuple *n = (TNtuple *)f->Get("pk");
   if (!n) { printf("no `pk` in %s\n", path.Data()); return; }
   float ke, th, vz, tc, ex, c2;
   n->SetBranchAddress("ke",&ke); n->SetBranchAddress("theta",&th); n->SetBranchAddress("vertexz",&vz);
   n->SetBranchAddress("thcm",&tc); n->SetBranchAddress("ex",&ex); n->SetBranchAddress("chi2ndf",&c2);
   std::vector<float> vk, vt;
   for (Long64_t i = 0; i < n->GetEntries(); i++) {
      n->GetEntry(i);
      if (c2 > chi2Cut || ke <= 0) continue;
      vk.push_back(ke); vt.push_back(th);
   }
   printf("\n=== %s  (Ebeam=%.0f, m3=%.6f, m4=%.6f) ===\n", gSystem->BaseName(path), Ebeam, mEjectAmu, mResidAmu);
   printf("tracks: %zu\n", vk.size());
   printf("%-12s %8s %9s %8s %8s %9s %10s  %s\n",
          "theta_lab", "N", "mu_Ex", "err", "sigma", "<KE>", "dKE_impl", "");
   printf("%s\n", TString('-', 78).Data());

   std::vector<double> gt, gm, ge;   // slices used for the trend
   for (double t = thLo; t < thHi; t += dth) {
      TH1F h("h", "", 300, -8, 16), hk("hk", "", 400, 0, 45);
      for (size_t i = 0; i < vk.size(); i++) {
         if (vt[i] < t || vt[i] >= t + dth) continue;
         double e = exOf(m1, m2, m3, m4, Ebeam, vt[i]*TMath::DegToRad(), vk[i]);
         if (!std::isnan(e)) { h.Fill(e); hk.Fill(vk[i]); }
      }
      double mu, er, sg, keM = hk.GetMean();
      bool nearThr = (keM < keFlag * keThresh);
      if (!anchoredPeak(&h, mu, er, sg, -1.5, 2.0)) {
         printf("%-12s %8.0f %9s %8s %8s %9.2f %10s  %s\n", Form("%.0f-%.0f", t, t + dth),
                h.GetEntries(), "-", "-", "-", keM, "-", nearThr ? "(!) near threshold" : "no fit");
         continue;
      }
      // dEx/dKE at this slice, then the KE bias that would remove the offset
      double thr = (t + dth/2) * TMath::DegToRad();
      double d = (exOf(m1,m2,m3,m4,Ebeam,thr,keM+0.25) - exOf(m1,m2,m3,m4,Ebeam,thr,keM-0.25)) / 0.5;
      double dke = (fabs(d) > 1e-6) ? -mu/d : std::nan("");
      printf("%-12s %8.0f %9.3f %8.3f %8.3f %9.2f %10.3f  %s\n", Form("%.0f-%.0f", t, t + dth),
             h.GetEntries(), mu, er, sg, keM, dke, nearThr ? "(!) near threshold, excluded" : "");
      if (!nearThr) { gt.push_back(t + dth/2); gm.push_back(mu); ge.push_back(er); }
   }

   if (gt.size() >= 3) {
      TGraphErrors g(gt.size());
      for (size_t i = 0; i < gt.size(); i++) {
         g.SetPoint(i, gt[i], gm[i]);
         g.SetPointError(i, 0, sqrt(ge[i]*ge[i] + 0.03*0.03));   // + a 30 keV systematic floor
      }
      TF1 lin("lin", "pol1", gt.front() - 2, gt.back() + 2);
      g.Fit(&lin, "QRN");
      printf("%s\n", TString('-', 78).Data());
      printf("straight-line fit over %zu clean slices (theta %.0f-%.0f):\n", gt.size(), gt.front(), gt.back());
      printf("   slope     = %+.4f +- %.4f MeV/deg\n", lin.GetParameter(1), lin.GetParError(1));
      printf("   at mid    = %+.4f MeV\n", lin.Eval(0.5*(gt.front()+gt.back())));
      printf("   chi2/ndf  = %.2f\n", lin.GetNDF() > 0 ? lin.GetChisquare()/lin.GetNDF() : 0.0);
      double sig = fabs(lin.GetParameter(1)) / std::max(1e-9, lin.GetParError(1));
      printf("   -> slope is %.1f sigma from zero: %s\n", sig,
             sig < 3 ? "consistent with FLAT" : "a real theta_lab trend");
      TCanvas *c = new TCanvas("ctd", "theta diag", 900, 650);
      g.SetMarkerStyle(20);
      g.SetTitle(Form("%s: reference peak vs #theta_{lab}, E_{beam}=%.0f;#theta_{lab} [deg];E_{x} [MeV]",
                      gSystem->BaseName(path), Ebeam));
      g.Draw("AP"); lin.SetLineColor(kRed); lin.Draw("same");
      auto *z = new TLine(gt.front()-2, 0, gt.back()+2, 0); z->SetLineStyle(2); z->Draw();
      c->SaveAs(dir + "/plots/thetadiag" + outTag + ".png");
   } else {
      printf("too few clean slices for a trend fit\n");
   }
   f->Close();
}
