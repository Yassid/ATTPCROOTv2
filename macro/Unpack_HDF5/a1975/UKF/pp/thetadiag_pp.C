/// @file thetadiag_pp.C
/// @brief Ex-versus-theta_LAB diagnostic for the a1975 16C(p,p') elastic peak.
///
/// Port of a2091's pp/thetadiag_C15.C to 16C + p (see that file for the full rationale).
/// The one rule worth repeating: bin in theta_LAB, never theta_cm. theta_cm is computed FROM Ex
/// inside the two-body kinematics, so binning in it partly sorts events by their own Ex and
/// manufactures a tilt out of perfect data. theta_lab is measured independently of Ex, so a trend
/// there is a real instrumental effect.
///
/// Reads the flat per-run caches written by pp/cache_pp_run.C (branches ke, theta[deg], vz,
/// chi2ndf, ic, run) -- NOT the a2091 ntuple layout, so Ex is recomputed here from Ebeam.
///
/// The anchor window defaults to [-1.5, 1.2] MeV, deliberately BELOW the 16C 2+ at 1.766 MeV:
/// a wider window lets the fit slide onto the 2+ in slices where inelastic dominates, which would
/// masquerade as a tilt. Widen it only if you know the elastic dominates the slice.
///
///   root -b -q 'pp/thetadiag_pp.C("/tmp/.../ppcache/pp_kin_genfit.root")'
///   root -b -q 'pp/thetadiag_pp.C("...pp_kin_ukf.root",192,30,100,5,5,950,1350,"_ukf")'

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

/// Masses default to 16C(p,p')16C; for the D2-target elastic 16C(d,d)16C pass
/// mTarg = mEject = 2.0135532.
void thetadiag_pp(TString cache, double Ebeam = 192.0, double thLo = 30, double thHi = 100,
                  double dth = 5, double chi2Cut = 5.0, double icMin = 950, double icMax = 1350,
                  TString outTag = "", double anchorLo = -1.5, double anchorHi = 1.2,
                  double keThresh = 0.79, double keFlag = 2.0,
                  double mBeamAmu = 16.0147013, double mTargAmu = 1.00782503,
                  double mEjectAmu = 1.00782503, double mResidAmu = 16.0147013)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = mBeamAmu*u, m2 = mTargAmu*u, m3 = mEjectAmu*u, m4 = mResidAmu*u;
   TString dir = gSystem->DirName(__FILE__);

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *n = (TTree *)f->Get("pk");
   if (!n) { printf("no `pk` in %s\n", cache.Data()); return; }
   // (p,p') caches use `vz` + `ic` + `run`; the D2 (d,d) cache uses `vertexz` and has no IC
   float ke, th, vz = 0, c2, ic = -1;
   n->SetBranchAddress("ke",&ke); n->SetBranchAddress("theta",&th);
   n->SetBranchAddress("chi2ndf",&c2);
   if (n->GetBranch("vz")) n->SetBranchAddress("vz",&vz);
   else if (n->GetBranch("vertexz")) n->SetBranchAddress("vertexz",&vz);
   const bool hasIC = n->GetBranch("ic") != nullptr;
   if (hasIC) n->SetBranchAddress("ic",&ic);
   if (!hasIC && icMin > 0) { printf("note: no `ic` branch -- IC gate disabled\n"); icMin = -1; }

   std::vector<float> vk, vt;
   Long64_t nIC = 0;
   for (Long64_t i = 0; i < n->GetEntries(); i++) {
      n->GetEntry(i);
      if (c2 > chi2Cut || ke <= 0) continue;
      if (icMin > 0 && (ic < icMin || ic > icMax)) { ++nIC; continue; }
      vk.push_back(ke); vt.push_back(th);
   }
   printf("\n=== %s  16C(p,p'), Ebeam=%.1f MeV ===\n", gSystem->BaseName(cache), Ebeam);
   printf("entries %lld -> %zu after chi2/ndf<%.1f + IC[%.0f,%.0f] (IC rejected %lld)\n",
          n->GetEntries(), vk.size(), chi2Cut, icMin, icMax, nIC);
   printf("%-12s %8s %9s %8s %8s %9s %10s  %s\n",
          "theta_lab", "N", "mu_Ex", "err", "sigma", "<KE>", "dKE_impl", "");
   printf("%s\n", TString('-', 86).Data());

   std::vector<double> gt, gm, ge;
   for (double t = thLo; t < thHi; t += dth) {
      TH1F h("h", "", 300, -8, 16), hk("hk", "", 400, 0, 45);
      for (size_t i = 0; i < vk.size(); i++) {
         if (vt[i] < t || vt[i] >= t + dth) continue;
         double e = exOf(m1, m2, m3, m4, Ebeam, vt[i]*TMath::DegToRad(), vk[i]);
         if (!std::isnan(e)) { h.Fill(e); hk.Fill(vk[i]); }
      }
      double mu, er, sg, keM = hk.GetMean();
      // Near the reconstruction threshold the measured ridge IS the threshold, not the physics,
      // so those slices are reported but kept out of the trend fit (a2091 convention).
      bool nearThr = (keM < keFlag * keThresh);
      if (!anchoredPeak(&h, mu, er, sg, anchorLo, anchorHi)) {
         printf("%-12s %8.0f %9s %8s %8s %9.2f %10s  %s\n", Form("%.0f-%.0f", t, t + dth),
                h.GetEntries(), "-", "-", "-", keM, "-", nearThr ? "(!) near threshold" : "no fit");
         continue;
      }
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
      printf("%s\n", TString('-', 70).Data());
      printf("straight-line fit over %zu slices (theta %.0f-%.0f):\n", gt.size(), gt.front(), gt.back());
      printf("   slope    = %+.4f +- %.4f MeV/deg\n", lin.GetParameter(1), lin.GetParError(1));
      printf("   at mid   = %+.4f MeV\n", lin.Eval(0.5*(gt.front()+gt.back())));
      printf("   chi2/ndf = %.2f\n", lin.GetNDF() > 0 ? lin.GetChisquare()/lin.GetNDF() : 0.0);
      double sig = fabs(lin.GetParameter(1)) / std::max(1e-9, lin.GetParError(1));
      printf("   -> slope is %.1f sigma from zero: %s\n", sig,
             sig < 3 ? "consistent with FLAT" : "a real theta_lab trend");
      TCanvas *c = new TCanvas("ctd", "theta diag", 900, 650);
      g.SetMarkerStyle(20);
      g.SetTitle(Form("16C(p,p') elastic peak vs #theta_{lab}, E_{beam}=%.0f;#theta_{lab} [deg];E_{x} [MeV]", Ebeam));
      g.Draw("AP"); lin.SetLineColor(kRed); lin.Draw("same");
      auto *z = new TLine(gt.front()-2, 0, gt.back()+2, 0); z->SetLineStyle(2); z->Draw();
      gSystem->Exec("mkdir -p " + dir + "/plots");
      c->SaveAs(dir + "/plots/thetadiag_pp" + outTag + ".png");
   } else {
      printf("too few slices for a trend fit\n");
   }
   f->Close();
}
