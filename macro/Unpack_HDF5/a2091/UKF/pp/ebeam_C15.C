/// @file ebeam_C15.C
/// @brief Calibrate the a2091 15C beam energy from a proton_kin ntuple, two independent ways.
///
/// Method 1 -- ELASTIC RECOIL RIDGE. For 15C(p,p') the recoil proton energy at fixed theta is
///   T_p = 2 m_p p1^2 cos^2(th) / [ (E1+m_p)^2 - p1^2 cos^2(th) ],  E1 = T1+m1, p1^2 = T1^2+2T1m1
/// Fitting the measured ridge over the theta acceptance constrains T1 through the whole cos^2
/// shape. The ridge is WALKED from high theta (unambiguous, high statistics) down, searching only
/// near where the previous points extrapolate -- taking the global max per theta slice instead
/// jumps onto a second population below ~58 deg and gives nonsense.
///
/// Method 2 -- THETA TILT OF THE ELASTIC PEAK. A wrong assumed Ebeam tilts Ex versus theta, so the
/// correct Ebeam makes the elastic peak sit at the same place, and at zero, in every theta bin.
/// Ex is recomputed from the stored (ke, theta), so no refitting is needed for any trial Ebeam.
///
/// >>> Both methods must AVOID two acceptance edges found on the ungated data:
///     (a) a low-energy reconstruction threshold near KE = 0.79 MeV: above theta ~ 80 deg the
///         measured ridge flattens there while the elastic locus keeps falling, so those bins
///         are meaningless. Hence thHi defaults to 79.
///     (b) below theta ~ 62-66 deg a second population overlaps the ridge (and statistics fall).
///     If the IC + PID gate has removed (b), widen thLo and the answer should not move -- that
///     stability is itself the check.
///
///   root -b -q 'pp/ebeam_C15.C("plots/proton_kin_gated.root")'

static double omega2(double x, double y, double z)
{ return std::sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z); }

static double exOf(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double Et1 = K + m1, Et3 = Ke + m3;
   double s = m1*m1 + m2*m2 + 2*m2*Et1;
   double u = m2*m2 + m3*m3 - 2*m2*Et3;
   double a = (std::cos(th)*omega2(s,m1*m1,m2*m2)*omega2(u,m2*m2,m3*m3)
               - (s-m1*m1-m2*m2)*(m2*m2+m3*m3-u)) / (2*m2*m2) + s + u - m2*m2;
   return a < 0 ? std::nan("") : std::sqrt(a) - m4;
}

static double Trec(double T1, double m1, double m2, double th)
{
   double E1 = T1 + m1, p = T1*T1 + 2*T1*m1, c = std::cos(th)*std::cos(th);
   return 2*m2*p*c / ((E1+m2)*(E1+m2) - p*c);
}

static bool pkFit(TH1 *h, double x0, double win, double &mu, double &er, double &sg)
{
   if (h->GetEntries() < 100) return false;
   TF1 f("f", "gaus(0)+pol1(3)", x0-win, x0+win);
   f.SetParameters(h->GetMaximum(), x0, 0.3*win, 0, 0);
   f.SetParLimits(2, 0.02, 1.5*win);
   f.SetParLimits(1, x0-win, x0+win);
   if (h->Fit(&f, "QRN") != 0) return false;
   mu = f.GetParameter(1); er = f.GetParError(1); sg = f.GetParameter(2);
   return er > 0;
}

void ebeam_C15(TString cache = "plots/proton_kin_gated.root", double chi2Cut = 5.0, double thLo = 66,
               double thHi = 79, double mEjectAmu = 1.007825, double mResidAmu = 15.0105993)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = 15.0105993*u, m2 = 1.007825*u, m3 = mEjectAmu*u, m4 = mResidAmu*u;

   TString dir = gSystem->DirName(__FILE__);
   TString path = cache.BeginsWith("/") ? cache : dir + "/" + cache;
   TFile *f = TFile::Open(path);
   if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", path.Data()); return; }
   TNtuple *n = (TNtuple *)f->Get("pk");
   if (!n) { printf("ERROR: no `pk` ntuple in %s\n", path.Data()); return; }

   float ke, theta, vz, thcm, ex, c2;
   n->SetBranchAddress("ke",&ke); n->SetBranchAddress("theta",&theta); n->SetBranchAddress("vertexz",&vz);
   n->SetBranchAddress("thcm",&thcm); n->SetBranchAddress("ex",&ex); n->SetBranchAddress("chi2ndf",&c2);

   std::vector<float> vk, vt;
   for (Long64_t i = 0; i < n->GetEntries(); i++) {
      n->GetEntry(i);
      if (c2 > chi2Cut || ke <= 0 || ke > 25) continue;
      vk.push_back(ke); vt.push_back(theta);
   }
   printf("\n=== ebeam_C15 on %s ===\n", path.Data());
   printf("tracks: %lld total, %zu usable\n", n->GetEntries(), vk.size());
   if (vk.size() < 500) { printf("too few tracks for a calibration\n"); return; }

   // ---------------- Method 1: walked ridge ----------------
   const double dt = 1.5, wHi = 88.0, wLo = 54.0;
   const int NB = (int)std::round((wHi-wLo)/dt);
   std::vector<TH1F*> hk(NB);
   for (int b = 0; b < NB; b++) hk[b] = new TH1F(Form("hk%d",b), "", 400, 0, 25);
   for (size_t i = 0; i < vk.size(); i++) {
      int b = (int)((vt[i]-wLo)/dt);
      if (b >= 0 && b < NB) hk[b]->Fill(vk[i]);
   }
   printf("\n-- ridge walk (high theta -> low) --\n%-8s %9s %9s %9s\n","theta","N","KE_ridge","err");
   std::vector<double> rt, rk, re;
   double prev = -1, prev2 = -1;
   for (int b = NB-1; b >= 0; b--) {
      double tc = wLo + dt*(b+0.5);
      if (hk[b]->GetEntries() < 100) continue;
      double guess, lo, hi;
      if (prev > 0 && prev2 > 0) { guess = prev + (prev-prev2); double w = std::max(0.8, 0.35*guess); lo = guess-w; hi = guess+w; }
      else if (prev > 0)         { guess = prev; lo = 0.5*prev; hi = 2.0*prev; }
      else                       { guess = hk[b]->GetBinCenter(hk[b]->GetMaximumBin()); lo = 0; hi = 25; }
      if (lo < 0.2) lo = 0.2;
      int b1 = hk[b]->FindBin(lo), b2 = hk[b]->FindBin(hi), bm = b1; double best = -1;
      for (int k = b1; k <= b2; k++) if (hk[b]->GetBinContent(k) > best) { best = hk[b]->GetBinContent(k); bm = k; }
      double x0 = hk[b]->GetBinCenter(bm), mu, er, sg;
      if (!pkFit(hk[b], x0, std::max(0.5, 0.20*x0), mu, er, sg)) continue;
      if (mu <= 0.2 || er > 1.0) continue;
      printf("%-8.1f %9.0f %9.3f %9.3f\n", tc, hk[b]->GetEntries(), mu, er);
      rt.push_back(tc); rk.push_back(mu); re.push_back(er);
      prev2 = prev; prev = mu;
   }

   double E_ridge = 0, E_ridge_err = 0;
   printf("\n-- ridge fits over theta windows (systematic = spread across windows) --\n");
   printf("%-16s %5s %9s %8s %9s\n","window","npts","Ebeam","+-","chi2/ndf");
   {
      double lo[4] = {thLo, thLo, thLo+4, thLo-4}, hi[4] = {thHi, thHi-3, thHi, thHi};
      std::vector<double> got;
      for (int w = 0; w < 4; w++) {
         TGraphErrors g; int k = 0;
         for (size_t i = 0; i < rt.size(); i++) {
            if (rt[i] < lo[w] || rt[i] > hi[w]) continue;
            g.SetPoint(k, rt[i], rk[i]);
            g.SetPointError(k, 0, std::sqrt(re[i]*re[i] + std::pow(0.03*rk[i],2) + 0.0025));
            k++;
         }
         if (k < 4) { printf("%-16s %5d   (too few points)\n", Form("%.0f-%.0f",lo[w],hi[w]), k); continue; }
         TF1 fit("fit", [&](double *x, double *p){ return Trec(p[0], m1, m2, x[0]*TMath::DegToRad()); }, lo[w]-1, hi[w]+1, 1);
         fit.SetParameter(0, 160.0);
         g.Fit(&fit, "QRN");
         printf("%-16s %5d %9.2f %8.2f %9.2f\n", Form("%.0f-%.0f",lo[w],hi[w]), k, fit.GetParameter(0),
                fit.GetParError(0), fit.GetNDF()>0 ? fit.GetChisquare()/fit.GetNDF() : 0.0);
         got.push_back(fit.GetParameter(0));
      }
      if (!got.empty()) {
         double s = 0; for (double v : got) s += v; E_ridge = s/got.size();
         double mn = *std::min_element(got.begin(),got.end()), mx = *std::max_element(got.begin(),got.end());
         E_ridge_err = 0.5*(mx-mn);
         printf("  RIDGE: Ebeam = %.1f +- %.1f MeV (%.2f MeV/u)\n", E_ridge, E_ridge_err, E_ridge/15.0);
      }
   }

   // ---------------- Method 2: theta tilt ----------------
   printf("\n-- theta tilt scan (Ex must be theta-independent AND zero) --\n");
   const int TB = 5;
   double tlo[TB], thi[TB];
   for (int b = 0; b < TB; b++) { tlo[b] = thLo + b*(thHi-thLo)/TB; thi[b] = thLo + (b+1)*(thHi-thLo)/TB; }
   printf("%-7s %8s %8s   per-bin mu\n","Ebeam","tilt","|mean|");
   double bestE = 0, bestScore = 1e9, bestTilt = 0, bestMean = 0;
   for (int i = 0; i <= 50; i++) {
      double E = 120 + 2.0*i;
      TH1F *hb[TB];
      for (int b = 0; b < TB; b++) hb[b] = new TH1F(Form("tb%d",b), "", 400, -8, 12);
      for (size_t k = 0; k < vk.size(); k++) {
         if (vt[k] < tlo[0] || vt[k] >= thi[TB-1]) continue;
         double e = exOf(m1,m2,m3,m4,E, vt[k]*TMath::DegToRad(), vk[k]);
         if (std::isnan(e)) continue;
         for (int b = 0; b < TB; b++) if (vt[k] >= tlo[b] && vt[k] < thi[b]) { hb[b]->Fill(e); break; }
      }
      int nok = 0; double sum = 0, sum2 = 0; TString det;
      for (int b = 0; b < TB; b++) {
         int bm = hb[b]->GetMaximumBin();
         double x0 = hb[b]->GetBinCenter(bm), mu, er, sg;
         if (std::fabs(x0) < 3.0 && pkFit(hb[b], x0, 0.9, mu, er, sg) && er < 0.2) {
            sum += mu; sum2 += mu*mu; nok++; det += Form(" %6.3f", mu);
         } else det += "    ---";
         delete hb[b];
      }
      if (nok < TB) continue;
      double mean = sum/nok, tilt = std::sqrt(std::max(0.0, sum2/nok - mean*mean));
      double score = std::sqrt(tilt*tilt + mean*mean);
      if (fmod(E,10) < 1e-6 || score < bestScore) printf("%-7.0f %8.4f %8.4f   |%s\n", E, tilt, std::fabs(mean), det.Data());
      if (score < bestScore) { bestScore = score; bestE = E; bestTilt = tilt; bestMean = mean; }
   }
   printf("  TILT: Ebeam = %.0f MeV (%.2f MeV/u), residual tilt %.4f MeV, offset %+.4f MeV\n",
          bestE, bestE/15.0, bestTilt, bestMean);

   // ---------------- consensus ----------------
   printf("\n================= RESULT =================\n");
   printf("  ridge : %.1f +- %.1f MeV\n", E_ridge, E_ridge_err);
   printf("  tilt  : %.0f MeV\n", bestE);
   if (E_ridge > 0) {
      double cons = 0.5*(E_ridge + bestE);
      printf("  CONSENSUS Ebeam = %.0f MeV  (%.2f MeV/u); methods differ by %.1f MeV\n",
             cons, cons/15.0, std::fabs(E_ridge - bestE));
      if (std::fabs(E_ridge - bestE) > 15)
         printf("  \033[1;33mWARNING: the two methods disagree by >15 MeV -- do not quote a number yet.\033[0m\n");
   }
   printf("==========================================\n");
}
