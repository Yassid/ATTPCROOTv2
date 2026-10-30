/// @file pid_dd_referee.C
/// @brief Referee the PID with KINEMATICS: tracks on the 16C(d,d) elastic locus under the deuteron
///        hypothesis with BOTH estimators' geometry (Spyral and 3D: theta_lab, Brho) are deuterons,
///        whatever their dE/dx says. The spread of each estimator's ln m on that sample -- tails
///        included -- is its real deuteron PID resolution. Same idea as the H2 elastic-proton referee.
///   root -l -b -q 'pid_dd_referee.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016.root")'
namespace {
double exDD(double Eb, double th, double brho)
{
   const double u = 931.49410242, me = 0.51099895;
   const double M = 16.014701 * u - 6 * me, md = 2.01410178 * u - me; // 16C atomic-6e, d nuclear
   double p = 299.792458 * brho, Ed = std::sqrt(p * p + md * md);
   double Et = Eb + M, pb = std::sqrt(Et * Et - M * M);
   double E4 = Et + md - Ed, px = -p * std::sin(th), pz = pb - p * std::cos(th);
   double m2 = E4 * E4 - px * px - pz * pz;
   return m2 > 0 ? std::sqrt(m2) - M : -999;
}
} // namespace
void pid_dd_referee(TString in, double Eb = 184.25, double exWin = 1.5)
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   int va, v3, nh;
   double br, bw, b3, sqa, sqw, sq3a, sq3w, sq3c, sq3cw, sq3f, sq3fw, hybB, hybBw, sq3t, sq3tw, sq3h, sq3hw, hybC, hybCw, pol, th3, sp3;
   int vw;
   w->SetBranchAddress("va", &va); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("nh", &nh);
   w->SetBranchAddress("br", &br); w->SetBranchAddress("b3", &b3); w->SetBranchAddress("sqa", &sqa);
   w->SetBranchAddress("sq3a", &sq3a); w->SetBranchAddress("sq3c", &sq3c); w->SetBranchAddress("pol", &pol);
   w->SetBranchAddress("th3", &th3); w->SetBranchAddress("sp3", &sp3);
   w->SetBranchAddress("bw", &bw); w->SetBranchAddress("sqw", &sqw); w->SetBranchAddress("sq3w", &sq3w);
   w->SetBranchAddress("sq3cw", &sq3cw); w->SetBranchAddress("sq3f", &sq3f); w->SetBranchAddress("sq3fw", &sq3fw);
   w->SetBranchAddress("hybB", &hybB); w->SetBranchAddress("hybBw", &hybBw); w->SetBranchAddress("sq3t", &sq3t);
   w->SetBranchAddress("sq3tw", &sq3tw); w->SetBranchAddress("sq3h", &sq3h); w->SetBranchAddress("sq3hw", &sq3hw);
   w->SetBranchAddress("hybC", &hybC); w->SetBranchAddress("hybCw", &hybCw); w->SetBranchAddress("vw", &vw);
   // medians for the same normalisation as the other macros
   auto med = [&](const char *v) {
      w->Draw(v, "va&&vw&&v3&&sp3<6&&nh>=30&&sq3c>0&&sq3cw>0&&sq3f>0&&sq3fw>0&&hybB>0&&hybBw>0&&sq3t>0&&sq3tw>0&&sq3h>0&&sq3hw>0&&hybC>0&&hybCw>0", "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      return x[x.size() / 2];
   };
   double mS = med("sqa"), m3 = med("sq3a"), mC = med("sq3c"), mSw = med("sqw"), m3w = med("sq3w"), mCw = med("sq3cw"), mF = med("sq3f"), mFw = med("sq3fw");
   double mB = med("hybB"), mBw = med("hybBw"), mT = med("sq3t"), mTw = med("sq3tw");
   double mH = med("sq3h"), mHw = med("sq3hw"), mC2 = med("hybC"), mC2w = med("hybCw");
   std::vector<double> LH, LHw, LC2, LC2w;
   std::vector<double> LS, L3, LC, LSw, L3w, LCw, LF, LFw, LA, LAw, LB, LBw, LT, LTw, TH;
   long nOnS = 0, nOn3 = 0;
   for (Long64_t i = 0; i < w->GetEntries(); ++i) {
      w->GetEntry(i);
      if (!va || !vw || !v3 || sp3 >= 6 || nh < 30 || sq3c <= 0 || sq3cw <= 0 || sq3f <= 0 || sq3fw <= 0 || hybB <= 0 || hybBw <= 0 ||
          sq3t <= 0 || sq3tw <= 0 || sq3h <= 0 || sq3hw <= 0 || hybC <= 0 || hybCw <= 0)
         continue;
      double thS = (180 - pol) * TMath::DegToRad();
      bool onS = std::abs(exDD(Eb, thS, br)) < exWin, on3 = std::abs(exDD(Eb, th3 * TMath::DegToRad(), b3)) < exWin;
      nOnS += onS; nOn3 += on3;
      if (!(onS && on3)) continue;
      LS.push_back(std::log(sqa / mS * br)); L3.push_back(std::log(sq3a / m3 * b3)); LC.push_back(std::log(sq3c / mC * b3));
      LSw.push_back(std::log(sqw / mSw * bw)); L3w.push_back(std::log(sq3w / m3w * b3)); LCw.push_back(std::log(sq3cw / mCw * b3));
      LF.push_back(std::log(sq3f / mF * b3)); LFw.push_back(std::log(sq3fw / mFw * b3));
      LA.push_back(std::log(sqa / mS * b3)); LAw.push_back(std::log(sqw / mSw * b3));
      LB.push_back(std::log(hybB / mB * b3)); LBw.push_back(std::log(hybBw / mBw * b3));
      LT.push_back(std::log(sq3t / mT * b3)); LTw.push_back(std::log(sq3tw / mTw * b3));
      LH.push_back(std::log(sq3h / mH * b3)); LHw.push_back(std::log(sq3hw / mHw * b3));
      LC2.push_back(std::log(hybC / mC2 * b3)); LC2w.push_back(std::log(hybCw / mC2w * b3));
      TH.push_back(th3);
   }
   printf("tracks on the (d,d) locus: Spyral geometry %ld, 3D geometry %ld, BOTH %zu (= the deuteron referee sample)\n",
          nOnS, nOn3, LS.size());
   auto rep = [&](const char *nm, std::vector<double> x) {
      int n = x.size();
      std::sort(x.begin(), x.end());
      double m = x[n / 2], hw = 0.5 * (x[n * 84 / 100] - x[n * 16 / 100]);
      int far = 0, prot = 0;
      for (double v : x) { if (std::abs(v - m) > 0.3) ++far; if (v < m - 0.45) ++prot; }
      printf("  %-26s median %+.3f  half-width %.3f  |ln m - med| > 0.3: %5.1f%%   would look like a PROTON (< med-0.45): %4.1f%%\n",
             nm, m, hw, 100.0 * far / n, 100.0 * prot / n);
   };
   rep("Spyral amp", LS); rep("Spyral width", LSw); rep("3D amp", L3); rep("3D width", L3w);
   rep("3D centreline", LC); rep("3D cl width", LCw); rep("3D first arc", LF); rep("3D 1st arc width", LFw);
   rep("HYB A: Spy dEdx, 3D Brho", LA); rep("HYB A width", LAw); rep("HYB B: Spy Q / 3D path", LB);
   rep("HYB B width", LBw); rep("3D trunc mean 70%", LT); rep("3D trunc mean width", LTw);
   rep("3D Q / HELIX path", LH); rep("3D Q / helix path width", LHw); rep("HYB C: Spy Q / helix path", LC2);
   rep("HYB C width", LC2w);
   printf("theta_lab of the referee deuterons: ");
   std::sort(TH.begin(), TH.end());
   int n = TH.size();
   if (n > 10) printf("p10 %.0f med %.0f p90 %.0f deg\n", TH[n / 10], TH[n / 2], TH[9 * n / 10]);
}
