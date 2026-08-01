/// @file excorr_pp.C
/// @brief Flatten the Ex-vs-angle drift of the a1975 16C(p,p') elastic peak -- and check that the
///        flattening is a real resolution gain rather than a self-referential collapse.
///
/// This merges the two halves of the lesson already learned elsewhere in this repo:
///
///   * from pp/excorr_pd.C : subtract the MEASURED drift interpolated through the per-angle points.
///     The drift here is an inverted-U with a back-angle drop (pp/thetadiag_pp.C: -0.7 MeV at 30-40
///     deg, +0.14 at 65-75, falling after; a straight line gives chi2/ndf ~ 17), so a polynomial is
///     the wrong model and the measured curve is the right one.
///
///   * from a2091's pp/tiltcorr_C15.C : WHICH angle you correct in decides whether the result means
///     anything. theta_lab is measured independently of Ex, so removing a trend in it removes real
///     instrumental smearing. theta_cm is computed FROM Ex, so subtracting a function of it
///     subtracts a function of each event's own Ex -- the peak narrows while the LEVEL SPACINGS
///     shrink. That is compression, not resolution.
///
/// 16C carries its own ruler: the 2+ at 1.766 MeV. So both schemes are judged on the g.s.-to-2+
/// SPACING, which a genuine correction leaves at 1.766 while a compressing one pulls in. Both are
/// computed and printed side by side; theta_lab is the one to trust when they disagree.
///
/// The elastic peak is isolated ON ITS KINEMATIC LOCUS (|KE - KE_el(theta)| < keWin) rather than
/// with a flat Ex window, so inelastic strength cannot leak into the reference and bend the drift.
///
///   root -b -q 'pp/excorr_pp.C("/tmp/.../pp_kin_genfit.root")'
///   root -b -q 'pp/excorr_pp.C("...pp_kin_ukf.root",192,5,950,1350,3.0,"_ukf")'

static double om(double x, double y, double z) { return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z); }

struct Kin { double ex, thcm; };
static Kin kine(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3, E4 = E1 + m2 - E3;
   double s = m1*m1 + m2*m2 + 2*m2*E1, u = m2*m2 + m3*m3 - 2*m2*E3;
   double a = (cos(th)*om(s,m1*m1,m2*m2)*om(u,m2*m2,m3*m3) - (s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)
              + s + u - m2*m2;
   if (a < 0) return {std::nan(""), std::nan("")};
   double m4x = sqrt(a), ex = m4x - m4, t = m2*m2 + m4x*m4x - 2*m2*E4;
   double g = (s*s + s*(2*t - m1*m1 - m2*m2 - m3*m3 - m4x*m4x) + (m1*m1 - m2*m2)*(m3*m3 - m4x*m4x))
              / (om(s,m1*m1,m2*m2)*om(s,m3*m3,m4x*m4x));
   if (g < -1) g = -1; if (g > 1) g = 1;
   return {ex, (TMath::Pi() - acos(g))*TMath::RadToDeg()};
}

/// KE at which Ex(KE) = exTarget for this theta_lab (Ex falls monotonically with KE -> bisect)
static double keAtEx(double m1, double m2, double m3, double m4, double E, double th, double exTarget)
{
   double lo = 0.2, hi = 120.0;
   double flo = kine(m1,m2,m3,m4,E,th,lo).ex, fhi = kine(m1,m2,m3,m4,E,th,hi).ex;
   if (std::isnan(flo) || std::isnan(fhi) || (flo-exTarget)*(fhi-exTarget) > 0) return std::nan("");
   for (int i = 0; i < 80; i++) {
      double mid = 0.5*(lo+hi), fm = kine(m1,m2,m3,m4,E,th,mid).ex;
      if (std::isnan(fm)) { hi = mid; continue; }
      if ((flo-exTarget)*(fm-exTarget) <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
   }
   return 0.5*(lo+hi);
}

static bool fitPk(TH1F *h, double c, double win, double &mu, double &er, double &sg)
{
   if (h->GetEntries() < 80) return false;
   TF1 g("g", "gaus(0)+pol1(3)", c-win, c+win);
   int bm = h->GetMaximumBin();
   g.SetParameters(h->GetBinContent(bm), c, 0.35, 0, 0);
   g.SetParLimits(1, c-win, c+win); g.SetParLimits(2, 0.03, 1.5);
   if (h->Fit(&g, "QRN") != 0) return false;
   mu = g.GetParameter(1); er = g.GetParError(1); sg = g.GetParameter(2);
   return er > 0 && er < 0.5;
}

void excorr_pp(TString cache, double Ebeam = 192.0, double chi2max = 5.0, double icMin = 950,
               double icMax = 1350, double keWin = 3.0, TString outTag = "",
               double ex2p = 1.766, double thLo = 30, double thHi = 95, int nProf = 13)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = 16.0147013*u, m2 = 1.00782503*u, m3 = m2, m4 = m1;   // 16C(p,p)16C
   TString dir = gSystem->DirName(__FILE__);

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree `pk` in %s\n", cache.Data()); return; }
   float ke, theta, vz, chi2ndf, ic;
   t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta); t->SetBranchAddress("vz",&vz);
   t->SetBranchAddress("chi2ndf",&chi2ndf); t->SetBranchAddress("ic",&ic);

   struct Ev { double ke, th, vz, ex, thcm; };
   std::vector<Ev> ev;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (chi2ndf > chi2max || ke <= 0) continue;
      if (icMin > 0 && (ic < icMin || ic > icMax)) continue;
      Kin k = kine(m1,m2,m3,m4,Ebeam,theta*TMath::DegToRad(),ke);
      if (std::isnan(k.ex)) continue;
      ev.push_back({ke, theta, vz, k.ex, k.thcm});
   }
   f->Close();
   printf("\n=== excorr_pp: %s ===\nselected %zu protons (chi2/ndf<%.1f, IC[%.0f,%.0f]), Ebeam %.1f\n",
          gSystem->BaseName(cache), ev.size(), chi2max, icMin, icMax, Ebeam);

   // ---- 1. measured elastic drift vs theta_lab, isolated ON THE LOCUS ----
   printf("\n-- elastic peak on its kinematic locus (|KE - KE_el(theta)| < %.1f MeV) --\n", keWin);
   printf("%-12s %8s %9s %9s %8s %8s\n", "theta_lab", "N", "KE_el", "mu_Ex", "err", "sigma");
   TGraph *gLab = new TGraph();
   TGraph *gKE  = new TGraph();   // implied KE bias vs KE, for the recalibration variant
   double bw = (thHi - thLo)/nProf;
   for (int b = 0; b < nProf; ++b) {
      double t0 = thLo + b*bw, t1 = t0 + bw, tmid = 0.5*(t0+t1);
      double keel = keAtEx(m1,m2,m3,m4,Ebeam,tmid*TMath::DegToRad(),0.0);
      if (std::isnan(keel)) continue;
      TH1F h("h","",300,-6,10); h.SetDirectory(nullptr);
      for (auto &e : ev) if (e.th >= t0 && e.th < t1 && fabs(e.ke - keel) < keWin) h.Fill(e.ex);
      double mu, er, sg;
      if (fitPk(&h, 0.0, 1.4, mu, er, sg)) {
         printf("%-12s %8.0f %9.2f %9.3f %8.3f %8.3f\n", Form("%.0f-%.0f",t0,t1), h.GetEntries(), keel, mu, er, sg);
         gLab->SetPoint(gLab->GetN(), tmid, mu);
         // the KE bias that would remove this slice's offset: dKE = -mu / (dEx/dKE)
         double tr = tmid*TMath::DegToRad();
         double d = (kine(m1,m2,m3,m4,Ebeam,tr,keel+0.25).ex - kine(m1,m2,m3,m4,Ebeam,tr,keel-0.25).ex)/0.5;
         if (fabs(d) > 1e-6) gKE->SetPoint(gKE->GetN(), keel, -mu/d);
      } else
         printf("%-12s %8.0f %9.2f %9s\n", Form("%.0f-%.0f",t0,t1), h.GetEntries(), keel, "(no fit)");
   }
   if (gLab->GetN() < 3) { printf("too few theta_lab slices for a drift curve\n"); return; }

   // ---- 2. same, vs theta_cm (the unsafe variable -- computed for comparison only) ----
   TGraph *gCm = new TGraph();
   for (double c = 20; c < 170; c += 12) {
      TH1F h("h","",300,-6,10); h.SetDirectory(nullptr);
      for (auto &e : ev) {
         double keel = keAtEx(m1,m2,m3,m4,Ebeam,e.th*TMath::DegToRad(),0.0);
         if (std::isnan(keel) || fabs(e.ke - keel) > keWin) continue;
         if (e.thcm >= c && e.thcm < c+12) h.Fill(e.ex);
      }
      double mu, er, sg;
      if (fitPk(&h, 0.0, 1.4, mu, er, sg)) gCm->SetPoint(gCm->GetN(), c+6, mu);
   }

   // measured drift, interpolated through the points and clamped outside their range
   auto mk = [](TGraph *g) -> std::function<double(double)> {
      double x0, x1, y;
      g->GetPoint(0, x0, y); g->GetPoint(g->GetN()-1, x1, y);
      return [g, x0, x1](double x) { return g->Eval(std::min(x1, std::max(x0, x))); };
   };
   std::function<double(double)> driftLab = mk(gLab);
   std::function<double(double)> driftCm;
   if (gCm->GetN() >= 3) driftCm = mk(gCm);
   // gKE was filled in theta order, i.e. DECREASING KE; TGraph::Eval interpolates assuming
   // ascending x, so it must be sorted before it is used as a lookup.
   gKE->Sort();
   std::function<double(double)> dkeOf;
   if (gKE->GetN() >= 3) dkeOf = mk(gKE);
   printf("\nmeasured drift: %d theta_lab points, %d theta_cm points\n", gLab->GetN(), gCm->GetN());

   // ---- 3. apply each and judge on the g.s.-to-2+ SPACING, not the width ----
   printf("\n-- effect of each correction (all events, no locus cut) --\n");
   printf("%-24s %9s %9s %9s %9s %11s\n", "scheme", "gs_mu", "gs_sigma", "2+_mu", "2+_sigma", "spacing");
   auto report = [&](const char *name, int mode) {
      TH1F h("h","",400,-8,16); h.SetDirectory(nullptr);
      for (auto &e : ev) {
         if (mode == 3) {   // recalibrate KE itself and recompute Ex -- the physical model
            Kin k = kine(m1,m2,m3,m4,Ebeam,e.th*TMath::DegToRad(), e.ke + dkeOf(e.ke));
            if (!std::isnan(k.ex)) h.Fill(k.ex);
            continue;
         }
         double c = 0;
         if (mode == 1) c = driftLab(e.th);
         if (mode == 2 && gCm->GetN() >= 3) c = driftCm(e.thcm);
         h.Fill(e.ex - c);
      }
      double g1,g1e,g1s, c1,c1e,c1s;
      bool og = fitPk(&h, 0.0, 1.0, g1,g1e,g1s);
      bool oc = fitPk(&h, ex2p, 0.8, c1,c1e,c1s);
      printf("%-24s %9s %9s %9s %9s %11s\n", name,
             og?Form("%+.3f",g1):"-", og?Form("%.3f",g1s):"-",
             oc?Form("%+.3f",c1):"-", oc?Form("%.3f",c1s):"-",
             (og&&oc)?Form("%.3f",c1-g1):"-");
   };
   report("none", 0);
   report("minus a(theta_lab)", 1);
   if (gCm->GetN() >= 3) report("minus g(theta_cm)", 2);
   if (gKE->GetN() >= 3) report("KE recalib dKE(KE)", 3);
   printf("\nJudge on SPACING: the 16C 2+ sits %.3f MeV above the g.s. A scheme that narrows the\n", ex2p);
   printf("peaks but shrinks the spacing is compressing the spectrum, not resolving it.\n");
   printf("NOTE a(theta_lab) and dKE(KE) are near-degenerate ON the elastic locus (KE_el is a\n");
   printf("monotonic function of theta there), so ONLY the inelastic 2+ can tell them apart.\n");

   // ---- 4. plots ----
   TCanvas *c = new TCanvas("cex","excorr_pp",1500,520); c->Divide(3,1);
   c->cd(1);
   gLab->SetTitle("elastic peak vs #theta_{lab} (subtracted);#theta_{lab} [deg];E_{x} [MeV]");
   gLab->SetMarkerStyle(20); gLab->SetMarkerColor(kBlue+1); gLab->Draw("APL");
   auto *z = new TLine(gLab->GetX()[0], 0, gLab->GetX()[gLab->GetN()-1], 0);
   z->SetLineStyle(2); z->Draw();
   const char *nm[2] = {"none", "minus a(#theta_{lab})"};
   for (int mode = 0; mode < 2; mode++) {
      TH1F *h = new TH1F(Form("hh%d",mode), Form("%s;E_{x}(^{16}C) [MeV];counts",nm[mode]), 300, -6, 12);
      for (auto &e : ev) h->Fill(e.ex - (mode == 1 ? driftLab(e.th) : 0.0));
      c->cd(mode+2); h->SetLineColor(mode ? kBlue+1 : kGray+2); h->SetLineWidth(2); h->Draw("hist");
      auto *l = new TLine(ex2p, 0, ex2p, h->GetMaximum()); l->SetLineColor(kRed); l->SetLineStyle(2); l->Draw();
   }
   gSystem->Exec("mkdir -p " + dir + "/plots");
   c->SaveAs(dir + "/plots/excorr_pp" + outTag + ".png");
   printf("wrote plots/excorr_pp%s.png\n", outTag.Data());
}
