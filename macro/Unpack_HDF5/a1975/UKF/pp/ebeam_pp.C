/// @file ebeam_pp.C
/// @brief Calibrate the a1975 16C beam energy from the 16C(p,p') spectrum itself, using TWO
///        conditions that a single parameter must satisfy at once.
///
/// pp/excorr_pp.C showed that the uncorrected spectrum puts the g.s. at +0.185 and the 2+ at
/// +1.862, i.e. offsets that SHRINK with Ex (+0.185 -> +0.096). A pure offset cannot do that; a
/// wrong beam energy can, because Ebeam sets the scale of the missing-mass reconstruction, not
/// just its zero. So scan Ebeam and ask for:
///
///     (1) the g.s. to sit at Ex = 0
///     (2) the g.s.-to-2+ spacing to equal the known 1.766 MeV
///
/// One parameter, two conditions -- the problem is over-determined, so agreement between the two
/// is a genuine test of the hypothesis rather than a fit. If they cross at the same Ebeam, the
/// scale error was the beam energy. If they do not, something else sets the scale and the residual
/// tells you how much.
///
/// Peaks are fitted uncorrected (no drift subtraction), because the drift schemes in excorr_pp.C
/// are themselves suspect until the scale is right.
///
///   root -b -q 'pp/ebeam_pp.C("/tmp/.../pp_kin_genfit.root")'
///   root -b -q 'pp/ebeam_pp.C("...pp_kin_ukf.root",180,205,2.5,"_ukf")'

static double om(double x, double y, double z) { return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z); }

static double exOf(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3;
   double s = m1*m1 + m2*m2 + 2*m2*E1, u = m2*m2 + m3*m3 - 2*m2*E3;
   double a = (cos(th)*om(s,m1*m1,m2*m2)*om(u,m2*m2,m3*m3) - (s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)
              + s + u - m2*m2;
   return a < 0 ? std::nan("") : sqrt(a) - m4;
}

static bool fitPk(TH1F *h, double c, double win, double &mu, double &er, double &sg)
{
   if (h->GetEntries() < 80) return false;
   TF1 g("g", "gaus(0)+pol1(3)", c-win, c+win);
   int bm = h->GetMaximumBin();
   g.SetParameters(h->GetBinContent(bm), c, 0.30, 0, 0);
   g.SetParLimits(1, c-win, c+win); g.SetParLimits(2, 0.03, 1.5);
   if (h->Fit(&g, "QRN") != 0) return false;
   mu = g.GetParameter(1); er = g.GetParError(1); sg = g.GetParameter(2);
   return er > 0 && er < 0.5;
}

/// Masses default to 16C(p,p')16C. For the D2-target elastic 16C(d,d)16C pass
/// mEject = mTarg = 2.0135532 and leave the rest; the 16C 2+ at 1.766 is the ruler either way,
/// since the residual is 16C in both channels.
void ebeam_pp(TString cache, double eLo = 180, double eHi = 205, double dE = 2.5, TString outTag = "",
              double chi2max = 5.0, double icMin = 950, double icMax = 1350,
              double ex2p = 1.766, double thLo = 30, double thHi = 95, double win2p = 1.2,
              double mBeamAmu = 16.0147013, double mTargAmu = 1.00782503,
              double mEjectAmu = 1.00782503, double mResidAmu = 16.0147013)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = mBeamAmu*u, m2 = mTargAmu*u, m3 = mEjectAmu*u, m4 = mResidAmu*u;
   TString dir = gSystem->DirName(__FILE__);

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree `pk`\n"); return; }
   // the (p,p') caches call the vertex column `vz` and carry an ion-chamber value; the D2 (d,d)
   // cache calls it `vertexz` and has no IC at all (the D2 reco has no _FRIB.root)
   float ke, theta, vz = 0, chi2ndf, ic = -1;
   t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta);
   t->SetBranchAddress("chi2ndf",&chi2ndf);
   if (t->GetBranch("vz")) t->SetBranchAddress("vz",&vz);
   else if (t->GetBranch("vertexz")) t->SetBranchAddress("vertexz",&vz);
   const bool hasIC = t->GetBranch("ic") != nullptr;
   if (hasIC) t->SetBranchAddress("ic",&ic);
   if (!hasIC && icMin > 0) {
      printf("note: no `ic` branch in this cache -- IC gate disabled\n");
      icMin = -1;
   }
   std::vector<float> vk, vt;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (chi2ndf > chi2max || ke <= 0) continue;
      if (icMin > 0 && (ic < icMin || ic > icMax)) continue;
      if (theta < thLo || theta >= thHi) continue;
      vk.push_back(ke); vt.push_back(theta);
   }
   f->Close();
   printf("\n=== ebeam_pp: %s ===\n%zu protons, theta[%.0f,%.0f]\n",
          gSystem->BaseName(cache), vk.size(), thLo, thHi);
   printf("%-9s %9s %9s %9s %9s %11s %11s\n",
          "Ebeam", "gs_mu", "gs_sig", "2+_mu", "2+_sig", "spacing", "spac-1.766");
   printf("%s\n", TString('-', 76).Data());

   TGraph *gGs = new TGraph(), *gSp = new TGraph();
   for (double E = eLo; E <= eHi + 1e-9; E += dE) {
      TH1F h("h","",400,-8,16); h.SetDirectory(nullptr);
      for (size_t i = 0; i < vk.size(); i++) {
         double e = exOf(m1,m2,m3,m4,E,vt[i]*TMath::DegToRad(),vk[i]);
         if (!std::isnan(e)) h.Fill(e);
      }
      double g1,g1e,g1s, c1,c1e,c1s;
      bool og = fitPk(&h, 0.0, 1.2, g1,g1e,g1s);
      // the 2+ MOVES as Ebeam is scanned, so the window must be wide enough that the fit is not
      // pinned to its edge -- an edge-pinned mu fakes a saturating spacing curve.
      bool oc = fitPk(&h, ex2p, win2p, c1,c1e,c1s);
      if (oc && fabs(c1 - ex2p) > 0.92*win2p) oc = false;   // pinned: reject
      printf("%-9.1f %9s %9s %9s %9s %11s %11s\n", E,
             og?Form("%+.3f",g1):"-", og?Form("%.3f",g1s):"-",
             oc?Form("%+.3f",c1):"-", oc?Form("%.3f",c1s):"-",
             (og&&oc)?Form("%.3f",c1-g1):"-", (og&&oc)?Form("%+.3f",c1-g1-ex2p):"-");
      if (og) gGs->SetPoint(gGs->GetN(), E, g1);
      if (og && oc) gSp->SetPoint(gSp->GetN(), E, c1-g1-ex2p);
   }

   auto zeroOf = [](TGraph *g) -> double {          // linear crossing of the last sign change
      for (int i = 0; i + 1 < g->GetN(); i++) {
         double x0,y0,x1,y1; g->GetPoint(i,x0,y0); g->GetPoint(i+1,x1,y1);
         if (y0 == 0) return x0;
         if (y0*y1 < 0) return x0 + (x1-x0)*(-y0)/(y1-y0);
      }
      return std::nan("");
   };
   double eGs = zeroOf(gGs), eSp = zeroOf(gSp);
   printf("%s\n", TString('-', 76).Data());
   printf("Ebeam putting the g.s. at zero        : %s MeV\n",
          std::isnan(eGs) ? "not bracketed in this range" : Form("%.2f", eGs));
   printf("Ebeam giving the correct 1.766 spacing: %s MeV\n",
          std::isnan(eSp) ? "not bracketed in this range" : Form("%.2f", eSp));
   if (!std::isnan(eGs) && !std::isnan(eSp))
      printf("  -> the two conditions %s (difference %.2f MeV)\n",
             fabs(eGs-eSp) < 2.0 ? "AGREE: a beam-energy error explains the scale"
                                 : "DISAGREE: Ebeam alone does NOT set this scale",
             eGs - eSp);

   TCanvas *c = new TCanvas("ceb","ebeam_pp",1200,520); c->Divide(2,1);
   c->cd(1); gGs->SetMarkerStyle(20); gGs->SetTitle("g.s. position vs E_{beam};E_{beam} [MeV];E_{x}^{g.s.} [MeV]");
   gGs->Draw("APL"); { auto*z=new TLine(eLo,0,eHi,0); z->SetLineStyle(2); z->SetLineColor(kRed); z->Draw(); }
   c->cd(2); gSp->SetMarkerStyle(20);
   gSp->SetTitle("spacing error vs E_{beam};E_{beam} [MeV];(2^{+}-g.s.) - 1.766 [MeV]");
   gSp->Draw("APL"); { auto*z=new TLine(eLo,0,eHi,0); z->SetLineStyle(2); z->SetLineColor(kRed); z->Draw(); }
   gSystem->Exec("mkdir -p " + dir + "/plots");
   c->SaveAs(dir + "/plots/ebeam_pp" + outTag + ".png");
}
