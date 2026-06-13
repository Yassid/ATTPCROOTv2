/// @file overlay_kin_pd.C
/// @brief 16C(p,d)15C kinematics with theoretical KE(theta_lab) curves overlaid for
///        the g.s. and excited states — clean open lines (the full kinematic loop).
///
/// Reads cached (ke, theta) from deuteron_kin.root; draws the data KE-vs-theta_lab
/// plus relativistic two-body curves for Ex = 0, 0.74, 3.10, 4.66 MeV.
///
///   root -b -q 'pd/overlay_kin_pd.C'

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

void overlay_kin_pd(double Ebeam = 192.0, TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double u = 931.49401, m1 = 16.0147 * u, m2 = 1.007825 * u, m3 = 2.01410178 * u, m4_0 = 15.0105993 * u;

   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   TH2F *h = new TH2F("h", "^{16}C(p,d)^{15}C kinematics + theory;#theta_{lab} [deg];KE_{d} [MeV]", 250, 8, 45, 250, 0,
                      50);
   t->Project("h", "ke:theta");

   // theoretical KE_d(theta_lab) for residual mass m4 = m4_0 + Ex: closed loop
   auto curve = [&](double Ex, int color, int style) -> TGraph * {
      double E1 = Ebeam + m1, p1 = std::sqrt(E1 * E1 - m1 * m1), Etot = E1 + m2, m4 = m4_0 + Ex;
      std::vector<double> up_t, up_k, lo_t, lo_k;
      for (double thd = 0.0; thd <= 50.0; thd += 0.1) {
         double c = std::cos(thd * TMath::DegToRad());
         double A = Etot * Etot - m4 * m4 - p1 * p1 + m3 * m3;
         double a = 4 * Etot * Etot - 4 * p1 * p1 * c * c, b = -4 * Etot * A,
                cc = A * A + 4 * p1 * p1 * c * c * m3 * m3;
         double disc = b * b - 4 * a * cc;
         if (disc < 0 || a == 0)
            continue;
         for (int s = -1; s <= 1; s += 2) {
            double E3 = (-b + s * std::sqrt(disc)) / (2 * a);
            if (E3 <= m3 || (2 * Etot * E3 - A) * c <= 0)
               continue;
            if (s > 0) {
               up_t.push_back(thd);
               up_k.push_back(E3 - m3);
            } else {
               lo_t.push_back(thd);
               lo_k.push_back(E3 - m3);
            }
         }
      }
      TGraph *g = new TGraph();
      for (size_t i = 0; i < up_t.size(); ++i)
         g->SetPoint(g->GetN(), up_t[i], up_k[i]); // upper branch forward
      for (size_t i = lo_t.size(); i-- > 0;)
         g->SetPoint(g->GetN(), lo_t[i], lo_k[i]); // lower branch back
      g->SetLineColor(color);
      g->SetLineWidth(2);
      g->SetLineStyle(style);
      g->SetFillStyle(0);
      return g;
   };

   TCanvas *c = new TCanvas("c", "overlay", 950, 720);
   gPad->SetLogz();
   h->Draw("colz");
   struct {
      double ex;
      int col;
      const char *lab;
   } st[] = {{0.0, kRed, "g.s."}, {0.74, kMagenta + 2, "0.74"}, {3.10, kGreen + 2, "3.10"}, {4.66, kOrange + 7, "4.66"}};
   TLegend *leg = new TLegend(0.6, 0.66, 0.88, 0.88);
   leg->SetHeader("^{15}C state (E_{x})");
   for (auto &s : st) {
      TGraph *g = curve(s.ex, s.col, 1);
      g->Draw("L same");
      leg->AddEntry(g, s.lab, "l");
   }
   leg->Draw();
   c->SaveAs("pd/plots/overlay_kin.png");
   printf("saved pd/plots/overlay_kin.png\n");
}
