/// @file overlay_pp_elastic.C
/// @brief Overlay the genfit-fitted proton KE-vs-theta_lab on the analytic 16C(p,p)
/// ELASTIC recoil-proton kinematic line. Validates the fitter's absolute calibration
/// (B, angle, energy): elastic protons must lie on the curve (excited states / inelastic
/// fall below it; vertex energy-loss shifts the whole band slightly down from nominal).
///
///   root -b -q 'pp/overlay_pp_elastic.C("/tmp/run_0106_genfitter_pp.root")'
///
/// Recoil of a target proton (m_p, at rest) struck elastically by the 16C beam
/// (m_C, lab KE Ebeam):  T_p(psi) = 2 m_p (p1c)^2 cos^2 psi / [ (E1+m_p)^2 - (p1c)^2 cos^2 psi ]
/// with E1 = Ebeam + m_C, p1c = sqrt(E1^2 - m_C^2).  (psi = proton lab angle.)

double Tproton(double psiDeg, double Ebeam, double mC, double mp)
{
   double psi = psiDeg * TMath::DegToRad();
   double E1 = Ebeam + mC;
   double p1c2 = E1 * E1 - mC * mC;
   double c2 = std::cos(psi) * std::cos(psi);
   double num = 2.0 * mp * p1c2 * c2;
   double den = (E1 + mp) * (E1 + mp) - p1c2 * c2;
   return num / den; // MeV
}

void overlay_pp_elastic(TString file = "/tmp/run_0106_genfitter_pp.root", double Ebeam = 192.0,
                        TString outPng = "/tmp/pp_elastic_overlay.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401, mC = 16.0147 * u, mp = 1.007825 * u;

   TFile f(file);
   TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);
   TH2F *h = new TH2F("h", "16C(p,p): fitted protons vs elastic line;#theta_{lab} [deg];KE [MeV]", 90, 0, 95, 110, 0, 45);
   long n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) continue;
      for (auto &fp : ev->GetFittedTracks()) {
         if (!fp) continue;
         auto &m = fp->GetTrackMetadata();
         if (!(m && m->GetGoodFit())) continue;
         const auto &k = fp->GetKinematics(0);
         h->Fill(k.theta * TMath::RadToDeg(), k.kineticEnergy);
         ++n;
      }
   }

   // analytic elastic recoil-proton curves: nominal Ebeam + a lower one to bracket
   // the vertex energy-loss spread (beam loses energy in the gas before reacting).
   auto makeCurve = [&](double eb, Color_t col, int style) {
      auto *g = new TGraph();
      for (double th = 1; th <= 90; th += 0.5) g->SetPoint(g->GetN(), th, Tproton(th, eb, mC, mp));
      g->SetLineColor(col); g->SetLineWidth(3); g->SetLineStyle(style);
      return g;
   };
   TGraph *gNom = makeCurve(Ebeam, kRed + 1, 1);
   TGraph *gLo = makeCurve(Ebeam * 0.75, kRed + 1, 2); // ~25% beam energy loss bracket

   TCanvas *c = new TCanvas("c", "pp elastic", 950, 680);
   h->Draw("colz");
   gNom->Draw("L same");
   gLo->Draw("L same");
   TLegend *lg = new TLegend(0.45, 0.74, 0.88, 0.88);
   lg->AddEntry(h, Form("genfit protons (good fit, n=%ld)", n), "p");
   lg->AddEntry(gNom, Form("elastic, E_{beam}=%.0f MeV", Ebeam), "l");
   lg->AddEntry(gLo, Form("elastic, E_{beam}=%.0f MeV (E-loss)", Ebeam * 0.75), "l");
   lg->Draw();
   c->SaveAs(outPng);
   printf("plotted %ld good protons; elastic T_p(0)=%.1f MeV, T_p(90)=%.1f MeV; wrote %s\n", n,
          Tproton(0, Ebeam, mC, mp), Tproton(90, Ebeam, mC, mp), outPng.Data());
}
