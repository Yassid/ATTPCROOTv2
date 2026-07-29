/// @file ket_C15.C
/// @brief KE vs theta_lab matrix for a2091 15C(p,p'), one fitter per panel, from the
///        proton_kin caches written by ex_C15.C. Reads only the small ntuples, so it is
///        instant and can be re-binned freely.
///
///   root -b -q 'pp/ket_C15.C()'                                  // ukf vs genfit_nomat, gated
///   root -b -q 'pp/ket_C15.C("plots/proton_kin_g_ukf.root","plots/proton_kin_g_genfit_nomat.root")'
///
/// The red curve is the ELASTIC RECOIL LOCUS at Ebeam, i.e. what a 15C(p,p) elastic proton must
/// follow:  T_p = 2 m_p p1^2 cos^2(th) / [ (E1+m_p)^2 - p1^2 cos^2(th) ].  Any KE-vs-theta ridge
/// that departs from it is either a wrong Ebeam or a reconstruction bias, so overlaying it is the
/// quickest sanity check on both fitters at once.
///
/// Marked on the plot, because both bit this analysis on the ungated data:
///   - the KE ~ 0.79 MeV reconstruction threshold (dashed): above theta ~ 80 deg the measured
///     ridge flattens onto it while the elastic locus keeps falling.
///   - theta 88-94 deg is a Spyral `polar` blind spot; the fitted theta here is a different
///     quantity, but the region deserves suspicion.

static double TrecKet(double T1, double m1, double m2, double th)
{
   double E1 = T1 + m1, p = T1*T1 + 2*T1*m1, c = std::cos(th)*std::cos(th);
   return 2*m2*p*c / ((E1+m2)*(E1+m2) - p*c);
}

static TH2F *fill(TString path, TString hname, TString title, double chi2Cut, double thMax, double keMax,
                  long &nUsed)
{
   nUsed = 0;
   if (path.IsNull() || path.EndsWith("/")) return nullptr; // one-panel mode: second cache omitted
   TFile *f = TFile::Open(path);
   if (!f || f->IsZombie()) { printf("  MISSING %s\n", path.Data()); return nullptr; }
   TNtuple *n = (TNtuple *)f->Get("pk");
   if (!n) { printf("  no `pk` ntuple in %s\n", path.Data()); f->Close(); return nullptr; }
   auto *h = new TH2F(hname, title + ";#theta_{lab} [deg];KE_{p} [MeV]", 200, 0, thMax, 200, 0, keMax);
   h->SetDirectory(nullptr);
   float ke, theta, vz, thcm, ex, c2;
   n->SetBranchAddress("ke",&ke); n->SetBranchAddress("theta",&theta); n->SetBranchAddress("vertexz",&vz);
   n->SetBranchAddress("thcm",&thcm); n->SetBranchAddress("ex",&ex); n->SetBranchAddress("chi2ndf",&c2);
   nUsed = 0;
   for (Long64_t i = 0; i < n->GetEntries(); i++) {
      n->GetEntry(i);
      if (c2 > chi2Cut || ke <= 0) continue;
      h->Fill(theta, ke); nUsed++;
   }
   printf("  %-34s %9ld tracks (of %lld)\n", path.Data(), nUsed, n->GetEntries());
   f->Close();
   return h;
}

void ket_C15(TString cacheA = "plots/proton_kin_g_ukf.root", TString cacheB = "plots/proton_kin_g_genfit_nomat.root",
             TString labelA = "UKF (gated)", TString labelB = "GENFIT no-matFX (gated)", double Ebeam = 157.0,
             double chi2Cut = 5.0, double thMax = 120, double keMax = 25, TString outTag = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double u = 931.49401, m1 = 15.0105993*u, m2 = 1.007825*u;

   TString dir = gSystem->DirName(__FILE__);
   auto full = [&](TString p) { return p.BeginsWith("/") ? p : dir + "/" + p; };

   printf("\n=== KE vs theta, Ebeam = %.0f MeV ===\n", Ebeam);
   long nA = 0, nB = 0;
   TH2F *hA = fill(full(cacheA), "hA", labelA, chi2Cut, thMax, keMax, nA);
   TH2F *hB = fill(full(cacheB), "hB", labelB, chi2Cut, thMax, keMax, nB);
   if (!hA && !hB) { printf("nothing to plot\n"); return; }

   // elastic locus + the threshold line
   auto *loc = new TF1("loc", [&](double *x, double *) {
      return TrecKet(Ebeam, m1, m2, x[0]*TMath::DegToRad()); }, 1, 89.5, 0);
   loc->SetNpx(500); loc->SetLineColor(kRed); loc->SetLineWidth(2);
   auto thrLine = [&]() {
      auto *l = new TLine(0, 0.79, thMax, 0.79);
      l->SetLineColor(kGray+2); l->SetLineStyle(2); l->Draw(); return l; };

   int npad = (hA && hB) ? 2 : 1;
   TCanvas *c = new TCanvas("cket", "KE vs theta", npad == 2 ? 1600 : 900, 700);
   if (npad == 2) c->Divide(2, 1);
   int ip = 1;
   for (TH2F *h : {hA, hB}) {
      if (!h) continue;
      c->cd(ip++); gPad->SetLogz(); gPad->SetRightMargin(0.13);
      h->Draw("colz"); loc->Draw("same"); thrLine();
      auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.030); tx->SetTextColor(kRed);
      tx->DrawLatex(0.16, 0.86, Form("red: elastic locus, E_{beam} = %.0f MeV", Ebeam));
      tx->SetTextColor(kGray+3);
      tx->DrawLatex(0.16, 0.82, "dashed: KE = 0.79 MeV reconstruction threshold");
   }
   TString png = dir + "/plots/ket_C15" + outTag + ".png";
   c->SaveAs(png);
   printf("wrote %s\n", png.Data());

   // the ridge as numbers, so the two fitters can be compared without eyeballing
   if (hA && hB) {
      printf("\n-- most-probable KE per theta slice (peak bin), both fitters --\n");
      printf("%-8s %12s %12s %10s   %s\n", "theta", labelA.Data(), labelB.Data(), "elastic", "A-B");
      for (double t = 60; t <= 92; t += 4) {
         int b1 = hA->GetXaxis()->FindBin(t), b2 = hA->GetXaxis()->FindBin(t + 4) - 1;
         auto *pA = hA->ProjectionY("pA", b1, b2);
         auto *pB = hB->ProjectionY("pB", b1, b2);
         double kA = pA->GetEntries() > 30 ? pA->GetBinCenter(pA->GetMaximumBin()) : -1;
         double kB = pB->GetEntries() > 30 ? pB->GetBinCenter(pB->GetMaximumBin()) : -1;
         double kel = TrecKet(Ebeam, m1, m2, (t + 2) * TMath::DegToRad());
         printf("%-8s %12s %12s %10.3f   %s\n", Form("%.0f-%.0f", t, t + 4),
                kA > 0 ? Form("%.3f", kA) : "-", kB > 0 ? Form("%.3f", kB) : "-", kel,
                (kA > 0 && kB > 0) ? Form("%+.3f", kA - kB) : "-");
         delete pA; delete pB;
      }
   }
}
