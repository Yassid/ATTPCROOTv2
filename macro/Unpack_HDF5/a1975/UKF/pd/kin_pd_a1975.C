/// @file kin_pd_a1975.C
/// @brief Deuteron KE vs lab angle for 16C(p,d)15C — data only, fully parameterized
///        so you can zoom / rebin freely.
///
/// Reads the cached ntuple deuteron_kin.root (written by ex_pd_a1975.C; columns
/// ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx).
///
/// Batch (saves a PNG):
///   root -b -q 'pd/kin_pd_a1975.C(10,40, 0,50, 450,500)'
///
/// Interactive (keeps the canvas so you can zoom/box with the mouse):
///   root -l 'pd/kin_pd_a1975.C(10,40, 0,50, 450,500)'
///
/// Args: thLo,thHi  = theta_lab window [deg]
///       keLo,keHi  = KE window [MeV]
///       nth,nke    = bin counts in theta, KE
///       logz       = log color scale (default true)
///       what       = the y:x expression to plot (default "ke:theta"; try
///                    "ex:thcm", "brho:theta", "ke:thcm", ...)
///       cut        = optional TTree selection, e.g. "chi2ndf<2" or "ex>-2&&ex<2"
///       cacheFile  = ntuple file
///       outPng     = output image ("" => auto name)

void kin_pd_a1975(double thLo = 10, double thHi = 40, double keLo = 0, double keHi = 50, int nth = 450, int nke = 500,
                  bool logz = true, TString what = "ke:theta", TString cut = "", TString cacheFile = "deuteron_kin.root",
                  TString outPng = "")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   if (gSystem->AccessPathName(cacheFile)) {
      printf("ERROR: %s not found (run pd/ex_pd_a1975.C first to build it).\n", cacheFile.Data());
      return;
   }
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   if (!t) {
      printf("ERROR: ntuple 'dk' not in %s\n", cacheFile.Data());
      return;
   }

   TH2F *h = new TH2F("h", "^{16}C(p,d)^{15}C kinematics", nth, thLo, thHi, nke, keLo, keHi);
   // axis titles follow the y:x expression
   TObjArray *yx = what.Tokenize(":");
   TString yexpr = ((TObjString *)yx->At(0))->GetString();
   TString xexpr = ((TObjString *)yx->At(1))->GetString();
   h->SetXTitle(xexpr);
   h->SetYTitle(yexpr);

   t->Project("h", what, cut);
   printf("plotted %.0f entries  (%s%s)\n", h->GetEntries(), what.Data(),
          cut.Length() ? Form("  cut: %s", cut.Data()) : "");

   TCanvas *c = new TCanvas("ckin", "kin_pd", 1100, 850);
   if (logz)
      gPad->SetLogz();
   h->Draw("colz");

   if (outPng.Length() == 0)
      outPng = Form("pd/plots/kin_pd_%s_%gto%g_%gto%g.png", yexpr.Data(), thLo, thHi, keLo, keHi);
   outPng.ReplaceAll(":", "_");
   c->SaveAs(outPng);
   printf("saved %s\n", outPng.Data());
}
