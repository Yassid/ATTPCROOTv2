/// @file closure_report.C
/// @brief Print the closure numbers for a set of proton_kin caches: g.s. centroid, width, and
///        the Ex-vs-theta_cm slope. Truth in the simulation is exactly Ebeam = 161.00 MeV, so
///        the centroid mu IS the reconstruction bias -- no reference spectrum needed.
///
///   root -b -q 'closure_report.C("negBgM,negBgP")'

void closure_report(TString tagsCSV = "negBgM,negBgP", TString plotDir = "", Double_t chi2Cut = 5.0)
{
   if (plotDir.IsNull())
      plotDir = TString(gSystem->Getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1954/UKF/pp/plots/";
   printf("\n%-10s %7s %8s %8s %8s %8s %11s\n", "tag", "N", "mu", "sigma", "FWHM", "|mu|", "dEx/dthcm");
   printf("%s\n", TString('-', 68).Data());

   TObjArray *tags = tagsCSV.Tokenize(",");
   for (int i = 0; i < tags->GetEntries(); ++i) {
      TString t = ((TObjString *)tags->At(i))->GetString().Strip(TString::kBoth);
      TString fn = plotDir + "proton_kin_" + t + ".root";
      TFile *f = TFile::Open(fn);
      if (!f || f->IsZombie()) {
         printf("%-10s  cannot open %s\n", t.Data(), fn.Data());
         continue;
      }
      TTree *tr = (TTree *)f->Get("pk");
      if (!tr) {
         printf("%-10s  no ntuple\n", t.Data());
         continue;
      }
      TString cut = TString::Format("chi2ndf<%g", chi2Cut);

      // g.s. peak: gaussian over |Ex| < 1 MeV, the same window the explorer uses
      // NOTE: do NOT detach the histogram before Draw. With SetDirectory(nullptr) the
      // ">>h" target cannot be found in gDirectory, so ROOT silently fills a DIFFERENT
      // histogram and the local pointer stays empty -- which read as N=0 for every config.
      tr->Draw("ex>>h(120,-3,3)", cut, "goff");
      TH1 *h = (TH1 *)gDirectory->Get("h");
      if (!h) { printf("%-10s  draw failed\n", t.Data()); f->Close(); continue; }
      double mu = 0, sg = 0;
      Long64_t n = (Long64_t)h->GetEntries();
      if (h->GetEntries() > 20) {
         TF1 g("g", "gaus", -1, 1);
         g.SetParameters(h->GetMaximum(), h->GetBinCenter(h->GetMaximumBin()), 0.3);
         h->Fit(&g, "QRN");
         mu = g.GetParameter(1);
         sg = std::fabs(g.GetParameter(2));
      }

      // trend: least-squares slope of Ex vs theta_cm over 15-45 deg, on the g.s. band only.
      // (Locating the Ebeam that minimises the swing is noise-dominated below ~20k tracks;
      // a slope at fixed Ebeam is well determined with far less.)
      // "prof" is REQUIRED: without it ">>p" builds a TH2, and casting that to TProfile*
      // gives meaningless GetBinEntries() -> the slope silently came out 0.00000.
      tr->Draw("ex:thcm>>p(30,15,45)", cut + "&&fabs(ex)<1", "prof goff");
      TProfile *p = (TProfile *)gDirectory->Get("p");
      if (!p) { printf("%-10s  profile failed\n", t.Data()); f->Close(); continue; }
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      int m = 0;
      for (int b = 1; b <= p->GetNbinsX(); ++b) {
         if (p->GetBinEntries(b) < 5)
            continue;
         double xx = p->GetBinCenter(b), yy = p->GetBinContent(b);
         sx += xx; sy += yy; sxx += xx * xx; sxy += xx * yy; ++m;
      }
      double slope = (m > 2 && (m * sxx - sx * sx) != 0) ? (m * sxy - sx * sy) / (m * sxx - sx * sx) : 0;

      printf("%-10s %7lld %8.3f %8.3f %8.3f %8.3f %11.5f\n", t.Data(), n, mu, sg, 2.355 * sg, std::fabs(mu), slope);
      f->Close();
   }
   printf("\ntruth is exactly 161.00 MeV -> the correct configuration is the one with mu ~ 0\n\n");
}
