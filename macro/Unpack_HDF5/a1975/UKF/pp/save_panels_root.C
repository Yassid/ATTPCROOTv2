/// @file save_panels_root.C
/// @brief Write the browser explorer's panels into a ROOT file, so they can be sent to someone
///        who wants to re-fit, restyle or re-cut them in plain ROOT.
///
/// The explorer's own "save figure / save panels" buttons produce PNGs and its "save data" gives
/// JSON -- fine for a report, useless for someone who wants to change the fit. This writes the
/// actual objects:
///
///   pk           TTree   the SELECTED events (ke, theta, vertexz, chi2ndf, ex, thcm) -- everything
///                        else can be rebuilt from this, including with different cuts
///   hEx          TH1D    excitation energy
///   hKEth        TH2D    ejectile KE vs theta_lab
///   hExThcm      TH2D    Ex vs theta_cm
///   hExVz        TH2D    Ex vs vertex z
///   locus_Ex*    TGraph  the kinematic locus of each reference level, for the KE-vs-theta panel
///   fGS          TF1     gaussian+linear fit to the reference peak (only if it converged)
///   config       TNamed  every number that produced the file, in the title
///   cPanels      TCanvas the 2x2 arrangement, drawn from the objects above
///
/// The TTree is the important one: a colleague can change Ebeam, cuts or binning and redo the lot,
/// which a canvas alone would not allow. Ex/thcm are stored as computed here AND recomputable,
/// since ke/theta are kept.
///
/// Channel-generic -- masses default to 16C(p,p')16C:
///   root -b -q 'pp/save_panels_root.C("cache.root","pp_panels.root","16C(p,p'\'')",195.5)'
///   (p,d): mEject 2.0135532, mResid 15.0105993, refEx "0.740,3.103,4.220,4.657"
///   (d,d): mTarg = mEject = 2.0135532
///
/// Reading it back needs no ATTPCROOT at all -- plain ROOT:  TFile f("pp_panels.root"); f.ls();

static double om_(double x, double y, double z) { return sqrt(x*x + y*y + z*z - 2*x*y - 2*y*z - 2*x*z); }

struct KinPair { double ex, thcm; };
static KinPair kine_(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3, E4 = E1 + m2 - E3;
   double s = m1*m1 + m2*m2 + 2*m2*E1, u = m2*m2 + m3*m3 - 2*m2*E3;
   double a = (cos(th)*om_(s,m1*m1,m2*m2)*om_(u,m2*m2,m3*m3) - (s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)
              + s + u - m2*m2;
   if (a < 0) return {std::nan(""), std::nan("")};
   double m4x = sqrt(a), ex = m4x - m4, t = m2*m2 + m4x*m4x - 2*m2*E4;
   double g = (s*s + s*(2*t - m1*m1 - m2*m2 - m3*m3 - m4x*m4x) + (m1*m1 - m2*m2)*(m3*m3 - m4x*m4x))
              / (om_(s,m1*m1,m2*m2)*om_(s,m3*m3,m4x*m4x));
   if (g < -1) g = -1; if (g > 1) g = 1;
   return {ex, (TMath::Pi() - acos(g))*TMath::RadToDeg()};
}

/// KE at which Ex(KE) = exTarget for this theta_lab (Ex falls monotonically with KE -> bisect)
static double keAtEx_(double m1, double m2, double m3, double m4, double E, double th, double exT)
{
   double lo = 0.05, hi = 200.0;
   double flo = kine_(m1,m2,m3,m4,E,th,lo).ex, fhi = kine_(m1,m2,m3,m4,E,th,hi).ex;
   if (std::isnan(flo) || std::isnan(fhi) || (flo-exT)*(fhi-exT) > 0) return std::nan("");
   for (int i = 0; i < 80; i++) {
      double mid = 0.5*(lo+hi), fm = kine_(m1,m2,m3,m4,E,th,mid).ex;
      if (std::isnan(fm)) { hi = mid; continue; }
      if ((flo-exT)*(fm-exT) <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
   }
   return 0.5*(lo+hi);
}

void save_panels_root(TString cache, TString outFile, TString tag = "16C(p,p')", double Ebeam = 195.5,
                      double mBeamAmu = 16.0147013, double mTargAmu = 1.00782503,
                      double mEjectAmu = 1.00782503, double mResidAmu = 16.0147013,
                      TString refExCSV = "1.766,3.027,3.986,4.142",
                      double chi2max = 5.0, double icMin = 950, double icMax = 1350,
                      double thLo = 0, double thHi = 180,
                      double exLo = -5, double exHi = 25, int exBins = 200,
                      double gsLo = -1.0, double gsHi = 1.0)
{
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const double m1 = mBeamAmu*u, m2 = mTargAmu*u, m3 = mEjectAmu*u, m4 = mResidAmu*u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree `pk` in %s\n", cache.Data()); return; }

   // caches differ: (p,p')/(p,d) use `vz` + `ic`; the D2 caches use `vertexz` and have no IC
   float ke, th, vz = 0, c2, ic = -1;
   t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&th); t->SetBranchAddress("chi2ndf",&c2);
   if (t->GetBranch("vz")) t->SetBranchAddress("vz",&vz);
   else if (t->GetBranch("vertexz")) t->SetBranchAddress("vertexz",&vz);
   const bool hasIC = t->GetBranch("ic") != nullptr;
   if (hasIC) t->SetBranchAddress("ic",&ic); else if (icMin > 0) { icMin = -1; }

   TFile out(outFile, "RECREATE");

   float oke, oth, ovz, oc2, oex, othcm;
   TTree *pk = new TTree("pk", Form("%s selected events (Ebeam=%.2f)", tag.Data(), Ebeam));
   pk->Branch("ke",&oke); pk->Branch("theta",&oth); pk->Branch("vertexz",&ovz);
   pk->Branch("chi2ndf",&oc2); pk->Branch("ex",&oex); pk->Branch("thcm",&othcm);

   auto *hEx = new TH1D("hEx", Form("%s;E_{x} [MeV];counts", tag.Data()), exBins, exLo, exHi);
   auto *hKEth = new TH2D("hKEth", Form("%s;#theta_{lab} [deg];ejectile KE [MeV]", tag.Data()),
                          100, 0, 95, 100, 0, 45);
   auto *hExThcm = new TH2D("hExThcm", Form("%s;#theta_{cm} [deg];E_{x} [MeV]", tag.Data()),
                            90, 0, 180, 100, exLo, exHi);
   auto *hExVz = new TH2D("hExVz", Form("%s;vertex z [mm];E_{x} [MeV]", tag.Data()),
                          100, -100, 1100, 100, exLo, exHi);

   Long64_t nSel = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ke <= 0 || c2 > chi2max) continue;
      if (icMin > 0 && (ic < icMin || ic > icMax)) continue;
      if (th < thLo || th >= thHi) continue;
      auto k = kine_(m1,m2,m3,m4,Ebeam,th*TMath::DegToRad(),ke);
      if (std::isnan(k.ex)) continue;
      oke = ke; oth = th; ovz = vz; oc2 = c2; oex = k.ex; othcm = k.thcm;
      pk->Fill(); ++nSel;
      hEx->Fill(k.ex); hKEth->Fill(th, ke); hExThcm->Fill(k.thcm, k.ex); hExVz->Fill(vz, k.ex);
   }
   printf("%s: %lld selected (chi2/ndf<%.1f%s)\n", tag.Data(), nSel, chi2max,
          icMin > 0 ? Form(", IC[%.0f,%.0f]", icMin, icMax) : ", no IC gate");

   // reference-peak fit, over the window the caller asked for
   TF1 *fGS = new TF1("fGS", "gaus(0)+pol1(3)", gsLo, gsHi);
   int bm = hEx->GetXaxis()->FindBin(0.5*(gsLo+gsHi));
   fGS->SetParameters(hEx->GetBinContent(bm), 0.5*(gsLo+gsHi), 0.3, 0, 0);
   fGS->SetParLimits(1, gsLo, gsHi); fGS->SetParLimits(2, 0.03, 2.0);
   bool ok = (hEx->Fit(fGS, "QRN") == 0);
   // A centroid sitting on the window edge is the minimiser running out of room, not a measurement.
   // Ship no fit at all rather than a number that looks real -- the caller must move the window.
   if (ok) {
      double mu = fGS->GetParameter(1), pad = 0.02*(gsHi - gsLo);
      if (mu < gsLo + pad || mu > gsHi - pad) {
         printf("  reference peak: PINNED at the window edge (mu=%+.3f in [%.2f,%.2f]) -- not written;"
                " move the window\n", mu, gsLo, gsHi);
         ok = false;
      }
   }
   if (ok)
      printf("  reference peak: mu = %+.3f +- %.3f, sigma = %.3f, FWHM = %.3f MeV\n",
             fGS->GetParameter(1), fGS->GetParError(1), fGS->GetParameter(2), 2.3548*fGS->GetParameter(2));
   else if (hEx->GetEntries() > 0)
      printf("  reference peak: no usable fit in [%.2f,%.2f]\n", gsLo, gsHi);

   // kinematic locus of each reference level, for overlaying on the KE-vs-theta panel
   std::vector<TGraph *> loci;
   TObjArray *toks = refExCSV.Tokenize(",");
   std::vector<double> refs{0.0};
   for (int i = 0; i < toks->GetEntries(); ++i)
      refs.push_back(((TObjString *)toks->At(i))->GetString().Atof());
   for (size_t r = 0; r < refs.size(); ++r) {
      auto *g = new TGraph();
      for (double a = 1; a < 95; a += 0.5) {
         double v = keAtEx_(m1,m2,m3,m4,Ebeam,a*TMath::DegToRad(),refs[r]);
         if (!std::isnan(v) && v > 0 && v < 200) g->SetPoint(g->GetN(), a, v);
      }
      g->SetName(Form("locus_Ex%.3f", refs[r]));
      g->SetTitle(Form("E_{x} = %.3f MeV locus", refs[r]));
      g->SetLineColor(r == 0 ? kRed+1 : kAzure+r); g->SetLineWidth(2);
      loci.push_back(g);
   }

   auto *cfg = new TNamed("config",
      Form("%s | Ebeam=%.3f MeV | masses(amu) beam=%.7f targ=%.7f eject=%.7f resid=%.7f | "
           "chi2/ndf<%.2f | IC[%.0f,%.0f] | theta[%.1f,%.1f] | refEx=%s | gsWindow[%.2f,%.2f] | "
           "cache=%s | selected=%lld",
           tag.Data(), Ebeam, mBeamAmu, mTargAmu, mEjectAmu, mResidAmu, chi2max, icMin, icMax,
           thLo, thHi, refExCSV.Data(), gsLo, gsHi, gSystem->BaseName(cache), nSel));

   auto *c = new TCanvas("cPanels", tag, 1400, 1000);
   c->Divide(2,2);
   c->cd(1); hEx->Draw("hist"); if (ok) { fGS->SetLineColor(kRed); fGS->Draw("same"); }
   c->cd(2); gPad->SetLogz(); hKEth->Draw("colz"); for (auto *g : loci) g->Draw("L same");
   c->cd(3); gPad->SetLogz(); hExThcm->Draw("colz");
   c->cd(4); gPad->SetLogz(); hExVz->Draw("colz");

   out.cd();
   pk->Write(); hEx->Write(); hKEth->Write(); hExThcm->Write(); hExVz->Write();
   for (auto *g : loci) g->Write();
   if (ok) fGS->Write();
   cfg->Write();
   c->Write();
   out.Close(); f->Close();
   printf("  wrote %s  (pk tree + 4 histos + %zu loci + config + cPanels)\n",
          outFile.Data(), loci.size());
}
