/// @file measure_gain_C15d.C
/// @brief Measure the per-run dE/dx gain factors from THIS analysis's own PID plane.
///
///   root -b -q 'measure_gain_C15d.C(0.30, 0.60)'          // Brho window, T*m
///   root -b -q 'measure_gain_C15d.C(0.30, 0.60, "/home/yassid/C15d_reco/", "gainmatch_C15d.csv")'
///
/// METHOD -- anchored on the dE/dx peak in a fixed Brho window.
///
/// Micromegas/GET gain drifts run to run, so the same particle deposits a different measured
/// charge in different runs. To separate that drift from physics, we need a quantity that
/// SHOULD be constant across runs and attribute any change in it to gain. Here that anchor is
/// the most-probable dE/dx of tracks in a fixed Brho slice:
///
///     sel(r)  = valid && brho in [bLo, bHi]
///     peak(r) = most-probable dEdx of sel(r)
///     f(r)    = peak_ref / peak(r)
///
/// The Brho window is what makes it species-anchored. Restricted to a slice where ONE species
/// dominates, the true dE/dx of the selected tracks is a property of the physics and does not
/// depend on how many of them there were -- so a run taken at a different rate, or with a
/// different mix of species overall, still yields the same peak unless the gain moved. That is
/// the whole reason not to anchor on the median of everything: a run with more deuterons has a
/// higher median for reasons that have nothing to do with gain, and that shift would be
/// absorbed into f(r) and silently distort the plane it was supposed to fix.
///
/// CHOOSE THE WINDOW BY LOOKING AT THE PLANE FIRST (mkpid_C15d.C). It must sit where one band
/// dominates; a window straddling two bands measures their ratio, not the gain, and the result
/// will look perfectly smooth while being wrong.
///
/// The reference is the POOLED peak over all runs (so the factors sit around 1 and no single
/// run's statistics define the scale), unless refRun is given. The overall normalisation is
/// arbitrary anyway -- only the run-to-run ratios matter -- but pooling avoids anchoring
/// everything to whichever run happened to be listed first.
///
/// Runs with too few tracks in the window get their factor INTERPOLATED from the neighbouring
/// runs in run-number order and are labelled `interp`; they are never silently assigned 1.0,
/// which would leave them unmatched inside a matched set and be invisible downstream.
///
/// Output is the CSV AtGainMatchTask reads: run,factor,source[,peak,n].

namespace {

/// Location estimator for the per-run dE/dx anchor.
///
/// ★ THE DEFAULT IS THE INTERQUARTILE TRUNCATED MEAN, NOT THE PEAK BIN, AND THAT MATTERS.
/// Measured on runs 17/19/20/21 (~790 tracks each, Brho 0.30-0.40):
///
///     estimator          run 17   run 19   run 20   run 21    spread
///     max bin (MPV)       303.6    308.3    284.8    341.3     20 %
///     median              305.1    302.9    301.4    309.4    2.6 %
///     trunc mean (IQR)    307.7    299.8    302.6    306.0    2.6 %
///
/// The true run-to-run drift over those adjacent runs is ~2.6 %. The 20 % from the peak bin is
/// pure estimator noise: the Landau peak is broad and flat, so on a single run's statistics the
/// highest bin wanders, and a 3-point parabola happily refines the wrong bin. Fed into f(run)
/// that noise becomes a 20 % per-run rescaling of dE/dx -- it would visibly smear the very bands
/// the gain match exists to sharpen, while looking like a perfectly reasonable measurement.
///
/// So: use an estimator that integrates many events. The interquartile truncated mean uses the
/// central 50 % of the sample, which is insensitive both to the Landau tail and to binning.
///
///   kind 0 = interquartile truncated mean (default)
///   kind 1 = median
///   kind 2 = max bin + 3-point parabola   -- kept for comparison; see above before using it
///   kind 3 = Landau fit                   -- proper MPV, but needs a good range and can fail quietly
double Anchor(TH1D *h, int kind, double &err)
{
   err = 0;
   const double N = h->GetEntries();
   if (N < 1)
      return -1;

   if (kind == 3) {
      const Int_t im = h->GetMaximumBin();
      const double lo = h->GetBinCenter(std::max(1, im - 20));
      const double hi = h->GetBinCenter(std::min(h->GetNbinsX(), im + 40));
      TF1 f("flan", "landau", lo, hi);
      f.SetParameters(h->GetMaximum(), h->GetBinCenter(im), 0.2 * h->GetBinCenter(im));
      if (h->Fit(&f, "QNR") == 0 && f.GetParameter(1) > 0) {
         err = f.GetParError(1);
         return f.GetParameter(1);
      }
      // fall through to the robust estimator rather than return a failed fit
   }

   if (kind == 2) {
      const Int_t im = h->GetMaximumBin();
      if (im <= 1 || im >= h->GetNbinsX())
         return h->GetBinCenter(im);
      const double y0 = h->GetBinContent(im - 1), y1 = h->GetBinContent(im), y2 = h->GetBinContent(im + 1);
      const double den = (y0 - 2 * y1 + y2);
      double shift = 0;
      if (std::abs(den) > 1e-12)
         shift = 0.5 * (y0 - y2) / den;
      if (std::abs(shift) > 1)
         shift = 0;
      const double w = h->GetBinWidth(im);
      const double n = y0 + y1 + y2;
      err = n > 0 ? w / std::sqrt(n) : w;
      return h->GetBinCenter(im) + shift * w;
   }

   double q[3];
   const double p[3] = {0.25, 0.5, 0.75};
   h->GetQuantiles(3, q, const_cast<double *>(p));
   if (kind == 1) {
      // Standard error of the median for a roughly normal core.
      err = 1.253 * (q[2] - q[0]) / 1.349 / std::sqrt(N);
      return q[1];
   }

   double s = 0, w = 0, s2 = 0;
   for (Int_t b = 1; b <= h->GetNbinsX(); ++b) {
      const double x = h->GetBinCenter(b);
      if (x < q[0] || x > q[2])
         continue;
      const double c = h->GetBinContent(b);
      s += x * c;
      s2 += x * x * c;
      w += c;
   }
   if (w <= 0)
      return q[1];
   const double mean = s / w;
   const double var = std::max(0.0, s2 / w - mean * mean);
   err = std::sqrt(var / w);
   return mean;
}

const char *kEstName[4] = {"interquartile truncated mean", "median", "max bin + parabola", "Landau fit"};

} // namespace

void measure_gain_C15d(Double_t bLo = 0.30, Double_t bHi = 0.60,
                       TString inDir = "/home/yassid/C15d_reco/", TString outCsv = "gainmatch_C15d.csv",
                       Int_t minTracks = 200, Int_t nbins = 200, Double_t dLo = 0.0, Double_t dHi = 0.0,
                       Int_t minClusters = 0, Int_t estimator = 0, Int_t refRun = -1,
                       TString plotDir = "plots/")
{
   gSystem->mkdir(plotDir, kTRUE);

   // ---- collect the per-run caches -------------------------------------------------------
   std::vector<TString> files;
   TSystemDirectory dir(inDir, inDir);
   TList *ls = dir.GetListOfFiles();
   if (ls == nullptr) {
      std::cout << "\033[1;31mERROR: cannot list " << inDir << "\033[0m\n";
      return;
   }
   TIter next(ls);
   while (auto *o = dynamic_cast<TSystemFile *>(next())) {
      TString n = o->GetName();
      if (!o->IsDirectory() && n.EndsWith("_pid.root"))
         files.push_back(inDir + n);
   }
   if (files.empty()) {
      std::cout << "\033[1;31mERROR: no *_pid.root in " << inDir << " -- run reco_batch.sh first.\033[0m\n";
      return;
   }
   std::sort(files.begin(), files.end());

   // ---- auto-range the dEdx axis from the pooled sample ------------------------------------
   // A hardcoded range would clip the peak of whichever runs sit at the extremes of the drift,
   // and a clipped peak reads as a perfectly good MPV at the edge bin.
   if (!(dHi > dLo)) {
      TChain pool("pid");
      for (const auto &f : files)
         pool.Add(f);
      TString sel = TString::Format("valid==1 && brho>%g && brho<%g", bLo, bHi);
      if (minClusters > 0)
         sel += TString::Format(" && nClusters>=%d", minClusters);
      pool.Draw("dEdx>>hauto(400,0,0)", sel, "goff");
      auto *ha = dynamic_cast<TH1 *>(gDirectory->Get("hauto"));
      if (ha == nullptr || ha->GetEntries() < 10) {
         std::cout << "\033[1;31mERROR: fewer than 10 tracks in Brho [" << bLo << "," << bHi
                   << "] across all runs. Wrong window?\033[0m\n";
         return;
      }
      // Range the axis on the MEDIAN, not on a high quantile. dE/dx is Landau-distributed with a
      // tail reaching several times the peak, so a 97 % quantile puts the upper edge ~4.5x the
      // peak: the bins then straddle the peak (+-12 on ~313, i.e. 4 %) and that coarseness lands
      // directly in the gain factor. 0 to 3x median keeps the whole peak with bins ~1.5 % of it.
      double q[1];
      const double p[1] = {0.5};
      ha->GetQuantiles(1, q, const_cast<double *>(p));
      dLo = 0.0;
      dHi = 3.0 * q[0];
   }

   std::cout << "\033[1;33m=== measure_gain_C15d ===\033[0m\n"
             << "  runs cached : " << files.size() << "\n"
             << "  anchor      : most-probable dEdx, Brho in [" << bLo << ", " << bHi << "] T*m\n"
             << "  dEdx range  : [" << dLo << ", " << dHi << "] in " << nbins << " bins\n"
             << "  estimator   : " << kEstName[estimator < 0 || estimator > 3 ? 0 : estimator] << "\n"
             << "  min tracks  : " << minTracks << " per run\n\n";

   // ---- per-run peak -----------------------------------------------------------------------
   struct RunPeak {
      int run;
      double peak;
      double err;
      long n;
      bool ok;
   };
   std::vector<RunPeak> peaks;

   for (const auto &fpath : files) {
      TFile *f = TFile::Open(fpath);
      if (f == nullptr || f->IsZombie())
         continue;
      auto *t = dynamic_cast<TTree *>(f->Get("pid"));
      if (t == nullptr) {
         f->Close();
         continue;
      }
      Int_t run = -1, valid = 0, nClusters = 0;
      Double_t dEdx = 0, brho = 0;
      t->SetBranchAddress("run", &run);
      t->SetBranchAddress("valid", &valid);
      t->SetBranchAddress("nClusters", &nClusters);
      t->SetBranchAddress("dEdx", &dEdx);
      t->SetBranchAddress("brho", &brho);

      TH1D h("hrun", ";dEdx;tracks", nbins, dLo, dHi);
      long n = 0;
      int theRun = -1;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         theRun = run;
         if (valid != 1)
            continue;
         if (minClusters > 0 && nClusters < minClusters)
            continue;
         if (brho <= bLo || brho >= bHi)
            continue;
         h.Fill(dEdx);
         ++n;
      }
      double err = 0;
      const double pk = (n >= minTracks) ? Anchor(&h, estimator, err) : -1;
      peaks.push_back({theRun, pk, err, n, n >= minTracks && pk > 0});
      f->Close();
   }
   std::sort(peaks.begin(), peaks.end(), [](const RunPeak &a, const RunPeak &b) { return a.run < b.run; });

   // ---- reference --------------------------------------------------------------------------
   double ref = -1;
   if (refRun > 0) {
      for (const auto &p : peaks)
         if (p.run == refRun && p.ok)
            ref = p.peak;
      if (ref <= 0) {
         std::cout << "\033[1;31mERROR: reference run " << refRun << " has no usable peak.\033[0m\n";
         return;
      }
   } else {
      // Pooled peak over every usable run, weighted by statistics.
      double sw = 0, swx = 0;
      for (const auto &p : peaks)
         if (p.ok) {
            sw += p.n;
            swx += p.n * p.peak;
         }
      if (sw <= 0) {
         std::cout << "\033[1;31mERROR: no run has >= " << minTracks << " tracks in the window.\033[0m\n";
         return;
      }
      ref = swx / sw;
   }
   std::cout << "  reference peak : " << ref << (refRun > 0 ? Form(" (run %d)", refRun) : " (pooled)")
             << "\n\n";

   // ---- factors, with interpolation for thin runs -------------------------------------------
   std::map<int, double> factor;
   std::map<int, TString> source;
   for (const auto &p : peaks)
      if (p.ok) {
         factor[p.run] = ref / p.peak;
         source[p.run] = "measured";
      }

   int nMeas = static_cast<int>(factor.size());
   for (const auto &p : peaks) {
      if (p.ok)
         continue;
      // Linear interpolation in run number between the nearest measured runs on each side.
      int lo = -1, hi = -1;
      for (const auto &kv : factor) {
         if (kv.first < p.run)
            lo = kv.first;
         else if (hi < 0 && kv.first > p.run)
            hi = kv.first;
      }
      double f;
      TString src;
      if (lo > 0 && hi > 0) {
         const double t = double(p.run - lo) / double(hi - lo);
         f = factor[lo] * (1 - t) + factor[hi] * t;
         src = "interp";
      } else if (lo > 0) {
         f = factor[lo];
         src = "held"; // no measured run after it: hold the last value rather than extrapolate
      } else if (hi > 0) {
         f = factor[hi];
         src = "held";
      } else {
         continue;
      }
      factor[p.run] = f;
      source[p.run] = src;
   }

   // ---- write -------------------------------------------------------------------------------
   std::ofstream out(outCsv.Data());
   out << "# C15d dE/dx gain-match factors, measured from this analysis's own PID plane\n";
   out << "# anchor      : most-probable dEdx for valid tracks with brho in [" << bLo << "," << bHi
       << "] T*m\n";
   out << "# reference   : " << ref << (refRun > 0 ? Form(" (run %d)", refRun) : " (statistics-weighted pooled peak)")
       << "\n";
   out << "# estimator   : " << kEstName[estimator < 0 || estimator > 3 ? 0 : estimator] << "\n";
   out << "# min tracks  : " << minTracks << " per run; thinner runs are interpolated in run number\n";
   out << "# apply       : dEdx *= factor, sqrt_dEdx *= sqrt(factor)\n";
   out << "run,factor,source,peak,ntracks\n";

   auto *gF = new TGraphErrors();
   auto *gP = new TGraphErrors();
   int ip = 0;
   double fmin = 1e30, fmax = -1e30;
   for (const auto &p : peaks) {
      auto it = factor.find(p.run);
      if (it == factor.end())
         continue;
      out << p.run << "," << it->second << "," << source[p.run] << "," << (p.ok ? p.peak : -1) << "," << p.n
          << "\n";
      fmin = std::min(fmin, it->second);
      fmax = std::max(fmax, it->second);
      if (p.ok) {
         gP->SetPoint(ip, p.run, p.peak);
         gP->SetPointError(ip, 0, p.err);
         gF->SetPoint(ip, p.run, it->second);
         gF->SetPointError(ip, 0, p.err * it->second / p.peak);
         ++ip;
      }
   }
   out.close();

   int nInterp = 0;
   for (const auto &kv : source)
      if (kv.second != "measured")
         ++nInterp;

   std::cout << "  \033[1;32mmeasured " << nMeas << " runs, " << nInterp << " interpolated/held\033[0m\n"
             << "  factor range   : " << fmin << " - " << fmax << "  (spread " << (fmin > 0 ? fmax / fmin : 0)
             << "x)\n"
             << "  wrote          : " << outCsv << "\n";

   TCanvas *c = new TCanvas("cgain", "gain match", 1100, 800);
   c->Divide(1, 2);
   c->cd(1);
   gP->SetTitle("anchor: most-probable dE/dx per run;run;MPV dEdx");
   gP->SetMarkerStyle(20);
   gP->Draw("AP");
   c->cd(2);
   gF->SetTitle("gain-match factor;run;factor");
   gF->SetMarkerStyle(20);
   gF->Draw("AP");
   c->SaveAs(plotDir + "gainmatch_C15d.png");
   std::cout << "  plot           : " << plotDir << "gainmatch_C15d.png\n";
}
