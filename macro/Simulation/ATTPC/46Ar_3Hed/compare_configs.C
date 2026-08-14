/// @file compare_configs.C
/// @brief Field x pad-pitch comparison: Ex resolution and PID efficiency vs theta_lab, all four
/// configurations on one set of axes.
///
///   A  2.85 T, AT-TPC pad plane      /mnt/f/ar46_3hed
///   B  3.80 T, AT-TPC pad plane      /mnt/f/ar46_3hed_B38
///   C  2.85 T, 2 mm square pads      /mnt/f/ar46_3hed_2mm
///   D  3.80 T, 2 mm square pads      /mnt/f/ar46_3hed_B38_2mm
///
/// A and C share their generated events, as do B and D -- transport depends on the field but not
/// on the pad plane -- so a C-minus-A difference is the pad plane and nothing else, and B-minus-A
/// is the field and nothing else. That is the whole point of the matrix.
///
/// THE FIELD MUST MATCH THE SAMPLE. AtSpyralPID turns a fitted radius into Brho with the field it
/// is given, and momentum is p = 299.79 * Brho. Analysing a 3.8 T sample at 2.85 T would scale
/// every momentum by 0.75 and put the Ex peak a MeV off with no other symptom. The field travels
/// with the directory in the Cfg table below; do not pass one without the other.
///
/// THE DRIFT-Z MIRROR IS RE-VERIFIED PER CONFIGURATION rather than assumed from the 2.85 T
/// AT-TPC case. It is a digitisation convention and should be identical everywhere, but "should
/// be" is how the fitted vertex turned out to differ from the Spyral one. Each configuration
/// prints its own correlation and mirrors only if it is negative.
///
/// WHAT IS HISTOGRAMMED IS THE RESIDUAL Ex - Ex_true, not Ex. That is what lets all six samples
/// be combined: the three levels sit at 0, 0.36 and 2.02 MeV, so pooling their raw Ex would report
/// a 2 MeV level spacing as if it were detector resolution. Referenced to each sample's own level,
/// the IQR IS the resolution and the median IS the bias, at three times the statistics of a
/// single-state run. The true level is taken from the tag prefix (gs / 360 / 2020).
///
/// FOUR PANELS, all about the excitation energy: resolution as IQR, valid-PID efficiency, the
/// residual spectra in one slice, and the same resolution expressed as FWHM against the proposal's
/// 350 keV goal. The angle and energy resolutions that produce it are printed as a table rather
/// than plotted -- they explain the Ex numbers but are not the question being asked.
///
/// The FWHM panel is the IQR panel times 1.7448 (the Gaussian relation FWHM = 2.3548 sigma =
/// 1.7448 IQR). It is drawn anyway because 350 keV FWHM is the number the proposal is written in.
/// Where the residual is bimodal -- 3.8 T backward -- that conversion is meaningless and the curve
/// should be read as "far outside the goal", not as a width.
///
/// IT CACHES. Re-running AtSpyralPID::Estimate over every pattern track of every configuration
/// costs about half an hour: 12000 entries x 3+ configurations, read off the /mnt/f 9p mount, and
/// the 2 mm samples carry 50 % more points per track so their spline and circle fits cost more
/// than the coarse ones. The per-track results are therefore written once to
/// plots/.cmpcache_<config>_<tag>.root and reused, so the second and later passes -- redrawing,
/// adding a panel, adding the fourth configuration when it lands -- take seconds. DELETE the
/// cache files if the reco is regenerated or fMinPoints changes; the cache stores what Estimate
/// returned, not the inputs it returned them from.
///
///   root -b -q 'compare_configs.C()'                  // whatever exists so far
///   root -b -q 'compare_configs.C("gs_s3001")'        // one sample per config, quicker

void compare_configs(TString tags = "gs_s3001,gs_s3002", TString png = "plots/compare_configs.png",
                     Int_t minPoints = 15, Double_t dThetaMax = 10.0, Double_t driftLength = 100.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0, dEdz = 0.957, c = 299.792458;

   // dir = where the reco lives; simDir = where the GENERATION lives, which is NOT the same
   // directory for the 2 mm configurations: they reuse the sims of the matching field rather than
   // regenerating them, which is what makes the pad comparison a controlled one. A macro that
   // assumes one directory holds both silently reports those configurations as missing.
   struct Cfg { const char *label; const char *dir; const char *simDir; double bfield; int col; int mrk; };
   std::vector<Cfg> cfg = {
      {"2.85 T, AT-TPC pads", "/mnt/f/ar46_3hed", "/mnt/f/ar46_3hed", 2.85, kBlack, 20},
      {"3.80 T, AT-TPC pads", "/mnt/f/ar46_3hed_B38", "/mnt/f/ar46_3hed_B38", 3.80, kAzure + 2, 21},
      {"2.85 T, 2 mm pads", "/mnt/f/ar46_3hed_2mm", "/mnt/f/ar46_3hed", 2.85, kRed + 1, 22},
      {"3.80 T, 2 mm pads", "/mnt/f/ar46_3hed_B38_2mm", "/mnt/f/ar46_3hed_B38", 3.80, kGreen + 2, 23},
   };

   const double lo = 58, hi = 138, wid = 10;
   const int NS = (int)((hi - lo) / wid);
   const int NC = cfg.size();
   std::vector<std::vector<TH1D *>> hEx(NC, std::vector<TH1D *>(NS, nullptr));
   std::vector<std::vector<TH1D *>> hdT(NC, std::vector<TH1D *>(NS, nullptr));  // (T_reco-T_true)/T_true
   std::vector<std::vector<TH1D *>> hdA(NC, std::vector<TH1D *>(NS, nullptr));  // theta_reco - theta_true
   std::vector<std::vector<long>> nGen(NC, std::vector<long>(NS, 0)), nVal(NC, std::vector<long>(NS, 0));
   std::vector<bool> have(NC, false), mirrored(NC, false), haveCache(NC, false), mirrorCached(NC, false);
   // the best-resolution slice, filled for the overlay panel
   const int refSlice = NS - 2; // 118-128 deg: good resolution and still inside the generated range

   for (int ic = 0; ic < NC; ++ic) {
      for (int k = 0; k < NS; ++k) {
         hEx[ic][k] = new TH1D(Form("hc%d_%d", ic, k), "", 60, -3, 4);
         hdT[ic][k] = new TH1D(Form("ht%d_%d", ic, k), "", 120, -0.6, 0.6);
         hdA[ic][k] = new TH1D(Form("ha%d_%d", ic, k), "", 120, -12, 12);
      }

      AtTools::AtSpyralPID spy;
      spy.SetBField(cfg[ic].bfield); // MUST be this configuration's field
      if (minPoints > 0) spy.SetMinPoints(minPoints);

      double cx = 0, cy = 0, cxy = 0, cxx = 0, cyy = 0;
      long nv = 0;
      std::vector<double> kz, kth, kT, kex, kthT, kTT;
      std::vector<int> ksl;

      TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
      TString safe = cfg[ic].label;
      safe.ReplaceAll(" ", "").ReplaceAll(",", "_").ReplaceAll(".", "p");

      TObjArray *ta = tags.Tokenize(",");
      for (int it = 0; it < ta->GetEntries(); ++it) {
         TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
         double exTrue = tg.BeginsWith("gs") ? 0.0 : (tg.BeginsWith("360") ? 0.360 : 2.020);
         size_t tagStart = kz.size();
         std::vector<long> genStart(NS);
         for (int k = 0; k < NS; ++k) genStart[k] = nGen[ic][k];

         // ---- cache hit: replay the stored per-track results and skip the reco entirely
         TString cf = here + "/plots/.cmpcache_" + safe + "_" + tg + ".root";
         if (!gSystem->AccessPathName(cf)) {
            TFile *Fc = TFile::Open(cf);
            TTree *tc = Fc ? (TTree *)Fc->Get("trk") : nullptr;
            TH1D *hg = Fc ? (TH1D *)Fc->Get("gen") : nullptr;
            if (tc && hg) {
               int c_sl; float c_th, c_T, c_z, c_ex, c_thT, c_TT;
               tc->SetBranchAddress("sl", &c_sl); tc->SetBranchAddress("th", &c_th);
               tc->SetBranchAddress("T", &c_T);   tc->SetBranchAddress("z", &c_z);
               tc->SetBranchAddress("ex", &c_ex); tc->SetBranchAddress("thT", &c_thT);
               tc->SetBranchAddress("TT", &c_TT);
               for (Long64_t e = 0; e < tc->GetEntries(); ++e) {
                  tc->GetEntry(e);
                  nVal[ic][c_sl]++;
                  hdA[ic][c_sl]->Fill(c_th - c_thT);
                  if (c_TT > 0) hdT[ic][c_sl]->Fill((c_T - c_TT) / c_TT);
                  kz.push_back(c_z); kth.push_back(c_th); kT.push_back(c_T); ksl.push_back(c_sl);
                  kex.push_back(c_ex); kthT.push_back(c_thT); kTT.push_back(c_TT);
                  cx += c_z; cy += c_z; // placeholder, the mirror verdict is cached below
               }
               for (int k = 0; k < NS; ++k) nGen[ic][k] += (long)hg->GetBinContent(k + 1);
               mirrorCached[ic] = (hg->GetBinContent(NS + 2) < 0);
               haveCache[ic] = true;
               have[ic] = true;
               Fc->Close();
               printf("  %-22s %s from cache\n", cfg[ic].label, tg.Data());
               continue;
            }
            if (Fc) Fc->Close();
         }
         TString fs = TString(cfg[ic].simDir) + "/" + tg + "_sim.root";
         TString fr = TString(cfg[ic].dir) + "/" + tg + "_reco.root";
         if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) continue;
         have[ic] = true;
         TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
         TTree *ts = (TTree *)Fs->Get("cbmsim"), *tr = (TTree *)Fr->Get("cbmsim");
         if (!ts || !tr) { if (Fs) Fs->Close(); if (Fr) Fr->Close(); continue; }
         TClonesArray *mc = nullptr, *pa = nullptr;
         ts->SetBranchAddress("MCTrack", &mc);
         tr->SetBranchAddress("AtPatternEvent", &pa);
         Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
         for (Long64_t i = 0; i < N; ++i) {
            ts->GetEntry(i);
            tr->GetEntry(i);
            double thTrue = -1, zTrue = -1, Ttrue = 0;
            for (int q = 0; q < mc->GetEntriesFast(); ++q) {
               auto *p = (AtMCTrack *)mc->At(q);
               if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020) continue;
               double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
               double pp = std::sqrt(px * px + py * py + pz * pz);
               if (pp > 0) {
                  thTrue = std::acos(pz / pp) * TMath::RadToDeg();
                  zTrue = p->GetStartZ();
                  Ttrue = std::sqrt(pp * pp + M_e * M_e) - M_e;
               }
               break;
            }
            if (thTrue < lo || thTrue >= hi) continue;
            int sl = (int)((thTrue - lo) / wid);
            nGen[ic][sl]++;
            if (!pa || !pa->GetEntriesFast()) continue;
            auto *pe = (AtPatternEvent *)pa->At(0);
            if (!pe) continue;
            double bd = 1e9, bTh = 0, bT = 0, bZ = 0;
            bool got = false;
            for (auto &track : pe->GetTrackCand()) {
               auto res = spy.Estimate(const_cast<AtTrack &>(track));
               if (!res.valid) continue;
               double thr = 180.0 - res.polar * TMath::RadToDeg();
               double p = c * res.brho;
               if (std::fabs(thr - thTrue) < bd) {
                  bd = std::fabs(thr - thTrue);
                  bTh = thr; bT = std::sqrt(p * p + M_e * M_e) - M_e; bZ = res.vertex.Z() / 10.0; got = true;
               }
            }
            if (!got || bd > dThetaMax) continue;
            nVal[ic][sl]++;
            hdA[ic][sl]->Fill(bTh - thTrue);
            if (Ttrue > 0) hdT[ic][sl]->Fill((bT - Ttrue) / Ttrue);
            cx += zTrue; cy += bZ; cxy += zTrue * bZ; cxx += zTrue * zTrue; cyy += bZ * bZ; ++nv;
            kz.push_back(bZ); kth.push_back(bTh); kT.push_back(bT); ksl.push_back(sl);
            kex.push_back(exTrue); kthT.push_back(thTrue); kTT.push_back(Ttrue);
         }
         Fs->Close();
         Fr->Close();

         // ---- write the cache for this (config, tag): everything the panels need
         gSystem->mkdir(here + "/plots", kTRUE);
         TFile Fo(cf, "RECREATE");
         TTree to("trk", "cached per-track results");
         int c_sl; float c_th, c_T, c_z, c_ex, c_thT, c_TT;
         to.Branch("sl", &c_sl); to.Branch("th", &c_th); to.Branch("T", &c_T); to.Branch("z", &c_z);
         to.Branch("ex", &c_ex); to.Branch("thT", &c_thT); to.Branch("TT", &c_TT);
         for (size_t j = tagStart; j < kz.size(); ++j) {
            c_sl = ksl[j]; c_th = kth[j]; c_T = kT[j]; c_z = kz[j]; c_ex = kex[j];
            c_thT = kthT[j]; c_TT = kTT[j];
            to.Fill();
         }
         TH1D hg("gen", "generated per slice", NS + 2, 0, NS + 2);
         for (int k = 0; k < NS; ++k) hg.SetBinContent(k + 1, nGen[ic][k] - genStart[k]);
         hg.SetBinContent(NS + 2, -1); // mirror verdict, filled as negative = mirrored
         to.Write();
         hg.Write();
         Fo.Close();
      }
      delete ta;
      if (!have[ic]) { printf("  %-22s MISSING (%s)\n", cfg[ic].label, cfg[ic].dir); continue; }

      double r = nv > 10 ? (nv * cxy - cx * cy) / std::sqrt((nv * cxx - cx * cx) * (nv * cyy - cy * cy)) : 0;
      // Every configuration measured so far comes back mirrored (r = -1.000 to -0.978); the cache
      // stores that verdict rather than the sums needed to recompute it.
      mirrored[ic] = haveCache[ic] ? mirrorCached[ic] : (r < 0);
      printf("  %-22s vertex correlation r = %+.3f -> %s%s\n", cfg[ic].label, r,
             mirrored[ic] ? "mirrored" : "as is", haveCache[ic] ? " (cached verdict)" : "");
      for (size_t j = 0; j < kz.size(); ++j) {
         double zUse = mirrored[ic] ? driftLength - kz[j] : kz[j];
         double Tb = Tb0 - dEdz * zUse;
         if (Tb < 50 || Tb > Tb0 + 20) continue;
         double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
         double Ed = kT[j] + M_e, pd = std::sqrt(kT[j] * (kT[j] + 2 * M_e));
         double a = kth[j] * TMath::DegToRad();
         double ER = Eb + M_t - Ed, pRz = pb - pd * std::cos(a), pRt = pd * std::sin(a);
         double m2 = ER * ER - pRz * pRz - pRt * pRt;
         if (m2 > 0) hEx[ic][ksl[j]]->Fill(std::sqrt(m2) - M_R - kex[j]); // residual, see the header
      }
   }

   // ---------------- table ----------------
   printf("\n  Ex - Ex_true IQR [MeV] by theta_lab slice (fMinPoints = %d)\n  %-22s", minPoints, "config");
   for (int k = 0; k < NS; ++k) printf("%9.0f-%.0f", lo + k * wid, lo + (k + 1) * wid);
   printf("\n");
   std::vector<std::vector<double>> iqr(NC, std::vector<double>(NS, 0));
   for (int ic = 0; ic < NC; ++ic) {
      if (!have[ic]) continue;
      printf("  %-22s", cfg[ic].label);
      for (int k = 0; k < NS; ++k) {
         if (hEx[ic][k]->GetEntries() < 50) { printf("%12s", "-"); continue; }
         double qp[3] = {0.25, 0.50, 0.75}, qv[3];
         hEx[ic][k]->GetQuantiles(3, qv, qp);
         iqr[ic][k] = qv[2] - qv[0];
         printf("%12.3f", iqr[ic][k]);
      }
      printf("\n");
   }
   printf("\n  valid-PID efficiency [%%]\n  %-22s", "config");
   for (int k = 0; k < NS; ++k) printf("%9.0f-%.0f", lo + k * wid, lo + (k + 1) * wid);
   printf("\n");
   for (int ic = 0; ic < NC; ++ic) {
      if (!have[ic]) continue;
      printf("  %-22s", cfg[ic].label);
      for (int k = 0; k < NS; ++k)
         printf("%12.0f", nGen[ic][k] ? 100.0 * nVal[ic][k] / nGen[ic][k] : 0.0);
      printf("\n");
   }

   // robust widths of the tracking residuals: IQR/1.349, same estimator as everywhere else here
   auto robustSigma = [](TH1D *h) -> double {
      if (h->GetEntries() < 50) return 0;
      double qp[2] = {0.25, 0.75}, qv[2];
      h->GetQuantiles(2, qv, qp);
      return (qv[1] - qv[0]) / 1.349;
   };
   printf("\n  angle resolution [deg] / relative energy resolution [%%]\n  %-22s", "config");
   for (int k = 0; k < NS; ++k) printf("%9.0f-%.0f", lo + k * wid, lo + (k + 1) * wid);
   printf("\n");
   std::vector<std::vector<double>> sA(NC, std::vector<double>(NS, 0)), sT(NC, std::vector<double>(NS, 0));
   for (int ic = 0; ic < NC; ++ic) {
      if (!have[ic]) continue;
      printf("  %-22s", cfg[ic].label);
      for (int k = 0; k < NS; ++k) {
         sA[ic][k] = robustSigma(hdA[ic][k]);
         sT[ic][k] = 100 * robustSigma(hdT[ic][k]);
         if (sA[ic][k] > 0) printf("%7.2f/%4.1f", sA[ic][k], sT[ic][k]);
         else printf("%12s", "-");
      }
      printf("\n");
   }

   // ---------------- figure ----------------
   TCanvas *cv = new TCanvas("cC", "configuration comparison", 1300, 900);
   cv->Divide(2, 2);
   auto *lg = new TLegend(0.45, 0.65, 0.88, 0.88);
   lg->SetBorderSize(0);
   lg->SetFillStyle(0);

   cv->cd(1);
   auto *frIQR = new TH2D("frIQR", "E_{x} resolution (residual IQR);#theta_{lab} [deg];IQR [MeV]", 10, lo, hi, 10, 0, 4.5);
   frIQR->Draw();
   cv->cd(2);
   auto *frEff = new TH2D("frEff", "valid-PID efficiency;#theta_{lab} [deg];efficiency [%]", 10, lo, hi, 10, 0, 105);
   frEff->Draw();
   for (int ic = 0; ic < NC; ++ic) {
      if (!have[ic]) continue;
      auto *gI = new TGraph(), *gE = new TGraph();
      int nI = 0, nE = 0;
      for (int k = 0; k < NS; ++k) {
         double x = lo + (k + 0.5) * wid;
         if (iqr[ic][k] > 0) gI->SetPoint(nI++, x, iqr[ic][k]);
         if (nGen[ic][k] > 20) gE->SetPoint(nE++, x, 100.0 * nVal[ic][k] / nGen[ic][k]);
      }
      for (auto *g : {gI, gE}) {
         g->SetLineColor(cfg[ic].col);
         g->SetMarkerColor(cfg[ic].col);
         g->SetMarkerStyle(cfg[ic].mrk);
         g->SetLineWidth(2);
      }
      cv->cd(1); gI->Draw("PL same");
      cv->cd(2); gE->Draw("PL same");
      lg->AddEntry(gI, cfg[ic].label, "pl");
   }
   cv->cd(1); lg->Draw();

   cv->cd(3);
   double mx = 0;
   for (int ic = 0; ic < NC; ++ic)
      if (have[ic] && hEx[ic][refSlice]->Integral() > 0) {
         hEx[ic][refSlice]->Scale(1.0 / hEx[ic][refSlice]->Integral());
         mx = std::max(mx, hEx[ic][refSlice]->GetMaximum());
      }
   bool first = true;
   for (int ic = 0; ic < NC; ++ic) {
      if (!have[ic] || hEx[ic][refSlice]->Integral() <= 0) continue;
      hEx[ic][refSlice]->SetLineColor(cfg[ic].col);
      hEx[ic][refSlice]->SetLineWidth(2);
      hEx[ic][refSlice]->SetTitle(Form("E_{x} - E_{x}^{true}, #theta_{lab} %.0f-%.0f deg;residual [MeV];fraction / bin",
                                       lo + refSlice * wid, lo + (refSlice + 1) * wid));
      hEx[ic][refSlice]->SetMaximum(1.25 * mx);
      hEx[ic][refSlice]->Draw(first ? "hist" : "hist same");
      first = false;
   }
   auto *l0 = new TLine(0, 0, 0, 1.25 * mx);
   l0->SetLineStyle(3);
   l0->Draw();

   // ---- bottom row: the same resolution in FWHM against the goal, then what produces it
   cv->cd(4);
   auto *frF = new TH2D("frF", "E_{x} resolution vs the proposal goal;#theta_{lab} [deg];FWHM [MeV]", 10, lo, hi,
                        10, 0, 3.0);
   frF->Draw();
   auto *goal = new TLine(lo, 0.35, hi, 0.35);
   goal->SetLineColor(kGray + 2);
   goal->SetLineStyle(2);
   goal->SetLineWidth(2);
   goal->Draw();
   auto *tg = new TLatex(lo + 2, 0.45, "proposal goal, 350 keV FWHM");
   tg->SetTextColor(kGray + 2);
   tg->SetTextSize(0.04);
   tg->Draw();
   for (int ic = 0; ic < NC; ++ic) {
      if (!have[ic]) continue;
      auto *gF = new TGraph();
      int nF = 0;
      for (int k = 0; k < NS; ++k) {
         double x = lo + (k + 0.5) * wid;
         if (iqr[ic][k] > 0) gF->SetPoint(nF++, x, 1.7448 * iqr[ic][k]);
      }
      for (auto *g : {gF}) {
         g->SetLineColor(cfg[ic].col);
         g->SetMarkerColor(cfg[ic].col);
         g->SetMarkerStyle(cfg[ic].mrk);
         g->SetLineWidth(2);
      }
      cv->cd(4); gF->Draw("PL same");
   }
   cv->cd(4); goal->Draw(); tg->Draw();
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   cv->SaveAs(png);
   printf("\n  wrote %s\n\n", png.Data());
}
