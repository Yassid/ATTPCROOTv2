/// @file mc_vs_reco_sim.C
/// @brief MC truth against reconstruction for 16C(p,p), with NO fitting anywhere.
///
/// Truth from <tag>_sim.root, reconstruction from <tag>_reco.root. AtSpyralPID::Estimate() is run
/// here in memory, exactly as AtPIDTask::Exec drives it, rather than read back from *_pid.root --
/// because fMinPoints has to be settable and the persisted files were written at the class default
/// of 30. Nothing here opens a genfitter file.
///
/// WHY fMinPoints MATTERS. Measured on these same six seeds, it is the only cut that rejects
/// anything: at 30 it accounts for 3451 of the 3453 rejections (99.94 %). Relaxing it recovers
/// short -- i.e. low-energy -- tracks, which is the population the forward-theta_cm acceptance is
/// missing:
///     30 -> 89.9 %   25 -> 92.4 %   20 -> 94.6 %   15 -> 96.7 %   10 -> 98.3 %   5 -> 98.8 %
/// Below 15 the direction checks (codes 8 and 9) start to fire, so the useful knee is 15-20.
/// NEITHER AtPIDTask NOR AtGenfitter EVER CALLS SetMinPoints, so the data runs at 30 as well --
/// changing this is not a simulation-only knob.
///
/// HANDEDNESS. The reconstructed polar comes out as 180 - theta_true (the simulation reverses
/// drift z in digitisation). Matching sp.polar directly against truth gives exactly 0 % in every
/// bin. The 2D panel plots the raw polar so this is visible rather than folded away.
///
/// WHAT "found" AND "valid" MEAN. AtSpyralPID assigns polar, radius, direction, nPoints and valid
/// in one block at the end of Estimate(), so an early exit leaves all of them at defaults and a
/// rejected track can never match a truth angle. "found" is therefore a lower limit on track
/// finding, not a measurement of it; the rejected count is printed alongside so the gap is visible.
///
///   root -b -q 'mc_vs_reco_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006",15)'

static double mr_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double mr_thetacm(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a4 = (std::cos(thetalab) * mr_omega2(s, m1 * m1, m2 * m2) * mr_omega2(u, m2 * m2, m3 * m3) -
                (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                  (2 * m2 * m2) +
               s + u - m2 * m2;
   if (a4 < 0) return std::nan("");
   double m4x = std::sqrt(a4);
   double t = m2 * m2 + m4x * m4x - 2 * m2 * Et4;
   double c = (s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4x * m4x) +
               (m1 * m1 - m2 * m2) * (m3 * m3 - m4x * m4x)) /
              (mr_omega2(s, m1 * m1, m2 * m2) * mr_omega2(s, m3 * m3, m4x * m4x));
   if (c < -1 || c > 1) return std::nan("");
   return (TMath::Pi() - std::acos(c)) * TMath::RadToDeg();
}

/// @param recoDir  where the <tag>_reco.root live, if not alongside the truth. Empty means the same
///                 directory as dir. The HDBSCAN comparison keeps its reco in its own directory
///                 while the truth stays with the TC production, so the two have to be separable.
void mc_vs_reco_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                    Int_t minPoints = 15, Double_t Ebeam = 188.0, Double_t bField = 2.85,
                    Double_t dThetaMax = 10.0, Int_t nBins = 18, TString tag = "mp15", TString recoDir = "",
                    TString gateJson = "")
{
   if (!recoDir.Length()) recoDir = dir;
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double u = 931.49401, m_C16 = 16.0147 * u, m_p = 1.007825 * u;

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0) spy.SetMinPoints(minPoints);
   printf("\n  fMinPoints = %d   (class default 30; data runs the default too)\n", minPoints > 0 ? minPoints : 30);

   // optional PID gate, applied to the SAME in-memory estimate the acceptance is measured from,
   // so the gate and the points it judges come from one fMinPoints rather than two.
   std::vector<double> gx, gy;
   if (gateJson.Length()) {
      std::ifstream in(gateJson.Data());
      std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      size_t p = all.find("vertices");
      while (p != std::string::npos) {
         size_t a = all.find('[', p + 1); if (a == std::string::npos) break;
         size_t b = all.find(']', a);     if (b == std::string::npos) break;
         double vx, vy; char cc;
         std::istringstream is(all.substr(a + 1, b - a - 1));
         if (is >> vx >> cc >> vy) { gx.push_back(vx); gy.push_back(vy); }
         p = b;
      }
      printf("  gate %s : %zu vertices%s\n", gateJson.Data(), gx.size(), gx.size() > 2 ? "" : "  <-- UNUSABLE");
      if (gx.size() < 3) { printf("\033[1;31mrefusing to run: an unusable gate silently means NO gate\033[0m\n"); return; }
   }
   auto inGate = [&](double qx, double qy) {
      if (gx.size() < 3) return true;
      bool in = false; size_t n = gx.size();
      for (size_t i = 0, j = n - 1; i < n; j = i++)
         if (((gy[i] > qy) != (gy[j] > qy)) && (qx < (gx[j] - gx[i]) * (qy - gy[i]) / (gy[j] - gy[i]) + gx[i]))
            in = !in;
      return in;
   };

   auto *hGen = new TH1D("hGen", "", nBins, 0, 180);
   auto *hFnd = new TH1D("hFnd", "", nBins, 0, 180);
   auto *hVal = new TH1D("hVal", "", nBins, 0, 180);
   hGen->Sumw2(); hFnd->Sumw2(); hVal->Sumw2();
   auto *h2 = new TH2D("h2", "16C(p,p) protons, no fit;#theta_{lab} generated [deg];#theta_{lab} reconstructed [deg]",
                       180, 0, 180, 180, 0, 180);
   long nG = 0, nF = 0, nV = 0, nTrk = 0, nRej = 0;
   double sum = 0, sum2 = 0; long nRes = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fs = dir + "/" + tg + "_sim.root", fr = recoDir + "/" + tg + "_reco.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) { printf("  skip %s\n", tg.Data()); continue; }
      TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
      TTree *ts = Fs ? (TTree *)Fs->Get("cbmsim") : nullptr;
      TTree *tr = Fr ? (TTree *)Fr->Get("cbmsim") : nullptr;
      if (!ts || !tr) { if (Fs) Fs->Close(); if (Fr) Fr->Close(); continue; }
      TClonesArray *mc = nullptr, *pa = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pa);

      Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i); tr->GetEntry(i);
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0) continue;
            keT = std::sqrt(pp * pp + m_p * m_p) - m_p;
            thT = std::acos(pz / pp);
            break;
         }
         if (keT <= 0) continue;
         double cm = mr_thetacm(m_C16, m_p, m_p, m_C16, Ebeam, thT, keT);
         if (std::isnan(cm)) continue;
         ++nG; hGen->Fill(cm);
         double thTdeg = thT * TMath::RadToDeg();

         bool found = false, valid = false;
         if (pa && pa->GetEntriesFast()) {
            auto *pe = (AtPatternEvent *)pa->At(0);
            if (pe)
               for (auto &track : pe->GetTrackCand()) {
                  AtTrack &t2 = const_cast<AtTrack &>(track);
                  auto res = spy.Estimate(t2);
                  ++nTrk;
                  if (res.polar <= 0) { ++nRej; continue; }
                  double thR = res.polar * TMath::RadToDeg();
                  h2->Fill(thTdeg, thR);
                  double d = (180.0 - thR) - thTdeg;
                  if (std::fabs(d) > dThetaMax) continue;
                  found = true;
                  sum += d; sum2 += d * d; ++nRes;
                  if (res.valid && inGate(res.sqrtdEdx, res.brho)) valid = true;
               }
         }
         if (found) { ++nF; hFnd->Fill(cm); }
         if (valid) { ++nV; hVal->Fill(cm); }
      }
      Fs->Close(); Fr->Close();
      printf("  %-8s done\n", tg.Data());
   }
   delete ta;
   if (!nG) { printf("\033[1;31mno truth protons found\033[0m\n"); return; }

   double mean = nRes ? sum / nRes : 0, rms = nRes ? std::sqrt(sum2 / nRes - mean * mean) : 0;
   printf("\n  generated protons          %ld\n", nG);
   printf("  pattern tracks estimated   %ld   (rejected by Spyral: %ld, %.1f %%)\n", nTrk, nRej,
          100.0 * nRej / std::max(1L, nTrk));
   printf("  matched to truth           %ld  (%.1f %%)\n", nF, 100.0 * nF / nG);
   printf("  ... and valid Spyral PID   %ld  (%.1f %%)\n", nV, 100.0 * nV / nG);
   printf("  angle residual about 180-x : mean %.2f deg, rms %.2f deg  (%ld entries)\n\n", mean, rms, nRes);

   auto mk = [&](TH1D *n, const char *nm, int col, int mrk) {
      auto *e = (TH1D *)n->Clone(nm);
      e->Divide(n, hGen, 1, 1, "B");
      e->SetDirectory(nullptr);
      e->SetTitle("^{16}C(p,p) reconstruction, no fit;#theta_{cm} [deg];reconstructed / generated");
      e->SetMinimum(0); e->SetMaximum(1.05);
      e->SetLineColor(col); e->SetMarkerColor(col); e->SetMarkerStyle(mrk); e->SetLineWidth(2);
      return e;
   };
   TH1D *aF = mk(hFnd, "hAcc_found", kGreen + 3, 24);
   TH1D *aV = mk(hVal, "hAcc_valid", kAzure + 2, 20);

   printf("  theta_cm      gen     matched          valid PID\n");
   for (int b = 1; b <= nBins; ++b) {
      double g = hGen->GetBinContent(b);
      if (g < 1) continue;
      printf("  %3.0f-%-3.0f %8.0f   %5.1f +- %3.1f %%    %5.1f +- %3.1f %%\n", hGen->GetBinLowEdge(b),
             hGen->GetBinLowEdge(b) + hGen->GetBinWidth(b), g, 100 * aF->GetBinContent(b), 100 * aF->GetBinError(b),
             100 * aV->GetBinContent(b), 100 * aV->GetBinError(b));
   }

   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TCanvas *c = new TCanvas("cMR", "MC vs reco", 1400, 640);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetGridy();
   aF->Draw("E1"); aV->Draw("E1 same");
   auto *lg = new TLegend(0.15, 0.15, 0.60, 0.30);
   lg->AddEntry(aF, Form("matched to truth (minPts %d)", minPoints > 0 ? minPoints : 30), "lp");
   lg->AddEntry(aV, gateJson.Length() ? "+ valid PID + gate" : "+ valid Spyral PID", "lp");
   lg->Draw();
   c->cd(2);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   h2->Draw("colz");
   auto *l = new TLine(0, 180, 180, 0);
   l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->SetLineStyle(2); l->Draw();
   auto *lid = new TLine(0, 0, 180, 180);
   lid->SetLineColor(kGray + 2); lid->SetLineWidth(1); lid->SetLineStyle(3); lid->Draw();

   TString png = "plots/mc_vs_reco_" + tag + ".png";
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(here + "/" + png);
   TFile fo(here + "/plots/mc_vs_reco_" + tag + ".root", "RECREATE");
   hGen->Write(); hFnd->Write(); hVal->Write(); h2->Write(); aF->Write(); aV->Write();
   fo.Close();
   printf("\n  wrote %s\n         plots/mc_vs_reco_%s.root\n\n", png.Data(), tag.Data());
}
