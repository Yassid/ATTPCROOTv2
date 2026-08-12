/// @file acceptance_fit_sim.C
/// @brief The a1975 16C(p,p') acceptance from the FITTED tracks, with and without the PID gate.
///
/// TWO NUMBERS, AND THE SECOND IS THE ONE THAT CORRECTS DATA. A fit-only acceptance describes a
/// production that fits every pattern track; the real (p,p) production fits only what the proton
/// gate admits, so a gated track that was never fitted has no Ex and cannot be recovered offline.
/// Both are reported: "fit" is reco + fit, "fit+gate" adds the gate the data production used.
///
/// The gate is matched to the fitted track by trackID, which is the id AtPIDTask stamps on the
/// AtSpyralResult and the fitter stamps on the AtFittedTrack -- NOT the array index, which differs
/// from the id about 38 % of the time.
///
/// theta_cm uses the non-relativistic elastic recoil relation theta_cm = 180 - 2*theta_lab. At
/// 12 MeV/u the relativistic correction is small but it is an approximation, not an identity.
///
/// THE ACCEPTANCE IS A RATIO, so the generator's flat theta_CM sampling largely cancels -- but the
/// sample is still flat in theta_CM, so each bin carries the statistics of a flat distribution and
/// not of the physical cross section. Bin-by-bin errors are quoted for that reason.
///
///   root -b -q 'acceptance_fit_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void acceptance_fit_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                        TString fitSuffix = "_genfitter_pp", TString gateJson = "", TString png = "plots/acceptance_fit_sim.png",
                        Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   const double mp = 938.272;
   if (!gateJson.Length())
      gateJson = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1975/UKF/pid/proton_band.json";

   // ---- the gate the data production used ----
   std::vector<double> gx, gy;
   {
      std::ifstream in(gateJson.Data());
      std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      size_t p = all.find("vertices");
      while (p != std::string::npos) {
         size_t a = all.find('[', p + 1);
         if (a == std::string::npos) break;
         size_t b = all.find(']', a);
         if (b == std::string::npos) break;
         double x, y; char c;
         std::istringstream is(all.substr(a + 1, b - a - 1));
         if (is >> x >> c >> y) { gx.push_back(x); gy.push_back(y); }
         p = b;
      }
   }
   printf("  gate %s : %zu vertices%s\n", gateJson.Data(), gx.size(), gx.size() > 2 ? "" : "  <-- NOT USABLE");
   auto inside = [&](double x, double y) {
      bool in = false;
      size_t n = gx.size();
      if (n < 3) return false;
      for (size_t i = 0, j = n - 1; i < n; j = i++)
         if (((gy[i] > y) != (gy[j] > y)) && (x < (gx[j] - gx[i]) * (y - gy[i]) / (gy[j] - gy[i]) + gx[i]))
            in = !in;
      return in;
   };

   const int NB = 18;
   auto *gen = new TH1D("gen", "", NB, 0, 90);
   auto *fit = new TH1D("fit", "", NB, 0, 90);
   auto *fgt = new TH1D("fgt", "", NB, 0, 90);
   // The SAME curve on a theta_cm axis, via the theta_cm = 180 - 2*theta_lab of the header. 18 bins
   // of 10 deg in theta_cm are the 18 bins of 5 deg in theta_lab, one for one and reversed, so this
   // is a relabelling of the axis and not a second measurement -- no event moves between bins.
   auto *genC = new TH1D("genC", "", NB, 0, 180);
   auto *fitC = new TH1D("fitC", "", NB, 0, 180);
   auto *fgtC = new TH1D("fgtC", "", NB, 0, 180);
   auto *genK = new TH1D("genK", "", 45, 0, 45);
   auto *fitK = new TH1D("fitK", "", 45, 0, 45);
   auto *fgtK = new TH1D("fgtK", "", 45, 0, 45);
   long nG = 0, nF = 0, nFG = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fsim = dir + "/" + tg + "_sim.root", ffit = dir + "/" + tg + fitSuffix + ".root";
      if (gSystem->AccessPathName(fsim) || gSystem->AccessPathName(ffit)) { printf("  skip %s\n", tg.Data()); continue; }
      TFile *fs = TFile::Open(fsim), *ff = TFile::Open(ffit);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      TTree *tf = ff ? (TTree *)ff->Get("cbmsim") : nullptr;
      if (!ts || !tf || tf->GetEntries() > ts->GetEntries()) { if (fs) fs->Close(); if (ff) ff->Close(); continue; }
      TClonesArray *mc = nullptr, *te = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtTrackingEvent", &te);
      tf->SetBranchAddress("AtPIDEvent", &pe);

      for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
         ts->GetEntry(i); tf->GetEntry(i);
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0) continue;
            keT = std::sqrt(pp * pp + mp * mp) - mp;
            thT = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (keT < 0) continue;
         ++nG; gen->Fill(thT); genK->Fill(keT); genC->Fill(180 - 2 * thT);

         // PID landscape of this event, keyed by trackID (NOT array index)
         std::map<int, bool> pass;
         if (pe && pe->GetEntriesFast()) {
            auto *ev = (AtPIDEvent *)pe->At(0);
            if (ev)
               for (auto &sp : ev->GetSpyral())
                  pass[sp.trackID] = sp.valid && inside(sp.sqrtdEdx, sp.brho);
         }

         double bd = 1e9; int bID = -1; bool got = false;
         if (te && te->GetEntriesFast()) {
            auto *ev = (AtTrackingEvent *)te->At(0);
            if (ev)
               for (auto &ft : ev->GetFittedTracks()) {
                  if (!ft) continue;
                  double th = ft->GetKinematics().theta * TMath::RadToDeg();
                  if (std::fabs(th - thT) < bd) { bd = std::fabs(th - thT); bID = ft->GetTrackID(); got = true; }
               }
         }
         if (got && bd < dThetaMax) {
            ++nF; fit->Fill(thT); fitK->Fill(keT); fitC->Fill(180 - 2 * thT);
            auto g = pass.find(bID);
            if (g != pass.end() && g->second) {
               ++nFG; fgt->Fill(thT); fgtK->Fill(keT); fgtC->Fill(180 - 2 * thT);
            }
         }
      }
      fs->Close(); ff->Close();
      printf("  %-8s done\n", tg.Data());
   }

   printf("\n  generated %ld,  fitted %ld (%.1f %%),  fitted+gated %ld (%.1f %%)\n\n", nG, nF, 100.0 * nF / nG, nFG,
          100.0 * nFG / nG);
   printf("  theta_lab    theta_cm      gen      fit acceptance      fit+gate acceptance\n");
   for (int b = 1; b <= NB; ++b) {
      double g = gen->GetBinContent(b);
      if (g < 1) continue;
      double lo = gen->GetBinLowEdge(b), hi = lo + gen->GetBinWidth(b);
      double f = fit->GetBinContent(b), q = fgt->GetBinContent(b);
      double ef = f / g, eq = q / g;
      printf("  %3.0f-%-3.0f     %4.0f-%-4.0f  %7.0f    %5.1f +- %3.1f %%        %5.1f +- %3.1f %%\n", lo, hi,
             180 - 2 * hi, 180 - 2 * lo, g, 100 * ef, 100 * std::sqrt(ef * (1 - ef) / g), 100 * eq,
             100 * std::sqrt(std::max(0., eq * (1 - eq)) / g));
   }
   printf("\n  KE band      gen      fit acceptance      fit+gate acceptance\n");
   const double edge[] = {0, 1, 2, 3, 5, 10, 20, 45};
   for (int b = 0; b + 1 < (int)(sizeof(edge) / sizeof(*edge)); ++b) {
      double g = genK->Integral(genK->FindBin(edge[b]), genK->FindBin(edge[b + 1]) - 1);
      double f = fitK->Integral(fitK->FindBin(edge[b]), fitK->FindBin(edge[b + 1]) - 1);
      double q = fgtK->Integral(fgtK->FindBin(edge[b]), fgtK->FindBin(edge[b + 1]) - 1);
      if (g < 1) continue;
      printf("  %5.0f-%-5.0f %7.0f    %5.1f %%             %5.1f %%\n", edge[b], edge[b + 1], g, 100 * f / g,
             100 * q / g);
   }

   auto mk = [](TH1D *n, TH1D *d, int col) {
      auto *e = (TH1D *)n->Clone(Form("%s_e", n->GetName()));
      e->Divide(n, d, 1, 1, "B");
      e->SetLineColor(col); e->SetLineWidth(2); e->SetMinimum(0); e->SetMaximum(1.05);
      return e;
   };
   TCanvas *c = new TCanvas("cA", "acceptance", 1200, 500);
   c->Divide(2, 1);
   c->cd(1);
   auto *eA = mk(fit, gen, kRed + 1);
   eA->SetTitle("acceptance vs #theta_{lab}: fit (red), fit+PID gate (blue);#theta_{lab} [deg];accepted / generated");
   eA->Draw("e"); mk(fgt, gen, kAzure + 2)->Draw("e same");
   c->cd(2);
   auto *eK = mk(fitK, genK, kRed + 1);
   eK->SetTitle("acceptance vs KE;KE [MeV];accepted / generated");
   eK->Draw("e"); mk(fgtK, genK, kAzure + 2)->Draw("e same");
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);

   // same curve, theta_cm axis
   TString pngCm = png;
   pngCm.ReplaceAll(".png", "_thcm.png");
   TCanvas *c2 = new TCanvas("cC", "acceptance vs theta_cm", 900, 650);
   gPad->SetGridy();
   auto *eCf = mk(fitC, genC, kRed + 1);
   eCf->SetTitle("^{16}C(p,p) acceptance;#theta_{cm} [deg];accepted / generated");
   eCf->SetMarkerStyle(24); eCf->SetMarkerColor(kRed + 1);
   eCf->Draw("e");
   auto *eCg = mk(fgtC, genC, kAzure + 2);
   eCg->SetMarkerStyle(20); eCg->SetMarkerColor(kAzure + 2);
   eCg->Draw("e same");
   auto *lg = new TLegend(0.15, 0.15, 0.55, 0.30);
   lg->AddEntry(eCf, "reco + fit", "lp");
   lg->AddEntry(eCg, "reco + fit + PID gate", "lp");
   lg->Draw();
   c2->SaveAs(pngCm);

   TString rootCm = pngCm; rootCm.ReplaceAll(".png", ".root");
   TFile fo(rootCm, "RECREATE");
   genC->Write(); fitC->Write(); fgtC->Write();
   eCf->Write("hAcc_fit"); eCg->Write("hAcc_gate");
   fo.Close();

   printf("\n  theta_cm       gen      fit acceptance     fit+gate acceptance\n");
   for (int b = 1; b <= NB; ++b) {
      double g = genC->GetBinContent(b);
      if (g < 1) continue;
      printf("  %3.0f-%-3.0f  %8.0f   %5.1f +- %3.1f %%       %5.1f +- %3.1f %%\n", genC->GetBinLowEdge(b),
             genC->GetBinLowEdge(b) + genC->GetBinWidth(b), g, 100 * eCf->GetBinContent(b),
             100 * eCf->GetBinError(b), 100 * eCg->GetBinContent(b), 100 * eCg->GetBinError(b));
   }
   printf("\n  wrote %s\n         %s\n         %s\n", png.Data(), pngCm.Data(), rootCm.Data());
}
