/// @file kinematics_sim.C
/// @brief Generated vs reconstructed kinematics for the a1975 16C(p,p') simulation -- what is lost.
///
/// TRUTH comes from MCTrack: the primary proton, KE and theta_lab from its momentum.
///
/// RECONSTRUCTED comes from AtPIDEvent, NOT from a fit -- the accumulation deliberately runs no
/// fitter. AtSpyralPID gives a rigidity from the circle through the hit cloud, so
///     p [MeV/c] = 299.792458 * Z * Brho [T.m]        (Z = 1)
/// and theta_lab = 180 - polar, which is the convention AtSpyralPID reports in.
/// The same proton mass is used on both sides so the comparison carries no offset of its own.
///
/// EFFICIENCY is truth-matched by angle within dThetaMax, exactly as pid_gate_from_sim.C does:
/// matching on rigidity as well would hide the badly reconstructed tracks, which are the ones the
/// question is about.
///
///   root -b -q 'kinematics_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void kinematics_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                    TString png = "plots/kinematics_sim.png", Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double mp = 938.272; // MeV

   auto *hT = new TH2D("hT", "GENERATED (MC truth);#theta_{lab} [deg];KE [MeV]", 90, 0, 90, 90, 0, 45);
   auto *hR = new TH2D("hR", "RECONSTRUCTED (Spyral B#rho);#theta_{lab} [deg];KE [MeV]", 90, 0, 90, 90, 0, 45);
   auto *keT = new TH1D("keT", "KE;KE [MeV];tracks", 90, 0, 45);
   auto *keR = new TH1D("keR", "KE", 90, 0, 45);
   auto *keM = new TH1D("keM", "KE matched", 90, 0, 45);
   auto *thT = new TH1D("thT", "#theta_{lab};#theta_{lab} [deg];tracks", 90, 0, 90);
   auto *thR = new TH1D("thR", "#theta", 90, 0, 90);
   auto *thM = new TH1D("thM", "#theta matched", 90, 0, 90);
   long nTruth = 0, nReco = 0, nMatch = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fsim = dir + "/" + tg + "_sim.root", fpid = dir + "/" + tg + "_pid.root";
      if (gSystem->AccessPathName(fsim) || gSystem->AccessPathName(fpid)) {
         printf("  skip %s (missing)\n", tg.Data());
         continue;
      }
      TFile *fs = TFile::Open(fsim), *fp = TFile::Open(fpid);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      TTree *tp = fp ? (TTree *)fp->Get("cbmsim") : nullptr;
      if (!ts || !tp || tp->GetEntries() > ts->GetEntries()) {
         printf("  %s: entry mismatch, skipping\n", tg.Data());
         if (fs) fs->Close();
         if (fp) fp->Close();
         continue;
      }
      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tp->SetBranchAddress("AtPIDEvent", &pe);

      for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
         ts->GetEntry(i);
         tp->GetEntry(i);

         // ---- truth: the primary proton -------------------------------------------------
         double keTrue = -1, thTrue = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0)
               continue;
            keTrue = std::sqrt(pp * pp + mp * mp) - mp;
            thTrue = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (keTrue >= 0) {
            ++nTruth;
            hT->Fill(thTrue, keTrue);
            keT->Fill(keTrue);
            thT->Fill(thTrue);
         }

         // ---- reconstructed: every valid Spyral entry -------------------------------------
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         double bd = 1e9;
         bool got = false;
         for (auto &sp : ev->GetSpyral()) {
            if (!sp.valid)
               continue;
            double pR = 299.792458 * sp.brho;
            double keRec = std::sqrt(pR * pR + mp * mp) - mp;
            double thRec = 180.0 - sp.polar * TMath::RadToDeg();
            ++nReco;
            hR->Fill(thRec, keRec);
            keR->Fill(keRec);
            thR->Fill(thRec);
            if (keTrue >= 0 && std::fabs(thRec - thTrue) < bd) {
               bd = std::fabs(thRec - thTrue);
               got = true;
            }
         }
         if (got && bd < dThetaMax && keTrue >= 0) {
            ++nMatch;
            keM->Fill(keTrue);
            thM->Fill(thTrue);
         }
      }
      fs->Close();
      fp->Close();
      printf("  %-8s done\n", tg.Data());
   }
   printf("\n  %ld generated protons, %ld reconstructed PID entries, %ld truth-matched (%.1f %%)\n\n", nTruth, nReco,
          nMatch, nTruth ? 100.0 * nMatch / nTruth : 0.0);

   // where the loss sits, in the units the acceptance argument is made in
   printf("  KE band      generated   matched    efficiency\n");
   const double edge[] = {0, 1, 2, 3, 5, 10, 20, 45};
   for (int b = 0; b + 1 < (int)(sizeof(edge) / sizeof(*edge)); ++b) {
      double g = keT->Integral(keT->FindBin(edge[b]), keT->FindBin(edge[b + 1]) - 1);
      double m = keM->Integral(keM->FindBin(edge[b]), keM->FindBin(edge[b + 1]) - 1);
      printf("  %5.0f-%-5.0f  %9.0f %9.0f      %5.1f %%\n", edge[b], edge[b + 1], g, m, g ? 100.0 * m / g : 0.0);
   }

   // Is the acceptance FLAT across the plateau, or does it carry structure of its own? This is the
   // question that decides whether a minimum in the measured angular distribution is physics. The
   // recoil relation for elastic scattering is theta_cm = 180 - 2*theta_lab, so a feature reported
   // at theta_cm 55-60 deg would sit at theta_lab 60-62.5 deg -- inside the plateau, where any dip
   // in the acceptance would masquerade as one in the cross section.
   printf("\n  theta_lab   theta_cm    generated   matched   efficiency\n");
   for (int b = 0; b < 18; ++b) {
      double lo = b * 5.0, hi = lo + 5.0;
      double g = thT->Integral(thT->FindBin(lo), thT->FindBin(hi) - 1);
      double m = thM->Integral(thM->FindBin(lo), thM->FindBin(hi) - 1);
      if (g <= 0)
         continue;
      double e = m / g, se = std::sqrt(e * (1 - e) / g);
      printf("  %3.0f-%-3.0f    %4.0f-%-4.0f  %9.0f %9.0f    %5.1f +- %.1f %%\n", lo, hi, 180 - 2 * hi, 180 - 2 * lo,
             g, m, 100 * e, 100 * se);
   }
   {  // flatness of the plateau, stated as a number rather than left to the eye
      double s = 0, s2 = 0, mn = 1e9, mx = -1e9;
      int n = 0;
      for (int i = thT->FindBin(25.0); i <= thT->FindBin(75.0); ++i) {
         double g = thT->GetBinContent(i);
         if (g < 50)
            continue;
         double e = thM->GetBinContent(i) / g;
         s += e; s2 += e * e; ++n;
         mn = std::min(mn, e); mx = std::max(mx, e);
      }
      if (n) {
         double mean = s / n, rms = std::sqrt(std::max(0.0, s2 / n - mean * mean));
         printf("\n  plateau 25-75 deg: mean %.3f, rms %.3f (%.1f %%), min %.3f, max %.3f over %d bins\n", mean, rms,
                100 * rms / mean, mn, mx, n);
      }
   }

   auto eff = [](TH1D *num, TH1D *den, const char *nm, const char *ti) {
      auto *e = (TH1D *)num->Clone(nm);
      e->Divide(num, den, 1, 1, "B");
      e->SetTitle(ti);
      e->SetMinimum(0);
      e->SetMaximum(1.05);
      e->SetLineColor(kRed + 1);
      e->SetLineWidth(2);
      return e;
   };

   TCanvas *c = new TCanvas("cK", "kinematics", 1500, 900);
   c->Divide(3, 2);
   c->cd(1); gPad->SetLogz(); hT->Draw("colz");
   c->cd(2); gPad->SetLogz(); hR->Draw("colz");
   c->cd(3);
   keT->SetLineColor(kBlack); keT->SetLineWidth(2);
   keR->SetLineColor(kRed + 1); keR->SetLineWidth(2);
   keT->SetTitle("KE: generated (black) vs reconstructed (red)");
   keT->Draw("hist"); keR->Draw("hist same");
   c->cd(4);
   thT->SetLineColor(kBlack); thT->SetLineWidth(2);
   thR->SetLineColor(kRed + 1); thR->SetLineWidth(2);
   thT->SetTitle("#theta_{lab}: generated (black) vs reconstructed (red)");
   thT->Draw("hist"); thR->Draw("hist same");
   c->cd(5);
   eff(keM, keT, "effKE", "efficiency vs generated KE;KE [MeV];matched / generated")->Draw("hist");
   c->cd(6);
   eff(thM, thT, "effTh", "efficiency vs generated #theta_{lab};#theta_{lab} [deg];matched / generated")->Draw("hist");

   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("\n  wrote %s\n", png.Data());
}
