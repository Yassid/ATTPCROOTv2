/// @file gate_view_sim.C
/// @brief What the production proton gate keeps and what it throws away, on simulated protons.
///
/// The gated acceptance collapses to 51 % and carries a dip at theta_cm 50-60 deg. This shows the
/// cause directly: the PID plane with the production gate on it, the protons that fall OUTSIDE it,
/// and the kinematics of exactly those protons -- so the loss is read in MeV and degrees rather
/// than inferred from a polygon.
///
/// Only truth-matched protons are drawn, so everything here is a real proton and every exclusion
/// is a loss, not a rejection of background.
///
///   root -b -q 'gate_view_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void gate_view_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                   TString fitSuffix = "_genfitter_pp", TString png = "plots/gate_view_sim.png",
                   Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   const double mp = 938.272;
   TString ukf = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1975/UKF/pid/";

   auto readGate = [](TString fn, std::vector<double> &ax, std::vector<double> &ay) {
      if (gSystem->AccessPathName(fn)) return;
      std::ifstream in(fn.Data());
      std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      size_t p = all.find("vertices");
      while (p != std::string::npos) {
         size_t a = all.find('[', p + 1);
         if (a == std::string::npos) break;
         size_t b = all.find(']', a);
         if (b == std::string::npos) break;
         double x, y; char c;
         std::istringstream is(all.substr(a + 1, b - a - 1));
         if (is >> x >> c >> y) { ax.push_back(x); ay.push_back(y); }
         p = b;
      }
   };
   std::vector<double> px, py, sx, sy;
   readGate(ukf + "proton_band.json", px, py);
   readGate(ukf + "proton_sim_6seeds.json", sx, sy);
   printf("  production gate: %zu vertices,  sim-built gate: %zu vertices\n", px.size(), sx.size());
   auto inside = [](double x, double y, const std::vector<double> &ax, const std::vector<double> &ay) {
      bool in = false; size_t n = ax.size();
      if (n < 3) return false;
      for (size_t i = 0, j = n - 1; i < n; j = i++)
         if (((ay[i] > y) != (ay[j] > y)) && (x < (ax[j] - ax[i]) * (y - ay[i]) / (ay[j] - ay[i]) + ax[i]))
            in = !in;
      return in;
   };

   auto *hAll = new TH2D("hAll", "truth-matched protons + gates;#sqrt{dE/dx};B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
   auto *hOut = new TH2D("hOut", "protons OUTSIDE the production gate;#sqrt{dE/dx};B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
   auto *kIn = new TH2D("kIn", "KEPT;#theta_{lab} [deg];KE [MeV]", 90, 0, 90, 90, 0, 45);
   auto *kOut = new TH2D("kOut", "LOST;#theta_{lab} [deg];KE [MeV]", 90, 0, 90, 90, 0, 45);
   long nIn = 0, nOut = 0, nSimIn = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fsim = dir + "/" + tg + "_sim.root", ffit = dir + "/" + tg + fitSuffix + ".root";
      if (gSystem->AccessPathName(fsim) || gSystem->AccessPathName(ffit)) continue;
      TFile *fs = TFile::Open(fsim), *ff = TFile::Open(ffit);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      TTree *tf = ff ? (TTree *)ff->Get("cbmsim") : nullptr;
      if (!ts || !tf || tf->GetEntries() > ts->GetEntries()) { if (fs) fs->Close(); if (ff) ff->Close(); continue; }
      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtPIDEvent", &pe);

      for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
         ts->GetEntry(i); tf->GetEntry(i);
         double keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212) continue;
            double px_ = p->GetPx() * 1000, py_ = p->GetPy() * 1000, pz_ = p->GetPz() * 1000;
            double pp = std::sqrt(px_ * px_ + py_ * py_ + pz_ * pz_);
            if (pp <= 0) continue;
            keT = std::sqrt(pp * pp + mp * mp) - mp;
            thT = std::acos(pz_ / pp) * TMath::RadToDeg();
            break;
         }
         if (keT < 0 || !pe || !pe->GetEntriesFast()) continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev) continue;
         double bd = 1e9, bx = 0, by = 0; bool got = false;
         for (auto &sp : ev->GetSpyral()) {
            if (!sp.valid) continue;
            double d = std::fabs(180.0 - sp.polar * TMath::RadToDeg() - thT);
            if (d < bd) { bd = d; bx = sp.sqrtdEdx; by = sp.brho; got = true; }
         }
         if (!got || bd > dThetaMax) continue;
         hAll->Fill(bx, by);
         if (inside(bx, by, sx, sy)) ++nSimIn;
         if (inside(bx, by, px, py)) { ++nIn; kIn->Fill(thT, keT); }
         else { ++nOut; hOut->Fill(bx, by); kOut->Fill(thT, keT); }
      }
      fs->Close(); ff->Close();
   }
   long tot = nIn + nOut;
   printf("\n  %ld truth-matched protons:  production gate keeps %ld (%.1f %%), loses %ld (%.1f %%)\n", tot, nIn,
          100.0 * nIn / tot, nOut, 100.0 * nOut / tot);
   printf("  the sim-built gate would keep %ld (%.1f %%)\n\n", nSimIn, 100.0 * nSimIn / tot);

   printf("  KE band     kept     lost    kept fraction\n");
   const double e[] = {0, 2, 5, 10, 20, 30, 45};
   for (int b = 0; b + 1 < (int)(sizeof(e) / sizeof(*e)); ++b) {
      double i1 = kIn->GetYaxis()->FindBin(e[b]), i2 = kIn->GetYaxis()->FindBin(e[b + 1]) - 1;
      double a = kIn->Integral(1, 90, i1, i2), c = kOut->Integral(1, 90, i1, i2);
      if (a + c < 1) continue;
      printf("  %4.0f-%-4.0f %8.0f %8.0f     %5.1f %%\n", e[b], e[b + 1], a, c, 100 * a / (a + c));
   }

   auto poly = [](const std::vector<double> &ax, const std::vector<double> &ay, int col, int sty) {
      auto *pl = new TPolyLine(ax.size() + 1);
      for (size_t i = 0; i < ax.size(); ++i) pl->SetPoint(i, ax[i], ay[i]);
      pl->SetPoint(ax.size(), ax[0], ay[0]);
      pl->SetLineColor(col); pl->SetLineWidth(3); pl->SetLineStyle(sty);
      return pl;
   };
   TCanvas *c = new TCanvas("cG", "gate", 1500, 900);
   c->Divide(2, 2);
   c->cd(1); gPad->SetLogz(); hAll->Draw("colz");
   poly(px, py, kBlack, 1)->Draw("L");
   if (sx.size() > 2) poly(sx, sy, kRed + 1, 2)->Draw("L");
   c->cd(2); gPad->SetLogz(); hOut->Draw("colz");
   poly(px, py, kBlack, 1)->Draw("L");
   c->cd(3); gPad->SetLogz(); kIn->Draw("colz");
   c->cd(4); gPad->SetLogz(); kOut->Draw("colz");
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("\n  wrote %s   (black = production proton_band.json, red dashed = sim-built)\n", png.Data());
}
