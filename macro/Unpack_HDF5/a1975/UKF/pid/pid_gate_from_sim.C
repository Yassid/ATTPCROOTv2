/// @file pid_gate_from_sim.C
/// @brief Draw the a1975 proton PID gate from simulated protons instead of by hand.
///
/// WHY NOT BY HAND. The gate in use (pid/proton_band.json) was drawn on the data plane, where the
/// proton band overlaps deuterons and the beam and its edges are a judgement call. The cutflow on a
/// truth-matched simulation showed what that costs: of 5999 generated protons, 92 % get a PID
/// entry, 87 % a valid one, but only 54 % fall inside the gate. The 38 % it throws away are not
/// random -- it keeps 94 % of the protons below 5 MeV and cuts hardest on the fast ones, clipping
/// the upper rigidity edge around sqrt(dE/dx) 4-8. That is exactly the region the elastic angular
/// distribution needs.
///
/// The simulation has no such ambiguity: every primary is a proton, so a truth-matched PID entry IS
/// a proton and the band it traces is the answer, edges included.
///
/// THE PREMISE, CHECKED FIRST. A gate built on the simulation is only usable on the data if the two
/// bands sit in the same place. Comparing the sqrt(dE/dx) ridge slice by slice in rigidity:
///
///     Brho     0.19   0.25   0.31   0.37   0.43   0.49   0.55   0.61   0.67   0.73   0.79   0.85
///     data/sim 1.04   0.98   0.94   0.90   0.89   1.00   0.98   1.07   1.00   1.00   0.97   0.91
///
/// agreement to 5-10 % with no trend, so the dE/dx scales match and the transfer is legitimate.
///
/// HOW THE BAND IS BUILT. In each slice of sqrt(dE/dx), the rigidities of the truth-matched protons
/// are ordered and the interval from the (1-keep)/2 to the (1+keep)/2 quantile is taken. keep is
/// deliberately high (0.98 by default): this gate is a particle identification, not a quality cut,
/// and every proton it drops has to be put back by the acceptance correction later. Slices with
/// fewer than nMin entries are skipped rather than extrapolated.
///
/// The polygon runs along the lower edge left to right and back along the upper edge, then closes.
///
/// TRUTH MATCHING. A PID entry counts as the proton when its polar angle is within dThetaMax of the
/// true emission angle. Angle alone is enough to separate it from the beam and the recoil ion,
/// which sit near zero degrees, and matching on rigidity too would bias the gate towards
/// well-reconstructed tracks -- the badly reconstructed ones are precisely what the edges must
/// contain. AtSpyralPID reports 180 - theta_lab; that convention is applied here.
///
///   root -b -q 'pid_gate_from_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003")'

void pid_gate_from_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003",
                       TString outName = "proton_sim.json", Double_t keep = 0.98, Int_t nSlice = 40,
                       Double_t xMin = 2.0, Double_t xMax = 30.0, Int_t nMin = 25,
                       Double_t dThetaMax = 10.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const double u = 931.49401, mp = 1.007825 * u;

   // one vector of rigidities per sqrt(dE/dx) slice, plus the same for every entry regardless of
   // truth, so the plot can show what the gate is being carved out of
   std::vector<std::vector<double>> slice(nSlice);
   auto *hAll = new TH2D("hAll", "simulated PID;#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
   auto *hPro = new TH2D("hPro", "truth-matched protons;#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
   long nTruth = 0, nMatched = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      // _pid.root is AtPIDTask run on its own (pidPass_a1975.C). Brho and dE/dx are computed from
      // the pattern track's hit cloud -- a circle fit for the radius, the charge integrated along
      // the first arc -- so no track fit enters the gate and none needs to be run to build it.
      TString fsim = dir + "/" + tg + "_sim.root", fpid = dir + "/" + tg + "_pid.root";
      if (gSystem->AccessPathName(fsim) || gSystem->AccessPathName(fpid)) {
         printf("  skip %s (missing)\n", tg.Data());
         continue;
      }
      TFile *fs = TFile::Open(fsim), *fp = TFile::Open(fpid);
      TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
      TTree *tp = fp ? (TTree *)fp->Get("cbmsim") : nullptr;
      // Pairing is by entry index, so the PID file may be SHORTER than the sim file but never
      // longer: reconstruction writes one entry per input event in order, so a run that was cut
      // short still aligns over its prefix. A longer PID file means the two are not the same
      // sample and no prefix can be trusted.
      if (!ts || !tp || tp->GetEntries() > ts->GetEntries()) {
         printf("\033[1;31m  %s: entry mismatch %lld sim vs %lld pid -- pairing is by index, skipping\033[0m\n",
                tg.Data(), ts ? ts->GetEntries() : -1, tp ? tp->GetEntries() : -1);
         if (fs) fs->Close();
         if (fp) fp->Close();
         continue;
      }
      if (tp->GetEntries() < ts->GetEntries())
         printf("\033[1;33m  %s: PID file is truncated, %lld of %lld events -- using the prefix\033[0m\n", tg.Data(),
                tp->GetEntries(), ts->GetEntries());
      TClonesArray *mc = nullptr, *pe = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tp->SetBranchAddress("AtPIDEvent", &pe);

      long nT = 0, nM = 0;
      for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
         ts->GetEntry(i);
         tp->GetEntry(i);
         double thT = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0)
               continue;
            thT = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (thT < 0)
            continue;
         ++nT;
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         bool got = false;
         double bd = 1e9, bx = 0, by = 0;
         for (auto &sp : ev->GetSpyral()) {
            hAll->Fill(sp.sqrtdEdx, sp.brho);
            if (!sp.valid)
               continue;
            double d = std::fabs(180.0 - sp.polar * TMath::RadToDeg() - thT);
            if (d < bd) {
               bd = d;
               bx = sp.sqrtdEdx;
               by = sp.brho;
               got = true;
            }
         }
         if (!got || bd > dThetaMax)
            continue;
         ++nM;
         hPro->Fill(bx, by);
         int s = int((bx - xMin) / (xMax - xMin) * nSlice);
         if (s >= 0 && s < nSlice)
            slice[s].push_back(by);
      }
      printf("  %-8s %6ld truth protons, %6ld truth-matched PID entries (%.1f %%)\n", tg.Data(), nT, nM,
             nT ? 100.0 * nM / nT : 0.0);
      nTruth += nT;
      nMatched += nM;
      fs->Close();
      fp->Close();
   }
   if (nMatched < 500) {
      printf("\033[1;31m  only %ld matched protons -- not enough to define edges, aborting\033[0m\n", nMatched);
      return;
   }
   printf("\n  total %ld truth protons, %ld matched (%.1f %%)\n\n", nTruth, nMatched, 100.0 * nMatched / nTruth);

   // the band, slice by slice
   std::vector<double> vx, vlo, vhi;
   printf("  sqrt(dE/dx)     n     Brho lo    median    Brho hi\n");
   for (int s = 0; s < nSlice; ++s) {
      if ((int)slice[s].size() < nMin)
         continue;
      auto &v = slice[s];
      std::sort(v.begin(), v.end());
      auto q = [&](double f) {
         double p = f * (v.size() - 1);
         int i = int(p);
         double w = p - i;
         return i + 1 < (int)v.size() ? v[i] * (1 - w) + v[i + 1] * w : v.back();
      };
      double xc = xMin + (s + 0.5) * (xMax - xMin) / nSlice;
      double lo = q((1 - keep) / 2), hi = q((1 + keep) / 2), md = q(0.5);
      vx.push_back(xc);
      vlo.push_back(lo);
      vhi.push_back(hi);
      printf("   %7.2f   %6zu    %7.3f   %7.3f   %7.3f\n", xc, v.size(), lo, md, hi);
   }
   if (vx.size() < 3) {
      printf("\033[1;31m  fewer than 3 populated slices, aborting\033[0m\n");
      return;
   }

   // polygon: along the bottom, back along the top
   std::vector<double> px, py;
   for (size_t i = 0; i < vx.size(); ++i) {
      px.push_back(vx[i]);
      py.push_back(vlo[i]);
   }
   for (int i = vx.size() - 1; i >= 0; --i) {
      px.push_back(vx[i]);
      py.push_back(vhi[i]);
   }

   // how much of the simulation each gate keeps -- the number this whole exercise is about
   auto inside = [&](double x, double y, const std::vector<double> &ax, const std::vector<double> &ay) {
      bool in = false;
      size_t n = ax.size();
      for (size_t i = 0, j = n - 1; i < n; j = i++)
         if (((ay[i] > y) != (ay[j] > y)) && (x < (ax[j] - ax[i]) * (y - ay[i]) / (ay[j] - ay[i]) + ax[i]))
            in = !in;
      return in;
   };
   auto readGate = [&](TString fn, std::vector<double> &ax, std::vector<double> &ay) {
      if (gSystem->AccessPathName(fn))
         return;
      std::ifstream in(fn.Data());
      std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      size_t p = all.find("vertices");
      while (p != std::string::npos) {
         size_t a = all.find('[', p + 1);
         if (a == std::string::npos)
            break;
         size_t b = all.find(']', a);
         if (b == std::string::npos)
            break;
         double x, y;
         char c;
         std::istringstream is(all.substr(a + 1, b - a - 1));
         if (is >> x >> c >> y) {
            ax.push_back(x);
            ay.push_back(y);
         }
         p = b;
      }
   };
   std::vector<double> ox, oy, dx, dy;
   readGate(here + "/proton_band.json", ox, oy);
   // the deuteron gate is drawn alongside for one reason: a proton gate widened on simulated
   // protons knows nothing about the neighbouring species, and if its new upper edge runs into the
   // deuteron band on the data then the extra yield it admits is contamination, not recovery
   readGate(here + "/deuteron_tight.json", dx, dy);
   long keptNew = 0, keptOld = 0;
   for (int bx = 1; bx <= hPro->GetNbinsX(); ++bx)
      for (int by = 1; by <= hPro->GetNbinsY(); ++by) {
         double c = hPro->GetBinContent(bx, by);
         if (c <= 0)
            continue;
         double x = hPro->GetXaxis()->GetBinCenter(bx), y = hPro->GetYaxis()->GetBinCenter(by);
         if (inside(x, y, px, py))
            keptNew += (long)c;
         if (ox.size() > 2 && inside(x, y, ox, oy))
            keptOld += (long)c;
      }
   printf("\n  of %ld matched protons:  new gate keeps %ld (%.1f %%)", nMatched, keptNew,
          100.0 * keptNew / nMatched);
   if (ox.size() > 2)
      printf(",  old gate keeps %ld (%.1f %%)", keptOld, 100.0 * keptOld / nMatched);
   printf("\n");

   // write it
   FILE *f = fopen(here + "/" + outName, "w");
   if (!f) {
      printf("\033[1;31m  cannot write %s\033[0m\n", (here + "/" + outName).Data());
      return;
   }
   fprintf(f, "{\n    \"name\": \"proton_sim\",\n    \"xaxis\": \"sqrtdEdx\",\n    \"yaxis\": \"brho\",\n");
   fprintf(f, "    \"vertices\": [\n");
   for (size_t i = 0; i < px.size(); ++i)
      fprintf(f, "        [\n            %.6f,\n            %.6f\n        ]%s\n", px[i], py[i],
              i + 1 < px.size() ? "," : "");
   fprintf(f, "    ]\n}\n");
   fclose(f);
   printf("  wrote %s (%zu vertices)\n", (here + "/" + outName).Data(), px.size());

   // and show it, on the simulation and on the data
   auto poly = [](const std::vector<double> &ax, const std::vector<double> &ay, int col) {
      auto *pl = new TPolyLine(ax.size() + 1);
      for (size_t i = 0; i < ax.size(); ++i)
         pl->SetPoint(i, ax[i], ay[i]);
      pl->SetPoint(ax.size(), ax[0], ay[0]);
      pl->SetLineColor(col);
      pl->SetLineWidth(3);
      return pl;
   };
   TCanvas *c1 = new TCanvas("cg", "gate from simulation", 1500, 500);
   c1->Divide(3, 1);
   c1->cd(1);
   gPad->SetLogz();
   hAll->Draw("colz");
   c1->cd(2);
   gPad->SetLogz();
   hPro->Draw("colz");
   poly(px, py, kRed + 1)->Draw("L");
   if (ox.size() > 2)
      poly(ox, oy, kGray + 2)->Draw("L");

   // the data plane, from the cached PID landscape if it is there
   c1->cd(3);
   gPad->SetLogz();
   auto *hD = new TH2D("hD", "data;#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
   int nRun = 0;
   for (int r = 106; r <= 130; ++r) {
      TString fn = Form("/mnt/f/a1975/reco/run_%04d_genfitter_pp.root", r);
      if (gSystem->AccessPathName(fn))
         continue;
      TFile *ff = TFile::Open(fn);
      TTree *t = ff ? (TTree *)ff->Get("cbmsim") : nullptr;
      if (!t) {
         if (ff) ff->Close();
         continue;
      }
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPIDEvent", &pe);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &sp : ev->GetSpyral())
            if (sp.valid)
               hD->Fill(sp.sqrtdEdx, sp.brho);
      }
      ff->Close();
      ++nRun;
   }
   hD->SetTitle(Form("data (%d runs);#sqrt{dE/dx} [arb];B#rho [T#upointm]", nRun));
   hD->Draw("colz");
   poly(px, py, kRed + 1)->Draw("L");
   if (ox.size() > 2)
      poly(ox, oy, kGray + 2)->Draw("L");
   if (dx.size() > 2)
      poly(dx, dy, kAzure + 2)->Draw("L");

   // how much extra data the new gate admits, and whether any of it is inside the deuteron gate
   long dOld = 0, dNew = 0, dBoth = 0;
   for (int bx = 1; bx <= hD->GetNbinsX(); ++bx)
      for (int by = 1; by <= hD->GetNbinsY(); ++by) {
         double c = hD->GetBinContent(bx, by);
         if (c <= 0)
            continue;
         double x = hD->GetXaxis()->GetBinCenter(bx), y = hD->GetYaxis()->GetBinCenter(by);
         bool iN = inside(x, y, px, py), iO = ox.size() > 2 && inside(x, y, ox, oy);
         if (iN)
            dNew += (long)c;
         if (iO)
            dOld += (long)c;
         if (iN && dx.size() > 2 && inside(x, y, dx, dy))
            dBoth += (long)c;
      }
   printf("\n  on the data:  old gate %ld tracks,  new gate %ld  (%+.1f %%)", dOld, dNew,
          dOld ? 100.0 * (dNew - dOld) / dOld : 0.0);
   if (dx.size() > 2)
      printf(",  of which %ld (%.2f %%) also sit inside the deuteron gate", dBoth,
             dNew ? 100.0 * dBoth / dNew : 0.0);
   printf("\n");

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/pid_gate_from_sim.png";
   c1->SaveAs(png);
   printf("  wrote %s   (red = new, grey = old)\n\n", png.Data());
}
