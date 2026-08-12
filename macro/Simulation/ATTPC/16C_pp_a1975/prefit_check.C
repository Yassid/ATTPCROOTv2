/// @file prefit_check.C
/// @brief Reconstructed kinematics against MC truth BEFORE any fitting, plus the PID.
///
/// Everything here comes from AtSpyralPID, which works off the pattern-recognition track: the
/// magnetic rigidity from the first-arc circle radius, the polar angle from a linear regression of
/// rho against z, and dE/dx from the charge over the first-arc inner segment. No genfit anywhere.
/// So this separates what the DETECTOR AND PATTERN RECOGNITION deliver from what the track fit
/// then does with it -- which matters, because the two failure modes look identical in the fitted
/// output and have completely different fixes.
///
/// Truth rigidity for a proton is p/(299.792458) with p in MeV/c, and truth polar is the emission
/// angle of the primary proton in the lab.
///
/// Panels:
///   1  Brho: reconstructed vs true          2  polar: reconstructed vs true
///   3  the residuals of both, against true Brho
///   4  the PID plane (sqrt(dE/dx) vs Brho) with the data's proton gate drawn on it
///
///   root -b -q 'prefit_check.C("/mnt/f/..._sim.root","/mnt/f/..._genfitter.root")'

void prefit_check(TString simFile, TString pidFile, TString gate = "", Double_t bField = 2.85, TString tag = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *fs = TFile::Open(simFile);
   TFile *fp = TFile::Open(pidFile);
   if (!fs || fs->IsZombie() || !fp || fp->IsZombie()) {
      printf("\033[1;31mcannot open %s or %s\033[0m\n", simFile.Data(), pidFile.Data());
      return;
   }
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TTree *tp = (TTree *)fp->Get("cbmsim");
   if (ts->GetEntries() != tp->GetEntries()) {
      printf("\033[1;31mENTRY MISMATCH %lld vs %lld -- pairing is by index\033[0m\n", ts->GetEntries(),
             tp->GetEntries());
      return;
   }
   TClonesArray *mc = nullptr, *pe = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tp->SetBranchAddress("AtPIDEvent", &pe);

   const double u = 931.49401, mp = 1.007825 * u;

   auto *hB = new TH2D("hB", "B#rho: reconstructed vs true;B#rho_{true} [T#upointm];B#rho_{PID} [T#upointm]", 100, 0,
                       1.0, 100, 0, 1.0);
   auto *hT = new TH2D("hT", "polar angle: reconstructed vs true;#theta_{true} [deg];#theta_{PID} [deg]", 100, 0, 100,
                       100, 0, 100);
   auto *rB = new TH2D("rB", "B#rho residual;B#rho_{true} [T#upointm];(B#rho_{PID} - B#rho_{true}) / B#rho_{true}",
                       50, 0, 1.0, 100, -1, 1);
   auto *rT = new TH2D("rT", "polar residual (convention applied);B#rho_{true} [T#upointm];#theta - #theta_{true} [deg]", 50, 0,
                       1.0, 100, -10, 10);
   auto *hP = new TH2D("hP", "PID (pre-fit);#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);

   long nTruth = 0, nAny = 0, nValid = 0;
   for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
      ts->GetEntry(i);
      tp->GetEntry(i);
      double keT = -1, thT = -1, brT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 2212)
            continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pp = std::sqrt(px * px + py * py + pz * pz);
         if (pp <= 0)
            continue;
         keT = std::sqrt(pp * pp + mp * mp) - mp;
         thT = std::acos(pz / pp) * TMath::RadToDeg();
         brT = pp / 299.792458; // T*m
         break;
      }
      if (keT < 0)
         continue;
      ++nTruth;
      if (!pe || pe->GetEntriesFast() == 0)
         continue;
      auto *ev = (AtPIDEvent *)pe->At(0);
      if (!ev)
         continue;
      // the PID entry whose polar angle is closest to the true one
      // Keep the best VALUES rather than a pointer: AtSpyralResult is only forward-declared in
      // the dictionary, so naming the type in a compiled macro fails to instantiate.
      bool found = false;
      double bd = 1e9, bBrho = 0, bPolar = 0;
      for (auto &sp : ev->GetSpyral()) {
         hP->Fill(sp.sqrtdEdx, sp.brho);
         if (!sp.valid)
            continue;
         double d = std::fabs(sp.polar * TMath::RadToDeg() - thT);
         if (d < bd) {
            bd = d;
            bBrho = sp.brho;
            bPolar = sp.polar;
            found = true;
         }
      }
      if (ev->GetSpyral().size())
         ++nAny;
      if (!found)
         continue;
      ++nValid;
      // CONVENTION: AtSpyralPID reports the polar angle in the opposite sense, so what it calls
      // polar is 180 - theta_lab. Verified against truth: 180 - polar - theta_true is -0.7 to
      // -2.1 deg across the whole range. Brho is unaffected, going as 1/|sin(polar)|.
      double thR = 180.0 - bPolar * TMath::RadToDeg();
      hB->Fill(brT, bBrho);
      hT->Fill(thT, thR);
      if (brT > 0)
         rB->Fill(brT, (bBrho - brT) / brT);
      rT->Fill(brT, thR - thT);
   }
   printf("\n  %ld truth protons, %ld events with any PID entry, %ld with a valid one\n", nTruth, nAny, nValid);

   // residuals in slices of true rigidity -- where does the PRE-FIT estimate go wrong?
   printf("\n  Brho_true   n     median dBrho/Brho   median dtheta [deg]\n");
   for (int b = 1; b <= rB->GetNbinsX(); b += 5) {
      TH1D *pb = rB->ProjectionY("pb", b, b + 4);
      TH1D *pt = rT->ProjectionY("pt", b, b + 4);
      if (pb->Integral() < 30) {
         delete pb;
         delete pt;
         continue;
      }
      double q = 0.5, mb, mt;
      pb->GetQuantiles(1, &mb, &q);
      pt->GetQuantiles(1, &mt, &q);
      printf("   %.3f     %6.0f   %+16.3f   %+18.2f\n", rB->GetXaxis()->GetBinCenter(b + 2), pb->Integral(), mb, mt);
      delete pb;
      delete pt;
   }

   TCanvas *c = new TCanvas("cpf", "pre-fit", 1500, 950);
   c->Divide(2, 2);
   c->cd(1);
   gPad->SetLogz();
   hB->Draw("colz");
   auto *d1 = new TLine(0, 0, 1, 1);
   d1->SetLineColor(kRed + 1);
   d1->SetLineWidth(2);
   d1->SetLineStyle(2);
   d1->Draw();
   c->cd(2);
   gPad->SetLogz();
   hT->Draw("colz");
   auto *d2 = new TLine(0, 0, 100, 100);
   d2->SetLineColor(kRed + 1);
   d2->SetLineWidth(2);
   d2->SetLineStyle(2);
   d2->Draw();
   c->cd(3);
   gPad->SetLogz();
   rB->Draw("colz");
   auto *z0 = new TLine(0, 0, 1, 0);
   z0->SetLineColor(kRed + 1);
   z0->SetLineWidth(2);
   z0->Draw();
   c->cd(4);
   gPad->SetLogz();
   hP->Draw("colz");
   // Read the gate straight from its JSON. AtTools::AtParticleID is only forward-declared in the
   // dictionary, so referring to it from a compiled macro fails to instantiate; the file is a flat
   // list of [x, y] pairs and parsing it here avoids the dependency entirely.
   if (gate.Length() && !gSystem->AccessPathName(gate)) {
      std::ifstream in(gate.Data());
      std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      std::vector<double> vx, vy;
      size_t q = all.find("vertices");
      while (q != std::string::npos) {
         size_t a = all.find('[', q + 1);
         if (a == std::string::npos)
            break;
         size_t b = all.find(']', a);
         if (b == std::string::npos)
            break;
         std::string in2 = all.substr(a + 1, b - a - 1);
         double x, y;
         char c;
         std::istringstream is(in2);
         if (is >> x >> c >> y) {
            vx.push_back(x);
            vy.push_back(y);
         }
         q = b;
      }
      if (vx.size() > 2) {
         auto *pl = new TPolyLine(vx.size() + 1);
         for (size_t k = 0; k < vx.size(); ++k)
            pl->SetPoint(k, vx[k], vy[k]);
         pl->SetPoint(vx.size(), vx[0], vy[0]);
         pl->SetLineColor(kRed + 1);
         pl->SetLineWidth(3);
         pl->Draw("L");
         printf("  gate drawn: %zu vertices\n", vx.size());
      }
   }

   gSystem->mkdir(here + "/plots", kTRUE);
   TString png = here + "/plots/prefit_check" + tag + ".png";
   c->SaveAs(png);
   printf("\n  wrote %s\n\n", png.Data());
}
