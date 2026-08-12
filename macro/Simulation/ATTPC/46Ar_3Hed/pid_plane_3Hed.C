/// @file pid_plane_3Hed.C
/// @brief The simulated PID plane for 46Ar(3He,d)47K. No gate, no cut, no truth by default.
///
/// This is the plane the deuteron gate gets drawn on BY HAND. It shows every pattern track the
/// reconstruction found, from every sample given, with nothing removed.
///
///   root -b -q 'pid_plane_3Hed.C()'                                    // all six default samples
///   root -b -q 'pid_plane_3Hed.C("/mnt/f/ar46_3hed","gs_s3001")'       // one sample
///
/// WHAT TO EXPECT ON IT. Measured on a 400-entry probe sample (g.s., seed 7001):
///   - the deuteron band at Brho 0.47 to 1.56 T*m, sqrt(dE/dx) 2.5 to 9.6. The rigidity range is
///     the kinematics and nothing else: 5.1 to 55.3 MeV of deuteron over theta_cm 15-80 deg is
///     p = 139 to 459 MeV/c, i.e. Brho = 0.46 to 1.53 T*m. If the band sits somewhere else, the
///     field or the charge scale is wrong, not the reaction.
///   - the 46Ar beam, Z = 18 and depositing ~9x the charge per unit length of a1975's 16C beam,
///     wherever it leaks past the inhibited beam hole. On the probe sample it does not appear at
///     all, which is the beam hole doing its job -- do not assume that survives a bigger sample.
/// The 47K recoil is slow and heavily ionising but travels only millimetres, so it mostly does
/// not survive pattern recognition at all.
///
/// ONLY VALID RESULTS ARE PLOTTED, and on the probe sample that was 89 of 187 pattern tracks
/// (48 %). The rest are AtSpyralPID early exits -- short, tightly curled tracks failing the
/// spline or circle fit are the expected population at forward theta_cm, where the helix is only
/// 24 cm across. Plotting them would pile ~half the sample onto the origin and distort the colour
/// scale. The count is printed per sample so the loss is visible rather than silent.
///
/// NOTE that AtSpyralResult::failCode is TRANSIENT (//! in the header), so it does NOT survive
/// into the _pid.root file -- every entry reads back as failCode 0. Attribute failures with
/// spyral_failcodes.C over the reco file instead; `valid` is what this macro can rely on.
///
/// NO TRUTH PANEL HERE, ON PURPOSE. The hits carry MC truth (run_reco_Ar46_TC.C calls
/// SetSaveMCInfo), so a truth-labelled version of this plane can be built -- but it belongs AFTER
/// a gate exists, as the check on it. A gate drawn by tracing a truth panel measures the truth
/// rather than the detector, and its efficiency then means nothing. Ask for that macro once the
/// gate is drawn; it has to reach back to AtPatternEvent for the per-hit truth, which is a
/// different input from the PID file this macro reads.

/// TWO SOURCES, one plane. Leave pointsFile empty and it reads the chain's own _pid.root files,
/// which AtPIDTask wrote at the fixed class default fMinPoints = 30. Give it a points file from
/// pid_points_3Hed.C and it reads that instead, at whatever cut that file was built with -- which
/// is the only way to see the theta_lab 77-104 deg band at all. The title records which.
///
///   root -b -q 'pid_plane_3Hed.C("/mnt/f/ar46_3hed","","plots/pid_mp15.png","plots/pid_points_mp15.root")'

void pid_plane_3Hed(TString dir = "/mnt/f/ar46_3hed",
                    TString tags = "gs_s3001,gs_s3002,360_s3011,360_s3012,2020_s3021,2020_s3022",
                    TString png = "plots/pid_plane_3Hed.png", TString pointsFile = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   // Axes are set for THIS reaction, not inherited from a1975: deuterons in 3He+CO2 reach
   // Brho 1.56 T*m, well past the 1.2 the proton/H2 plane was drawn to.
   auto *h = new TH2D("hPID", "simulated PID, 46Ar(^{3}He,d)^{47}K;#sqrt{dE/dx} [arb];B#rho [T#upointm]",
                      300, 0, 20, 250, 0, 2.0);
   long n = 0;

   if (!pointsFile.IsNull()) {
      TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
      TString fp = pointsFile.BeginsWith("/") ? pointsFile : here + "/" + pointsFile;
      TFile *f = TFile::Open(fp);
      TTree *t = f ? (TTree *)f->Get("pts") : nullptr;
      if (!t) {
         printf("  no 'pts' tree in %s -- build it with pid_points_3Hed.C first\n", fp.Data());
         return;
      }
      float x, y;
      t->SetBranchAddress("x", &x);
      t->SetBranchAddress("y", &y);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         h->Fill(x, y);
         ++n;
      }
      // The cut is in the file name, not in the tree, so the title takes it from there rather than
      // claiming a number this macro cannot verify.
      TString cut = fp(TRegexp("mp[0-9]+"));
      h->SetTitle(TString::Format("simulated PID, 46Ar(^{3}He,d)^{47}K  [%s];#sqrt{dE/dx} [arb];B#rho [T#upointm]",
                                  cut.IsNull() ? "points file" : cut.Data()));
      printf("  %ld points from %s\n", n, fp.Data());
      f->Close();
   } else {

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      if (tg.IsNull())
         continue;
      TString fn = dir + "/" + tg + "_pid.root";
      if (gSystem->AccessPathName(fn)) {
         printf("  skip %-12s (missing)\n", tg.Data());
         continue;
      }
      TFile *f = TFile::Open(fn);
      TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t) {
         if (f)
            f->Close();
         continue;
      }
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPIDEvent", &pe);
      long nt = 0, nInvalid = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &sp : ev->GetSpyral()) {
            if (!sp.valid) {
               ++nInvalid;
               continue;
            }
            h->Fill(sp.sqrtdEdx, sp.brho);
            ++nt;
         }
      }
      printf("  %-12s %7ld valid   %6ld rejected by AtSpyralPID (%.0f%%)\n", tg.Data(), nt, nInvalid,
             (nt + nInvalid) ? 100. * nInvalid / (nt + nInvalid) : 0.);
      n += nt;
      f->Close();
   }
   delete ta;
   } // end of the _pid.root source

   printf("\n  %ld valid PID entries plotted\n", n);
   if (n == 0) {
      printf("  nothing to draw -- run the accumulation first\n");
      return;
   }

   TCanvas *c = new TCanvas("cPID", "simulated PID", 900, 700);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   h->Draw("colz");
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("  wrote %s\n", png.Data());
}
