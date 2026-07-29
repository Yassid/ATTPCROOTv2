/// @file mkpid_C15.C
/// @brief Build the full-statistics a2091 15C PID plane from the MERGED per-run Spyral
///        ntuples (pidall_C15.sh -> pid_C15_all.root), and write pid_C15_data.csv in the
///        exact column order draw_gate_C15.C / apply_gate_C15.C expect:
///            sqrtdedx,brho,dedx,polar,arclen,npts,ic
///        There is no ion-chamber value available yet (the 11 raw-merger runs have no
///        /frib group at all, and the FRIB extraction has not been run for the legacy
///        ones), so `ic` is written as -1 and the gate must be drawn with icLo<0
///        ("ALL beam"). Re-run this once <run>_FRIB.root exist to get a real IC gate.
///
///   root -b -q 'pid/mkpid_C15.C("pid/plots/pid_C15_all.root")'
/// brhoMax defaults to 3.0, not 1.6: at 1.6 the plane clipped 80,826 tracks (9.8%), of which
/// 73,317 survive arclen>200 -- so the arclen cut does NOT remove them. 98.8% of that high-Brho
/// population is polar<25 or >155 deg, i.e. the Brho = B*radius/|sin(polar)| divergence near the
/// beam axis; only 887 tracks above 1.6 T m are at mid-polar. Cut on POLAR to remove it, not arclen.
void mkpid_C15(TString inFile = "pid/plots/pid_C15_all.root", double sdMax = 60, double brhoMax = 3.0,
               double dedxMax = 3000, double arclenMin = 0)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TFile *f = TFile::Open(inFile);
   if (!f || f->IsZombie()) { printf("ERROR: cannot open %s\n", inFile.Data()); return; }
   TNtuple *nt = (TNtuple *)f->Get("spid");
   if (!nt) { printf("ERROR: no `spid` ntuple in %s\n", inFile.Data()); return; }

   float dedx, sqrtdedx, brho, polar, arclen, npts, vtxr, vtxz, radius;
   nt->SetBranchAddress("dedx", &dedx);       nt->SetBranchAddress("sqrtdedx", &sqrtdedx);
   nt->SetBranchAddress("brho", &brho);       nt->SetBranchAddress("polar", &polar);
   nt->SetBranchAddress("arclen", &arclen);   nt->SetBranchAddress("npts", &npts);
   nt->SetBranchAddress("vtxr", &vtxr);       nt->SetBranchAddress("vtxz", &vtxz);
   nt->SetBranchAddress("radius", &radius);

   TString dir = gSystem->DirName(__FILE__);
   TString plotDir = dir + "/plots/";
   TString csvPath = dir + "/pid_C15_data.csv";
   std::ofstream csv(csvPath.Data());
   csv << "sqrtdedx,brho,dedx,polar,arclen,npts,ic\n";

   // main PID plane, plus the two cuts most likely to clean the band
   TH2F *hS   = new TH2F("hS",  "15C Spyral PID, all tracks;#sqrt{dEdx};B#rho [T m]", 500, 0, sdMax, 500, 0, brhoMax);
   // DIAGNOSTIC ONLY -- kept to show what the cut does, not as a gating plane. The 200 mm cut
   // carves a false gap through the continuous proton band (see the file header of draw_gate_C15.C).
   TH2F *hSa  = new TH2F("hSa", "DIAGNOSTIC: arclen>200 carves a FALSE gap -- do not gate on this;#sqrt{dEdx};B#rho [T m]",
                         500, 0, sdMax, 500, 0, brhoMax);
   TH2F *hSd  = new TH2F("hSd", "15C Spyral PID;dEdx [counts];B#rho [T m]", 500, 0, dedxMax, 500, 0, brhoMax);
   TH2F *hKt  = new TH2F("hKt", "dEdx vs #theta_{lab};#theta_{lab} [deg];dEdx [counts]", 180, 0, 180, 400, 0, dedxMax);
   TH1F *hAl  = new TH1F("hAl", "track arclength;arclen [mm];tracks", 200, 0, 1000);
   TH1F *hNp  = new TH1F("hNp", "points per track;npts;tracks", 150, 0, 300);

   Long64_t N = nt->GetEntries(), kept = 0;
   for (Long64_t i = 0; i < N; i++) {
      nt->GetEntry(i);
      hAl->Fill(arclen); hNp->Fill(npts);
      if (arclenMin > 0 && arclen <= arclenMin) continue;
      hS->Fill(sqrtdedx, brho);
      hSd->Fill(dedx, brho);
      hKt->Fill(polar, dedx);
      if (arclen > 200) hSa->Fill(sqrtdedx, brho);
      csv << sqrtdedx << "," << brho << "," << dedx << "," << polar << ","
          << arclen << "," << (int)npts << ",-1\n";
      kept++;
   }
   csv.close();

   printf("\n===== a2091 15C PID, FULL STATISTICS =====\n");
   printf("tracks in merged ntuple : %lld\n", N);
   printf("written to CSV          : %lld  -> %s\n", kept, csvPath.Data());
   printf("arclen>200              : %.0f (%.1f%%)\n", hSa->GetEntries(), 100.0 * hSa->GetEntries() / N);
   printf("arclen  : mean %.1f mm  |  npts: mean %.1f\n", hAl->GetMean(), hNp->GetMean());
   printf("Brho    : mean %.3f  rms %.3f T m\n", hS->GetMean(2), hS->GetRMS(2));
   printf("polar   : %.1f - %.1f deg (1-99%%)\n",
          hKt->GetXaxis()->GetBinCenter(hKt->ProjectionX()->FindFirstBinAbove(0.01 * hKt->ProjectionX()->GetMaximum())),
          hKt->GetXaxis()->GetBinCenter(hKt->ProjectionX()->FindLastBinAbove(0.01 * hKt->ProjectionX()->GetMaximum())));

   TCanvas *c = new TCanvas("c", "pid", 1700, 950);
   c->Divide(3, 2);
   c->cd(1); gPad->SetLogz(); hS->Draw("colz");
   c->cd(2); gPad->SetLogz(); hSa->Draw("colz");
   c->cd(3); gPad->SetLogz(); hSd->Draw("colz");
   c->cd(4); gPad->SetLogz(); hKt->Draw("colz");
   c->cd(5); gPad->SetLogy(); hAl->Draw();
   c->cd(6); gPad->SetLogy(); hNp->Draw();
   c->SaveAs(plotDir + "pid_C15_FULL.png");

   // the plane on its own, large, for gate drawing by eye
   TCanvas *c2 = new TCanvas("c2", "pid plane", 1100, 850);
   c2->SetLogz(); hS->SetTitle("a2091 15C Spyral PID -- draw the gate here;#sqrt{dEdx};B#rho [T m]");
   hS->Draw("colz");
   c2->SaveAs(plotDir + "pid_C15_PLANE.png");
   TCanvas *c3 = new TCanvas("c3", "pid plane arclen", 1100, 850);
   c3->SetLogz(); hSa->Draw("colz");
   c3->SaveAs(plotDir + "pid_C15_PLANE_arclen200.png");

   printf("\nplots: %spid_C15_FULL.png\n       %spid_C15_PLANE.png\n       %spid_C15_PLANE_arclen200.png\n",
          plotDir.Data(), plotDir.Data(), plotDir.Data());
   printf("\nTo draw the gate interactively (needs a GUI, do NOT use -b):\n");
   printf("  root -l 'macro/Unpack_HDF5/a2091/UKF/pid/draw_gate_C15.C(\"proton_15C\")'\n");
   printf("  (defaults are now correct: ALL beam, NO arclen cut, axes 60 x 1.5)\n");
   printf("  NOTE: do NOT add an arclen cut -- it carves a false gap in the proton band,\n");
   printf("        see pid_C15_gap_diag.png. The arclen200 plane below is a DIAGNOSTIC only.\n");
}
