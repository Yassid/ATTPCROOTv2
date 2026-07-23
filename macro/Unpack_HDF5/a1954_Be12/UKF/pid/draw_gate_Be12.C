/// @file draw_gate_Be12.C
/// @brief INTERACTIVE PID gate drawing for a1954 12Be. Shows the (IC-gated) PID plane and
///        lets you draw a polygon with the mouse; saves it as an AtParticleID JSON gate that
///        the analysis macros read (dump_pid_Be12.C, ic_gated_pid_Be12.C, etc).
///
/// USAGE (run interactively, NOT with -b):
///   root -l 'draw_gate_Be12.C("proton_Be12")'                 // proton gate (Brho axis to 2.5)
///   root -l 'draw_gate_Be12.C("deuteron_Be12",625,750,0,1,2)' // deuteron gate (Z=1,A=2)
///
/// HOW TO DRAW: when the canvas appears, LEFT-CLICK to place each polygon vertex around the
/// band; DOUBLE-CLICK (or click the first point) to CLOSE the polygon. The gate is then saved
/// automatically to <name>.json in this folder. Re-run to redo.
///
/// Args: name      = gate name -> <name>.json  (Z,A from the name: proton=1,1 deuteron=1,2 ...)
///       icLo,icHi = IC beam gate (12Be = 625-750; set icLo<0 to use ALL beam)
///       arclenMin = min track arclength cut (0 = none; 200 recommended for a cleaner band)
///       Z,A       = charge/mass number stored in the JSON (for downstream particle mass)
///       sdMax,brMax = axis ranges (sqrtdEdx, Brho)
///       dataFile  = the persisted PID csv (default max-based; use *_integral.csv for integral scale)
void draw_gate_Be12(TString name = "proton_Be12", double icLo = 625, double icHi = 750, double arclenMin = 0,
                    int Z = 1, int A = 1, double sdMax = 55, double brMax = 2.5,
                    TString dataFile = "pid_Be12_data.csv")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString dir = gSystem->DirName(__FILE__);
   TString csv = dir + "/" + dataFile;
   if (gSystem->AccessPathName(csv)) {
      printf("\033[1;31mERROR: %s not found.\033[0m\n", csv.Data());
      return;
   }

   TTree *t = new TTree("pid", "pid");
   t->ReadFile(csv, "", ','); // first line = branch names (sqrtdedx,brho,dedx,polar,arclen,npts,ic)

   TH2F *h = new TH2F("h", TString::Format("a1954 12Be PID  -  draw the %s gate;#sqrt{dEdx};B#rho [T m]", name.Data()),
                      350, 0, sdMax, 350, 0, brMax);
   TString sel = "1";
   if (icLo >= 0)
      sel = TString::Format("ic>=%g && ic<=%g", icLo, icHi);
   if (arclenMin > 0)
      sel += TString::Format(" && arclen>%g", arclenMin);
   t->Draw("brho:sqrtdedx>>h", sel, "goff");

   TCanvas *c = new TCanvas("cGate", "draw PID gate", 950, 750);
   c->SetLogz();
   c->SetRightMargin(0.13);
   h->Draw("colz");
   c->Update();

   printf("\n\033[1;33m========================================================================\033[0m\n");
   printf("\033[1;32m Draw the %s gate:\033[0m\n", name.Data());
   printf("   * LEFT-CLICK to place each vertex around the band\n");
   printf("   * DOUBLE-CLICK (or click back on the first point) to CLOSE the polygon\n");
   printf("   IC gate: %s   arclen>%g   (%lld tracks shown)\n",
          (icLo >= 0 ? Form("%g-%g", icLo, icHi) : "ALL beam"), arclenMin, (Long64_t)h->GetEntries());
   printf("\033[1;33m========================================================================\033[0m\n\n");

   TCutG *cutg = (TCutG *)gPad->WaitPrimitive("CUTG", "CutG");
   if (!cutg || cutg->GetN() < 3) {
      printf("\033[1;31mNo valid polygon drawn (need >=3 points). Nothing saved.\033[0m\n");
      return;
   }

   // write AtParticleID JSON: {name, xaxis, yaxis, Z, A, vertices:[[x,y],...]}
   TString out = dir + "/" + name + ".json";
   std::ofstream j(out.Data());
   j << "{\n  \"name\": \"" << name << "\",\n  \"xaxis\": \"sqrtdEdx\",\n  \"yaxis\": \"brho\",\n";
   j << "  \"Z\": " << Z << ",\n  \"A\": " << A << ",\n  \"vertices\": [\n";
   int n = cutg->GetN();
   for (int i = 0; i < n; ++i) {
      double x, y;
      cutg->GetPoint(i, x, y);
      j << "    [" << x << ", " << y << "]" << (i < n - 1 ? "," : "") << "\n";
   }
   j << "  ]\n}\n";
   j.close();

   // draw the saved gate on top and report
   cutg->SetLineColor(kRed);
   cutg->SetLineWidth(3);
   cutg->Draw("L");
   c->Update();
   long inside = 0, total = 0;
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
   }
   // count via TTree::Draw with the cut
   TString incut = sel + TString::Format(" && (sqrtdedx>0)");
   double nInGate = t->Draw("brho:sqrtdedx", sel + " && CUTG", "goff");
   double nShown = h->GetEntries();
   printf("\n\033[1;32mSaved gate -> %s  (%d vertices)\033[0m\n", out.Data(), n);
   printf("In-gate: %.0f / %.0f = %.1f%% of the shown sample\n", nInGate, nShown,
          nShown > 0 ? 100.0 * nInGate / nShown : 0);
   printf("The analysis macros can now use it via this JSON.\n");
}
