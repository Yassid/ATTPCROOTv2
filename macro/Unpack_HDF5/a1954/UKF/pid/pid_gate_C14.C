/// @file pid_gate_C14.C
/// @brief a1954 14C PID gate: the PID plane ALONE, with the drawn proton_14C gate overlaid.
///        This is the left panel of the two-panel figure the write-up ships as
///        figures/pid_gate_C14.png (identical to pid/plots/gate_drawn_check.png, produced
///        one-off on 2026-07-26 by no macro in the tree -- hence this one). The right panel
///        of that figure, polar angle all vs in-gate, is NOT reproduced here.
///
///        Plane = the persisted PID data pid_C14_data.csv written by pid_C14.C
///        (columns sqrtdedx, brho, dedx, polar, arclen, npts, ic), gate = the hand-drawn
///        proton_14C.json used in production. Nothing is recomputed, so the plane is exactly
///        the one the gate was drawn on.
///
///   root -b -q 'pid_gate_C14.C'                            // as shipped: no IC gate
///   root -b -q 'pid_gate_C14.C("proton_14C.json",950,1350)' // + the production IC window
///   root -b -q 'pid_gate_C14.C("proton_14C_auto.json",-1,-1,0,"_auto")'
void pid_gate_C14(TString gateJson = "proton_14C.json", double icLo = -1, double icHi = -1, double arclenMin = 0,
                  TString outTag = "", TString dataFile = "pid_C14_data.csv", double sdMax = 55, double brMax = 1.2)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   // Inputs resolve beside the macro first, then in ../data/ -- so the same file works both in
   // the ATTPCROOT tree and in the self-contained a1954_C14_absnorm repo, which vendors them.
   TString dir = gSystem->DirName(__FILE__);
   auto resolve = [&dir](TString f) {
      if (f.BeginsWith("/"))
         return f;
      TString beside = dir + "/" + f;
      return gSystem->AccessPathName(beside) ? TString(dir + "/../data/" + f) : beside;
   };
   TString jpath = resolve(gateJson), csv = resolve(dataFile);
   for (TString p : {jpath, csv})
      if (gSystem->AccessPathName(p)) {
         printf("\033[1;31mmissing: %s\033[0m\n", p.Data());
         return;
      }

   // --- the gate: parse {"vertices":[[x,y],...]} anchored on "vertices", never on the whole
   //     file -- the header carries name/Z/A and a blind number scan prepends 1 and 1.
   std::ifstream in(jpath.Data());
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto lb = s.find('[', s.find("vertices"));
   auto rb = s.rfind(']');
   if (lb == std::string::npos || rb == std::string::npos) {
      printf("\033[1;31mno vertices block in %s\033[0m\n", jpath.Data());
      return;
   }
   std::vector<double> nums;
   for (const char *p = s.c_str() + lb; p < s.c_str() + rb;) {
      char *np = nullptr;
      double v = strtod(p, &np);
      if (np != p) {
         nums.push_back(v);
         p = np;
      } else
         ++p;
   }
   int npt = nums.size() / 2;
   if (npt < 3) {
      printf("\033[1;31mgate has <3 vertices\033[0m\n");
      return;
   }
   auto *cut = new TCutG("CUTG", npt + 1);
   for (int i = 0; i < npt; ++i)
      cut->SetPoint(i, nums[2 * i], nums[2 * i + 1]);
   cut->SetPoint(npt, nums[0], nums[1]); // close it
   cut->SetVarX("sqrtdedx");
   cut->SetVarY("brho");
   cut->SetLineColor(kRed);
   cut->SetLineWidth(3);

   // --- the plane, straight off the csv
   auto *t = new TTree("pid", "pid");
   t->ReadFile(csv, "", ',');
   TString sel = (icLo >= 0) ? TString::Format("ic>=%g && ic<=%g", icLo, icHi) : TString("1");
   if (arclenMin > 0)
      sel += TString::Format(" && arclen>%g", arclenMin);

   auto *h = new TH2F("h", TString::Format("a1954 14C PID + the drawn %s gate;#sqrt{dEdx};B#rho [T m]",
                                           TString(gateJson).ReplaceAll(".json", "").Data()),
                      350, 0, sdMax, 350, 0, brMax);
   t->Draw("brho:sqrtdedx>>h", sel, "goff");
   double nSel = t->GetEntries(sel);
   double nIn = t->GetEntries(sel + " && CUTG");

   TCanvas *c = new TCanvas("cPidGate", "PID gate", 800, 700);
   c->SetLogz();
   c->SetRightMargin(0.13);
   c->SetTicks(1, 1);
   h->Draw("colz");
   cut->Draw("L");
   auto *tx = new TLatex();
   tx->SetNDC();
   tx->SetTextSize(0.032);
   tx->DrawLatex(0.42, 0.86, TString::Format("in gate: %.0f / %.0f = %.1f %%", nIn, nSel, nSel > 0 ? 100. * nIn / nSel : 0.));

   TString plotDir = dir + "/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE); // the absnorm copy has no plots/ dir; SaveAs fails mutely without this
   TString png = plotDir + "pid_gate_C14" + outTag + ".png";
   c->SaveAs(png);
   if (gSystem->AccessPathName(png)) {
      printf("\033[1;31mSaveAs did not write %s\033[0m\n", png.Data());
      return;
   }
   printf("\nselection: %s\nin gate: %.0f / %.0f = %.1f %%  (plane %.0f rows, %s)\nsaved %s\n", sel.Data(), nIn, nSel,
          nSel > 0 ? 100. * nIn / nSel : 0., (double)t->GetEntries(), gSystem->BaseName(csv), png.Data());
}
