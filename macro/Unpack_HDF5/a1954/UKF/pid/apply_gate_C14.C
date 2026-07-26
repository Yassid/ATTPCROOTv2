/// @file apply_gate_C14.C
/// @brief Apply/verify a PID gate JSON (drawn with draw_gate_C14.C) on the a1954 14C PID.
///        Overlays the gate, counts in-gate tracks, saves a PNG. Runs headless (root -b).
///
///   root -b -q 'apply_gate_C14.C("proton_C14.json")'
///   root -b -q 'apply_gate_C14.C("proton_C14.json",625,750,200)'   // + arclen cut
///
/// Parses the AtParticleID JSON {vertices:[[x,y],...]} into a TCutG and applies it to the
/// persisted PID data (pid_C14_data.csv). This is how the drawn gate gets "used".
void apply_gate_C14(TString gateJson = "proton_C14.json", double icLo = 625, double icHi = 750,
                     double arclenMin = 0, TString dataFile = "pid_C14_data.csv")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString dir = gSystem->DirName(__FILE__);
   TString jpath = (gateJson.BeginsWith("/") ? gateJson : dir + "/" + gateJson);
   TString csv = dir + "/" + dataFile;
   if (gSystem->AccessPathName(jpath)) {
      printf("\033[1;31mgate not found: %s\033[0m\n", jpath.Data());
      return;
   }

   // --- parse vertices from the JSON (simple, no JSON lib) ---
   std::ifstream in(jpath.Data());
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto vpos = s.find("vertices");
   auto lb = s.find('[', vpos);
   auto rb = s.rfind(']');
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
   TCutG *cut = new TCutG("CUTG", npt + 1);
   for (int i = 0; i < npt; ++i)
      cut->SetPoint(i, nums[2 * i], nums[2 * i + 1]);
   cut->SetPoint(npt, nums[0], nums[1]); // close
   cut->SetVarX("sqrtdedx");
   cut->SetVarY("brho");
   cut->SetLineColor(kRed);
   cut->SetLineWidth(3);

   TTree *t = new TTree("pid", "pid");
   t->ReadFile(csv, "", ',');
   TString sel = (icLo >= 0) ? TString::Format("ic>=%g && ic<=%g", icLo, icHi) : TString("1");
   if (arclenMin > 0)
      sel += TString::Format(" && arclen>%g", arclenMin);

   double sdMax = 55, brMax = 2.5;
   for (int i = 0; i < npt; ++i) {
      sdMax = std::max(sdMax, nums[2 * i] * 1.15);
      brMax = std::max(brMax, nums[2 * i + 1] * 1.15);
   }
   TH2F *h = new TH2F("h", TString::Format("14C PID + %s;#sqrt{dEdx};B#rho [T m]", gateJson.Data()), 350, 0, sdMax,
                      350, 0, brMax);
   t->Draw("brho:sqrtdedx>>h", sel, "goff");
   double nShown = h->GetEntries();
   double nIn = t->Draw("brho:sqrtdedx", sel + " && CUTG", "goff");

   TCanvas *c = new TCanvas("c", "gate", 950, 750);
   c->SetLogz();
   c->SetRightMargin(0.13);
   h->Draw("colz");
   cut->Draw("L");
   TString png = dir + "/plots/applied_" + gateJson;
   png.ReplaceAll(".json", ".png");
   c->SaveAs(png);

   printf("\n\033[1;32m=== %s ===\033[0m\n", gateJson.Data());
   printf("selection: %s\n", sel.Data());
   printf("in-gate: %.0f / %.0f shown = %.1f%%\n", nIn, nShown, nShown > 0 ? 100.0 * nIn / nShown : 0);
   printf("saved %s\n", png.Data());
}
