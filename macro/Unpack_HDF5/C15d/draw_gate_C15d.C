/// @file draw_gate_C15d.C
/// @brief Draw a PID gate interactively on the ATTPCROOT gain-matched plane.
///
/// MUST run WITHOUT -b (it needs a canvas you can click on):
///
///   root -l 'draw_gate_C15d.C("proton_C15d", 1, 1)'
///   root -l 'draw_gate_C15d.C("deuteron_C15d", 1, 2)'
///
/// Left-click each vertex around the band, then DOUBLE-CLICK to close the polygon. The cut is
/// written as spyral_utils Cut2D JSON (plus Z and A), so it loads with AtCut2D/AtParticleID and
/// with Spyral itself.
///
/// ★ GATES FOR THIS ANALYSIS ARE DRAWN HERE, ON THIS PLANE. Do not import a gate from another
/// pipeline: at matched gain the two agree on brho to better than a percent but NOT on
/// sqrt(dE/dx), which is ~2.8x smaller here (a factor ~8 in dE/dx). Rigidity is geometry and
/// reproduces; dE/dx is charge over arclength and charge is not measured the same way. A
/// foreign gate landing on this plane selects whatever happens to lie there and still returns
/// a plausible-looking count, which is the worst kind of wrong.
///
/// spyralOverlay/spyralScaleX exist only to sketch a reference outline while drawing. Whatever
/// they draw is a dashed grey line on the canvas and nothing more -- it is never written to the
/// gate, and it is not a gate for this plane.

void draw_gate_C15d(TString gateName = "proton_C15d", Int_t Z = 1, Int_t A = 1,
                    TString planeFile = "plots/pid_C15d.root", TString outDir = "gates/",
                    TString spyralOverlay = "", Double_t spyralScaleX = 1.0)
{
   gSystem->Load("libAtTools.so");

   if (gROOT->IsBatch()) {
      std::cout << "\033[1;31mERROR: this macro is interactive -- run it WITHOUT -b.\033[0m\n";
      return;
   }
   if (gSystem->AccessPathName(planeFile.Data())) {
      std::cout << "\033[1;31mERROR: " << planeFile << " not found. Run mkpid_C15d.C first.\033[0m\n";
      return;
   }

   TFile *f = TFile::Open(planeFile);
   auto *h = dynamic_cast<TH2 *>(f->Get("hpid"));
   if (h == nullptr) {
      std::cout << "\033[1;31mERROR: no 'hpid' in " << planeFile << "\033[0m\n";
      return;
   }

   auto *c = new TCanvas("cgate", "draw the gate: left-click vertices, DOUBLE-CLICK to close", 1100, 850);
   c->SetLogz();
   c->SetRightMargin(0.13);
   h->Draw("colz");

   if (spyralOverlay.Length()) {
      // gSystem->ExpandPathName has TWO overloads: the char* one RETURNS the expanded
      // path, the TString one expands IN PLACE and returns a Bool_t error flag.

      // Assigning that Bool_t back to the TString blanks it, and the caller then
      // reports "cannot load" against an empty filename.
      gSystem->ExpandPathName(spyralOverlay);
      auto ref = AtTools::AtCut2D::LoadJSON(spyralOverlay.Data());
      if (ref.IsValid()) {
         const auto &v = ref.GetVertices();
         auto *g = new TGraph(static_cast<Int_t>(v.size()) + 1);
         for (size_t k = 0; k < v.size(); ++k)
            g->SetPoint(static_cast<Int_t>(k), v[k].first * spyralScaleX, v[k].second);
         g->SetPoint(static_cast<Int_t>(v.size()), v[0].first * spyralScaleX, v[0].second);
         g->SetLineColor(kGray + 2);
         g->SetLineStyle(2);
         g->SetLineWidth(2);
         g->Draw("L same");
         std::cout << "  reference overlay: " << ref.GetName() << " with x scaled by " << spyralScaleX
                   << " (dashed grey) -- REFERENCE ONLY, not a gate for this plane\n";
      }
   }

   std::cout << "\033[1;33m=== draw the '" << gateName << "' gate ===\033[0m\n"
             << "  left-click each vertex around the band, DOUBLE-CLICK to close.\n";

   c->Update();
   auto *cut = dynamic_cast<TCutG *>(c->WaitPrimitive("CUTG", "CutG"));
   if (cut == nullptr) {
      std::cout << "\033[1;31mNo cut drawn.\033[0m\n";
      return;
   }

   std::vector<std::pair<double, double>> verts;
   for (Int_t i = 0; i < cut->GetN(); ++i) {
      Double_t x, y;
      cut->GetPoint(i, x, y);
      verts.emplace_back(x, y);
   }
   // TCutG repeats the first point to close itself; Cut2D does not want the duplicate.
   if (verts.size() > 2 && std::abs(verts.front().first - verts.back().first) < 1e-9 &&
       std::abs(verts.front().second - verts.back().second) < 1e-9)
      verts.pop_back();

   if (verts.size() < 3) {
      std::cout << "\033[1;31mNeed at least 3 vertices; got " << verts.size() << ".\033[0m\n";
      return;
   }

   gSystem->mkdir(outDir, kTRUE);
   AtTools::AtCut2D c2d(gateName.Data(), verts, "sqrt_dEdx", "brho");
   AtTools::AtParticleID pid(c2d, Z, A);
   TString outPath = outDir + gateName + ".json";
   if (!pid.WriteJSON(outPath.Data())) {
      std::cout << "\033[1;31mERROR: could not write " << outPath << "\033[0m\n";
      return;
   }

   // Report what it holds straight away: a gate whose count you only learn later is a gate
   // you will re-draw later.
   Long64_t in = 0, tot = 0;
   for (Int_t ix = 1; ix <= h->GetNbinsX(); ++ix)
      for (Int_t iy = 1; iy <= h->GetNbinsY(); ++iy) {
         const double n = h->GetBinContent(ix, iy);
         if (n <= 0)
            continue;
         tot += static_cast<Long64_t>(n);
         if (c2d.IsInside(h->GetXaxis()->GetBinCenter(ix), h->GetYaxis()->GetBinCenter(iy)))
            in += static_cast<Long64_t>(n);
      }

   std::cout << "\033[1;32mwrote\033[0m " << outPath << "  (" << verts.size() << " vertices, Z=" << Z
             << " A=" << A << ")\n"
             << "  holds " << in << " of " << tot << " plane entries = " << (tot ? 100.0 * in / tot : 0.)
             << "%\n"
             << "  apply it with:\n"
             << "    root -b -q 'apply_gate_C15d.C(\"" << outPath << "\")'\n";
}
