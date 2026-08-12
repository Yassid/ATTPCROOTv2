/// @file make_sim_gate.C
/// @brief Build the simulated-proton PID gate automatically, no hand drawing.
///
/// Every point in the cache is a TRUTH-MATCHED simulated proton, so there is nothing to separate
/// here -- no contamination to exclude, no second band to avoid. The gate just has to enclose the
/// locus. That makes hand drawing unnecessary: slice in sqrt(dE/dx), take a lower and upper
/// quantile of Brho in each slice, and join the two envelopes into a closed polygon. The result is
/// reproducible and its efficiency is set by construction rather than by where the mouse went.
///
/// keepFrac sets the quantiles: 0.98 takes the 1st and 99th percentile of Brho per slice, so about
/// 98 % of protons in every slice fall inside, and the efficiency is flat in sqrt(dE/dx) instead of
/// being high in the dense part of the locus and poor in the tails. A gate drawn by eye tends to do
/// the opposite, because the sparse ends are where the eye has least to go on -- and on this data
/// set the sparse end at low Brho is exactly the forward-theta_cm population that matters.
///
/// Slices with fewer than minPerSlice entries are dropped: a quantile over 3 points is noise, and
/// letting it set a vertex puts a spike in the boundary.
///
///   root -b -q 'make_sim_gate.C(0.98,"pid_gate_sim_mp15.json","plots/sim_proton_points_mp15.root")'

void make_sim_gate(Double_t keepFrac = 0.98, TString outJson = "pid_gate_sim_mp15.json",
                   TString cache = "plots/sim_proton_points_mp15.root", Int_t nSlice = 40,
                   Int_t minPerSlice = 40, Double_t pad = 0.02, TString png = "plots/sim_gate_mp15.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   TFile *f = TFile::Open(here + "/" + cache);
   TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("pts") : nullptr;
   if (!t) { printf("\033[1;31mno tree 'pts' in %s\033[0m\n", cache.Data()); return; }
   float x, y;
   t->SetBranchAddress("x", &x);
   t->SetBranchAddress("y", &y);
   std::vector<double> X, Y;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); X.push_back(x); Y.push_back(y); }
   f->Close();
   if (X.size() < 100) { printf("\033[1;31monly %zu points\033[0m\n", X.size()); return; }

   double xlo = *std::min_element(X.begin(), X.end()), xhi = *std::max_element(X.begin(), X.end());
   printf("\n  %zu simulated protons, sqrt(dE/dx) in [%.2f, %.2f]\n", X.size(), xlo, xhi);

   // per-slice Brho quantiles
   const double qLo = 0.5 * (1.0 - keepFrac), qHi = 1.0 - qLo;
   double w = (xhi - xlo) / nSlice;
   std::vector<double> cx, clo, chi;
   for (int s = 0; s < nSlice; ++s) {
      double a = xlo + s * w, b = a + w;
      std::vector<double> ys;
      for (size_t i = 0; i < X.size(); ++i)
         if (X[i] >= a && X[i] < b) ys.push_back(Y[i]);
      if ((int)ys.size() < minPerSlice) continue;
      std::sort(ys.begin(), ys.end());
      double lo = ys[(size_t)(qLo * (ys.size() - 1))];
      double hi = ys[(size_t)(qHi * (ys.size() - 1))];
      double m = pad * (hi - lo);            // a small margin so the boundary is not ON the points
      cx.push_back(0.5 * (a + b));
      clo.push_back(std::max(0.0, lo - m));
      chi.push_back(hi + m);
   }
   if (cx.size() < 3) { printf("\033[1;31monly %zu usable slices\033[0m\n", cx.size()); return; }

   // closed polygon: left to right along the bottom, right to left along the top
   std::vector<double> px, py;
   for (size_t i = 0; i < cx.size(); ++i) { px.push_back(cx[i]); py.push_back(clo[i]); }
   for (size_t i = cx.size(); i-- > 0;)     { px.push_back(cx[i]); py.push_back(chi[i]); }

   auto inside = [&](double qx, double qy) {
      bool in = false;
      size_t n = px.size();
      for (size_t i = 0, j = n - 1; i < n; j = i++)
         if (((py[i] > qy) != (py[j] > qy)) && (qx < (px[j] - px[i]) * (qy - py[i]) / (py[j] - py[i]) + px[i]))
            in = !in;
      return in;
   };
   long in = 0;
   for (size_t i = 0; i < X.size(); ++i) if (inside(X[i], Y[i])) ++in;
   printf("  %zu slices used, %zu vertices\n", cx.size(), px.size());
   printf("  keeps %ld / %zu simulated protons = %.2f %%\n", in, X.size(), 100.0 * in / X.size());

   // efficiency where it is easiest to lose it: the low-Brho end
   long nlo = 0, inlo = 0;
   for (size_t i = 0; i < X.size(); ++i)
      if (Y[i] < 0.25) { ++nlo; if (inside(X[i], Y[i])) ++inlo; }
   printf("  below Brho 0.25: keeps %ld / %ld = %.2f %%\n\n", inlo, nlo, nlo ? 100.0 * inlo / nlo : 0.0);

   FILE *o = fopen((here + "/" + outJson).Data(), "w");
   if (!o) { printf("\033[1;31mcannot write %s\033[0m\n", outJson.Data()); return; }
   fprintf(o, "{\n    \"name\": \"proton_sim_auto\",\n    \"xaxis\": \"sqrtdedx\",\n    \"yaxis\": \"brho\",\n"
              "    \"vertices\": [\n");
   for (size_t i = 0; i < px.size(); ++i)
      fprintf(o, "        [%.3f, %.3f]%s\n", px[i], py[i], (i + 1 == px.size()) ? "" : ",");
   fprintf(o, "    ]\n}\n");
   fclose(o);

   auto *h = new TH2F("hg", "simulated protons + auto gate;#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0,
                      1.2);
   for (size_t i = 0; i < X.size(); ++i) h->Fill(X[i], Y[i]);
   TCanvas *c = new TCanvas("cg", "sim gate", 950, 750);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   h->Draw("colz");
   auto *pl = new TPolyLine(px.size() + 1);
   for (size_t i = 0; i < px.size(); ++i) pl->SetPoint(i, px[i], py[i]);
   pl->SetPoint(px.size(), px[0], py[0]);
   pl->SetLineColor(kGreen + 2);
   pl->SetLineWidth(3);
   pl->Draw("L");
   c->SaveAs(here + "/" + png);
   printf("  wrote %s\n         %s\n\n", outJson.Data(), png.Data());
}
