/// @file make_union_gate.C
/// @brief Union of the hand-drawn production gate and the simulation-built one.
///
/// WHY A UNION AND NOT A REPLACEMENT. Neither gate can be used on its own:
///   proton_band.json      -- drawn by hand on the data. Keeps 61.2 % of simulated protons; its
///                            upper edge cuts through the high-rigidity band (1.6 % kept above
///                            30 MeV) and a notch near sqrt(dE/dx) 8-11 costs 18 % of the 5-10 MeV
///                            protons, which lands on theta_cm 50-60 deg -- on top of the
///                            diffraction minimum.
///   proton_sim_*.json     -- quantiles of truth-matched simulated protons. Keeps 97.3 %, but the
///                            simulation runs out of statistics at sqrt(dE/dx) ~ 23 while the data
///                            ridge continues to 40, so on data it admits 38.8 % FEWER tracks.
/// Taking, in each sqrt(dE/dx) slice, lo = min(lo) and hi = max(hi) gives a band that is never
/// narrower than what the production runs today and wider exactly where the simulation shows real
/// protons. It cannot lose anything currently kept, which is what makes it safe to adopt.
///
/// The polygon extent at a given x is read by crossing its edges rather than from the vertex list,
/// so a gate whose vertices are unevenly spaced is still sampled correctly.
///
///   root -b -q 'pid/make_union_gate.C()'

void make_union_gate(TString gateA = "pid/proton_band.json", TString gateB = "pid/proton_sim_6seeds.json",
                     TString outName = "pid/proton_union.json", Int_t nSlice = 80, Double_t xMin = 2.0,
                     Double_t xMax = 40.0)
{
   auto readGate = [](TString fn, std::vector<double> &ax, std::vector<double> &ay) {
      if (gSystem->AccessPathName(fn)) { printf("\033[1;31m  %s NOT FOUND\033[0m\n", fn.Data()); return; }
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
   std::vector<double> ax, ay, bx, by;
   readGate(gateA, ax, ay);
   readGate(gateB, bx, by);
   if (ax.size() < 3 || bx.size() < 3) { printf("  need both gates\n"); return; }
   printf("  %s: %zu vertices\n  %s: %zu vertices\n", gateA.Data(), ax.size(), gateB.Data(), bx.size());

   // vertical extent of a closed polygon at x, from edge crossings
   auto extent = [](const std::vector<double> &X, const std::vector<double> &Y, double x, double &lo, double &hi) {
      lo = 1e9; hi = -1e9;
      size_t n = X.size();
      for (size_t i = 0, j = n - 1; i < n; j = i++) {
         double x1 = X[j], y1 = Y[j], x2 = X[i], y2 = Y[i];
         if ((x1 <= x && x2 > x) || (x2 <= x && x1 > x)) {
            double t = (x - x1) / (x2 - x1);
            double y = y1 + t * (y2 - y1);
            lo = std::min(lo, y); hi = std::max(hi, y);
         }
      }
      return hi > lo;
   };

   std::vector<double> vx, vlo, vhi;
   printf("\n  sqrt(dE/dx)   hand lo/hi        sim lo/hi         union lo/hi\n");
   for (int s = 0; s < nSlice; ++s) {
      double x = xMin + (s + 0.5) * (xMax - xMin) / nSlice;
      double alo, ahi, blo, bhi;
      bool inA = extent(ax, ay, x, alo, ahi);
      bool inB = extent(bx, by, x, blo, bhi);
      if (!inA && !inB) continue;
      double lo = inA ? alo : blo, hi = inA ? ahi : bhi;
      if (inB) { lo = std::min(lo, blo); hi = std::max(hi, bhi); }
      vx.push_back(x); vlo.push_back(lo); vhi.push_back(hi);
      if (s % 8 == 0)
         printf("   %6.2f     %s   %s   %.3f-%.3f\n", x,
                inA ? Form("%.3f-%.3f", alo, ahi) : "    --      ",
                inB ? Form("%.3f-%.3f", blo, bhi) : "    --      ", lo, hi);
   }
   if (vx.size() < 3) { printf("  nothing to write\n"); return; }

   std::vector<double> px, py;
   for (size_t i = 0; i < vx.size(); ++i) { px.push_back(vx[i]); py.push_back(vlo[i]); }
   for (int i = vx.size() - 1; i >= 0; --i) { px.push_back(vx[i]); py.push_back(vhi[i]); }

   FILE *f = fopen(outName, "w");
   if (!f) { printf("\033[1;31m  cannot write %s\033[0m\n", outName.Data()); return; }
   fprintf(f, "{\n    \"name\": \"proton_union\",\n    \"xaxis\": \"sqrtdEdx\",\n    \"yaxis\": \"brho\",\n");
   fprintf(f, "    \"vertices\": [\n");
   for (size_t i = 0; i < px.size(); ++i)
      fprintf(f, "        [\n            %.6f,\n            %.6f\n        ]%s\n", px[i], py[i],
              i + 1 < px.size() ? "," : "");
   fprintf(f, "    ]\n}\n");
   fclose(f);
   printf("\n  wrote %s (%zu vertices, sqrt(dE/dx) %.1f to %.1f)\n", outName.Data(), px.size(), vx.front(), vx.back());
}
