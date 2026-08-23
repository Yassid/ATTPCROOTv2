/// @file score_gate_pt_sim.C
/// @brief Score a hand-drawn (p,t) PID gate against MC truth, PER STATE.
///
/// The data plane cannot do this: there is no truth there, so a gate is trusted rather than
/// measured. Here every track carries `istriton` (species A=3, Z=1, majority vote over its hits),
/// so the polygon can be given the two numbers that matter:
///
///   EFFICIENCY  kept tritons / all tritons        -- what fraction of the signal survives
///   PURITY      kept tritons / all kept           -- how much of what survives is signal
///
/// PER STATE IS THE POINT. A gate with 95% efficiency everywhere costs statistics and nothing
/// else; a gate with 95% on the g.s. and 80% on the 7.012 silently rescales the level ratios,
/// which is the quantity the analysis exists to measure. On (p,d) the data-plane polygon applied
/// to simulation ran 75.7 / 77.2 / 80.5 / 85.9 across four states -- a 10-point spread that would
/// have distorted exactly that.
///
///   root -b -q 'score_gate_pt_sim.C("triton_pt_sim.json","data/pid_plane_pt_sim.root")'
void score_gate_pt_sim(TString gateFile = "triton_pt_sim.json",
                       TString cache = "data/pid_plane_pt_sim.root")
{
   // --- read the polygon straight from the JSON: no loader, no chance of the wrong gate ---
   std::ifstream in(gateFile.Data());
   if (!in) { printf("cannot open %s\n", gateFile.Data()); return; }
   std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   size_t vp = all.find("\"vertices\"");
   if (vp == std::string::npos) { printf("no vertices in %s\n", gateFile.Data()); return; }
   std::vector<double> vx, vy;
   size_t p = all.find('[', vp + 10);
   while (true) {
      size_t a = all.find('[', p + 1);
      if (a == std::string::npos) break;
      size_t b = all.find(']', a);
      if (b == std::string::npos) break;
      double x, y;
      if (sscanf(all.substr(a + 1, b - a - 1).c_str(), "%lf , %lf", &x, &y) == 2) { vx.push_back(x); vy.push_back(y); }
      p = b;
      size_t nxt = all.find_first_not_of(" \n\r\t,", b + 1);
      if (nxt == std::string::npos || all[nxt] == ']') break;
   }
   if (vx.size() < 3) { printf("only %zu vertices parsed\n", vx.size()); return; }
   TCutG cut("gate", vx.size());
   for (size_t i = 0; i < vx.size(); ++i) cut.SetPoint(i, vx[i], vy[i]);
   printf("gate %s: %zu vertices\n\n", gateFile.Data(), vx.size());

   TFile *F = TFile::Open(cache);
   TTree *t = F ? (TTree *)F->Get("pts") : nullptr;
   if (!t) { printf("no tree pts in %s\n", cache.Data()); return; }
   float x, y; int istri, st;
   t->SetBranchAddress("sqrtdedx", &x);
   t->SetBranchAddress("brho", &y);
   t->SetBranchAddress("istriton", &istri);
   t->SetBranchAddress("state", &st);

   const char *nm[5] = {"gs  0.000", "ex1 0.740", "ex2 3.103", "ex3 4.657", "ex4 7.012"};
   long tri[5] = {0}, keptTri[5] = {0}, kept[5] = {0}, tot[5] = {0};
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (st < 0 || st > 4) continue;
      ++tot[st];
      if (istri) ++tri[st];
      if (cut.IsInside(x, y)) { ++kept[st]; if (istri) ++keptTri[st]; }
   }

   printf("%-11s %8s %8s %9s %9s %9s\n", "state", "tracks", "tritons", "kept", "eff", "purity");
   long T = 0, KT = 0, K = 0;
   double emin = 1e9, emax = -1e9;
   for (int s = 0; s < 5; ++s) {
      if (!tot[s]) continue;
      double eff = tri[s] ? 100.0 * keptTri[s] / tri[s] : 0.0;
      double pur = kept[s] ? 100.0 * keptTri[s] / kept[s] : 0.0;
      printf("%-11s %8ld %8ld %9ld %8.1f%% %8.1f%%\n", nm[s], tot[s], tri[s], kept[s], eff, pur);
      T += tri[s]; KT += keptTri[s]; K += kept[s];
      emin = std::min(emin, eff); emax = std::max(emax, eff);
   }
   printf("%-11s %8s %8ld %9ld %8.1f%% %8.1f%%\n", "ALL", "", T, K,
          T ? 100.0 * KT / T : 0.0, K ? 100.0 * KT / K : 0.0);
   printf("\nefficiency spread across states: %.1f points (%.1f%% to %.1f%%)\n", emax - emin, emin, emax);
   printf(emax - emin < 2.0
             ? "  -> flat. The gate does not rescale one level against another.\n"
             : "  -> STATE DEPENDENT. This biases the level ratios and must be carried as a\n"
               "     correction or the gate redrawn; it is not a pure statistics loss.\n");
   F->Close();
}
