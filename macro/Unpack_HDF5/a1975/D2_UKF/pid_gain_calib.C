/// @file pid_gain_calib.C
/// @brief Write the per-run dE/dx GAIN TABLE (AtTools::AtPIDGainTable format) for a production.
///
/// For every run: the position of the dominant PID band on a FIXED-KINEMATICS slice, separately for
/// the Spyral and the 3D estimator. On a1975 D2 the slice theta_lab 78-86 deg, Brho 0.15-0.35 T m is
/// dominated by elastic deuterons, so momentum does not enter and only the gain moves the band.
/// Factor = exp(ref - position) multiplies sqrt(dE/dx); the reference run gets 1.000.
/// Position = charge-weighted mean of a +-0.065 window around the mode of ln(sqrt(dE/dx)*Brho)
/// (robust against the tails); the run's slice count and the window's RMS are written as comments
/// so a thinly populated run is visible in the table itself.
/// Input: PID trees from pid_variants_dp.C (Spyral sqa/br, 3D s0/b0), file <prefix><run>.root.
///
///   root -l -b -q 'pid_gain_calib.C("/mnt/f/a1975/pid3d_gaincal/pv_","pid/gain_d2_dv1104.txt","0016")'
void pid_gain_calib(TString prefix, TString tableOut, TString refRun = "0016",
                    TString slice = "va&&v3&&sp3<6&&nh>=30&&th3>78&&th3<86", double brLo = 0.15, double brHi = 0.35)
{
   // runs = every <prefix><digits>.root
   TString dir = gSystem->GetDirName(prefix), base = gSystem->BaseName(prefix);
   std::vector<TString> runs;
   void *d = gSystem->OpenDirectory(dir);
   while (const char *e = gSystem->GetDirEntry(d)) {
      TString s(e);
      if (s.BeginsWith(base) && s.EndsWith(".root")) {
         TString r = s(base.Length(), s.Length() - base.Length() - 5);
         if (r.IsDigit()) runs.push_back(r);
      }
   }
   gSystem->FreeDirectory(d);
   std::sort(runs.begin(), runs.end());
   struct Pos { double m; int n; double rms; };
   auto band = [&](TTree *w, const char *sq, const char *br) {
      TH1D h("h", "", 400, -2, 6);
      w->Draw(Form("log(%s*%s)>>h", sq, br), slice + Form("&&%s>%g&&%s<%g&&%s>0", br, brLo, br, brHi, sq), "goff");
      int ib = h.GetMaximumBin();
      double s = 0, s2 = 0, sw = 0;
      for (int b = ib - 6; b <= ib + 6; ++b) {
         double x = h.GetBinCenter(b), c = h.GetBinContent(b);
         s += x * c; s2 += x * x * c; sw += c;
      }
      Pos p{sw > 0 ? s / sw : 0, (int)h.GetEntries(), sw > 0 ? std::sqrt(std::max(0.0, s2 / sw - (s / sw) * (s / sw))) : 0};
      return p;
   };
   std::map<TString, std::pair<Pos, Pos>> pos;
   for (auto &r : runs) {
      TFile f(prefix + r + ".root");
      auto *w = (TTree *)f.Get("w");
      if (!w) continue;
      pos[r] = {band(w, "sqa", "br"), band(w, "s0", "b0")};
   }
   if (!pos.count(refRun)) { printf("ERROR: reference run %s not found\n", refRun.Data()); return; }
   double refS = pos[refRun].first.m, ref3 = pos[refRun].second.m;
   std::ofstream out(tableOut.Data());
   out << "# a1975 D2 per-run dE/dx gain table (AtTools::AtPIDGainTable): run gainSpyral gain3D\n"
       << "# factor multiplies sqrt(dE/dx); reference run " << refRun.Data() << " = 1.000\n"
       << "# made by pid_gain_calib.C on " << TDatime().AsSQLString() << " from " << prefix.Data() << "<run>.root\n"
       << "# band = mode of ln(sqrt(dE/dx)*Brho) on: " << slice.Data() << " && Brho " << brLo << "-" << brHi << "\n"
       << "# columns after '#': slice entries, band position and window RMS for Spyral | 3D\n";
   printf("%-6s %8s %8s   %7s %7s %6s\n", "run", "gSpyral", "g3D", "posS", "pos3D", "n");
   for (auto &[r, p] : pos) {
      double gS = std::exp(refS - p.first.m), g3 = std::exp(ref3 - p.second.m);
      out << Form("%d %.4f %.4f   # n=%d posS=%.4f rmsS=%.4f | pos3D=%.4f rms3D=%.4f\n", r.Atoi(), gS, g3, p.second.n,
                  p.first.m, p.first.rms, p.second.m, p.second.rms);
      printf("%-6s %8.4f %8.4f   %7.4f %7.4f %6d%s\n", r.Data(), gS, g3, p.first.m, p.second.m, p.second.n,
             p.second.n < 500 ? "  <-- thin" : "");
   }
   printf("wrote %s (%zu runs)\n", tableOut.Data(), pos.size());
}
