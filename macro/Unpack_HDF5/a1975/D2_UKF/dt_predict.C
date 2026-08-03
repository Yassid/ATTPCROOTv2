/// @file dt_predict.C
/// @brief What the p_z rescaling predicts the 16C(d,t)15C levels do as the drift velocity moves.
///
/// This exists to be WRONG in public. dt_dvscan.C sweeps dv by remapping already-fitted tracks
/// (p_T fixed by the pad plane, p_z scaled by k = dv_true/dv_used) instead of re-running the
/// reconstruction, so it ignores how a changed z would have altered clustering, seeding and the
/// fit's own convergence. Printing the predicted level positions BEFORE the re-reco at 1.25 lands
/// makes the shortcut testable: if the re-reco reproduces this table the remap can be trusted for
/// scanning, and if it does not, every dv number that came out of the scan has to be dropped.
///
/// Two observables, chosen for what thin statistics can actually support:
///   the g.s. complex -- the strongest structure in the spectrum, so it survives a 3-run subset
///   the 3.103 level  -- isolated and therefore the honest ruler, but needs the full 47 runs
///
///   root -b -q 'dt_predict.C("/mnt/f/a1975/dt_kin_full.root")'

static double omP(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// smoothed tallest bin inside a window, plus the windowed mean -- no seeded fit, which is what
/// fell apart on the thin (d,t) peak when the scale moved underneath it
static void peakIn(TH1F *h, double lo, double hi, double &pk, double &mean)
{
   TH1F s(*h);
   s.SetDirectory(nullptr);
   s.Smooth(2);
   int bm = 0;
   double mx = -1, sw = 0, sx = 0;
   for (int b = 1; b <= s.GetNbinsX(); ++b) {
      double c = s.GetBinCenter(b);
      if (c < lo || c > hi)
         continue;
      if (s.GetBinContent(b) > mx) { mx = s.GetBinContent(b); bm = b; }
      sw += h->GetBinContent(b);
      sx += h->GetBinContent(b) * c;
   }
   pk = bm ? s.GetBinCenter(bm) : std::nan("");
   mean = sw > 0 ? sx / sw : std::nan("");
}

void dt_predict(TString cache = "/mnt/f/a1975/dt_kin_full.root", double dvUsed = 1.15, double Ebeam = 180.0,
                double icLo = 900, double icHi = 1300, double chi2max = 5, double keMin = 5, double vzLo = 50,
                double vzHi = 700, TString runsOnly = "")
{
   const double u = 931.49401;
   const double m1 = 16.0147013 * u, m2 = 2.0135532 * u, m3 = 3.01550072 * u, m4 = 15.0105993 * u;

   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   float ke, th, vz = 0, c2, ic = -1, run = -1;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("vertexz", &vz);
   t->SetBranchAddress("chi2ndf", &c2);
   if (t->GetBranch("ic")) t->SetBranchAddress("ic", &ic);
   if (t->GetBranch("run")) t->SetBranchAddress("run", &run);

   std::set<int> keep;
   if (runsOnly.Length()) {
      TObjArray *a = runsOnly.Tokenize(",");
      for (int i = 0; i < a->GetEntries(); ++i) keep.insert(((TObjString *)a->At(i))->GetString().Atoi());
   }

   std::vector<double> pT, pz, vzv;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (c2 > chi2max || ke <= 0) continue;
      if (icLo > 0 && (ic < icLo || ic > icHi)) continue;
      if (!keep.empty() && !keep.count((int)(run + 0.5))) continue;
      double p = std::sqrt((ke + m3) * (ke + m3) - m3 * m3), a = th * TMath::DegToRad();
      pT.push_back(p * std::sin(a));
      pz.push_back(p * std::cos(a));
      vzv.push_back(vz);
   }
   printf("\n=== dt_predict: %zu tracks%s, Ebeam %.0f, dv in use %.2f ===\n", pT.size(),
          keep.empty() ? "" : Form(" (runs %s only)", runsOnly.Data()), Ebeam, dvUsed);
   printf("%-8s %7s %10s %10s %12s %11s\n", "dv", "N", "g.s. peak", "g.s. mean", "3.103 peak", "shift[MeV]");
   std::vector<double> ref;

   for (double dv : {1.15, 1.20, 1.25, 1.30}) {
      const double k = dv / dvUsed;
      TH1F h("h", "", 90, -3, 6);
      h.SetDirectory(nullptr);
      long n = 0;
      for (size_t i = 0; i < pT.size(); ++i) {
         double qz = k * pz[i], qT = pT[i], z = k * vzv[i];
         if (z < vzLo || z > vzHi) continue;
         double kep = std::sqrt(qT * qT + qz * qz + m3 * m3) - m3;
         if (kep < keMin) continue;
         double thp = std::atan2(qT, qz);
         double E1 = Ebeam + m1, E3 = kep + m3;
         double s = m1 * m1 + m2 * m2 + 2 * m2 * E1, uu = m2 * m2 + m3 * m3 - 2 * m2 * E3;
         double ar = (std::cos(thp) * omP(s, m1 * m1, m2 * m2) * omP(uu, m2 * m2, m3 * m3) -
                      (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) / (2 * m2 * m2) + s + uu - m2 * m2;
         if (ar < 0) continue;
         h.Fill(std::sqrt(ar) - m4);
         ++n;
      }
      double gp, gm, tp, tm;
      peakIn(&h, -1.5, 1.5, gp, gm);
      peakIn(&h, 2.0, 4.5, tp, tm);
      // How far the WHOLE spectrum moved, by cross-correlation against the dv = 1.15 shape.
      // With ~3.7k events the tallest-bin position hops between structures and a windowed mean
      // is pinned by its own window as soon as the shift approaches the window width, so both
      // are useless here; the correlation shift uses every bin and degrades gracefully.
      if (ref.empty()) {
         ref.assign(h.GetNbinsX() + 2, 0.0);
         for (int b = 1; b <= h.GetNbinsX(); ++b) ref[b] = h.GetBinContent(b);
      }
      double bestS = 0, bestC = -1;
      const double bw = h.GetBinWidth(1);
      for (int sh = -40; sh <= 40; ++sh) {
         double c = 0, na = 0, nb = 0;
         for (int b = 1; b <= h.GetNbinsX(); ++b) {
            int bs = b + sh;
            if (bs < 1 || bs > h.GetNbinsX()) continue;
            c += ref[b] * h.GetBinContent(bs);
            na += ref[b] * ref[b];
            nb += h.GetBinContent(bs) * h.GetBinContent(bs);
         }
         double cc = (na > 0 && nb > 0) ? c / std::sqrt(na * nb) : 0;
         if (cc > bestC) { bestC = cc; bestS = sh * bw; }
      }
      printf("%-8.2f %7ld %10.3f %10.3f %12.3f %11.2f\n", dv, n, gp, gm, tp, bestS);
   }
   printf("\nshift is of the whole spectrum vs the dv = 1.15 shape (positive = E_x moves up)\n");
   printf("the re-reco at 1.25 has to reproduce the 1.25 row, or the remap is not usable\n");
}
