/// @file ring_cluster_test.C
/// @brief Test "cluster per ring": on a resistive pad plane the charge spreads over
///        several AZIMUTHAL pads within the ring the track crosses. Group each track's
///        hits by ring and take the CHARGE-WEIGHTED centroid per ring -> one sub-pad
///        hit per ring -> fit a circle -> momentum. If sigma_p/p drops toward the MC
///        floor (2%), the resistive sharing IS exploitable with ring-centroiding.
/// Run: root -b -q 'ring_cluster_test.C("./data/output_digi_prf0.5.root")'
/// WARNING: the Kasa algebraic circle fit here is BIASED (~-79% median) on shallow
/// noisy arcs — deceptively tight sigma, wildly wrong radius. Use
/// resistive_centroid_test.C (framework Taubin, unbiased) for real momentum numbers.
/// Kept to DOCUMENT the Kasa-vs-Taubin gap. ALWAYS check the MEDIAN, not just sigma.
double iqr(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void ring_cluster_test(TString digiFile = "./data/output_digi_prf0.5.root", double Bfield = 4.0, double p0 = 374.9)
{
   gSystem->Load("libAtReconstruction.so");
   const double kRin = 62.9, kRout = 121.1; const int kNring = 16;
   std::vector<double> redge(kNring + 1); redge[0] = kRin;
   double dA = (kRout * kRout - kRin * kRin) / kNring;
   for (int i = 0; i < kNring; ++i) redge[i + 1] = std::sqrt(redge[i] * redge[i] + dA);
   auto ringOf = [&](double r) { for (int i = 0; i < kNring; ++i) if (r >= redge[i] && r < redge[i + 1]) return i; return (r < kRin) ? -1 : kNring - 1; };

   TFile f(digiFile); TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pe = new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent", &pe);

   auto kasa = [](std::vector<std::pair<double, double>> &pp, double &R) -> bool {
      int N = pp.size(); if (N < 5) return false;
      double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
      for (auto &q : pp) { double x = q.first, y = q.second, z = x * x + y * y;
         Sx += x; Sy += y; Sxx += x * x; Syy += y * y; Sxy += x * y; Sxz += x * z; Syz += y * z; Sz += z; }
      double a11 = Sxx, a12 = Sxy, a13 = Sx, a22 = Syy, a23 = Sy, a33 = N, b1 = -Sxz, b2 = -Syz, b3 = -Sz;
      double det = a11 * (a22 * a33 - a23 * a23) - a12 * (a12 * a33 - a23 * a13) + a13 * (a12 * a23 - a22 * a13);
      if (std::abs(det) < 1e-9) return false;
      double D = (b1 * (a22 * a33 - a23 * a23) - a12 * (b2 * a33 - a23 * b3) + a13 * (b2 * a23 - a22 * b3)) / det;
      double E = (a11 * (b2 * a33 - a23 * b3) - b1 * (a12 * a33 - a23 * a13) + a13 * (a12 * b3 - b2 * a13)) / det;
      double F = (a11 * (a22 * b3 - b2 * a23) - a12 * (a12 * b3 - b2 * a13) + b1 * (a12 * a23 - a22 * a13)) / det;
      double cx = -D / 2, cy = -E / 2, arg = cx * cx + cy * cy - F; if (arg <= 0) return false;
      R = std::sqrt(arg); return R > 0 && R < 1e5;
   };

   std::vector<double> dpRing, dpAll, dpPRA, dpTaubin;
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      if (!pe->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pe->At(0))->GetTrackCand()) {
         double Rpra = tr.GetGeoRadius();
         if (Rpra > 0 && Rpra < 1e5) dpPRA.push_back((0.299792458 * Bfield * Rpra - p0) / p0);
         // ring-grouped charge-weighted centroids
         std::vector<double> sx(kNring, 0), sy(kNring, 0), sq(kNring, 0);
         std::vector<std::pair<double, double>> allHits;
         std::vector<std::array<double,3>> allHitsQ;
         for (auto &h : tr.GetHitArray()) {
            auto p = h->GetPosition(); double q = std::max(1e-6, (double)h->GetCharge());
            double r = std::hypot(p.X(), p.Y()); int ri = ringOf(r); if (ri < 0) continue;
            sx[ri] += q * p.X(); sy[ri] += q * p.Y(); sq[ri] += q;
            allHits.emplace_back(p.X(), p.Y());
            allHitsQ.push_back({p.X(), p.Y(), q});
         }
         std::vector<std::pair<double, double>> ringPts;
         for (int i = 0; i < kNring; ++i) if (sq[i] > 0) ringPts.emplace_back(sx[i] / sq[i], sy[i] / sq[i]);
         double R;
         if (kasa(ringPts, R)) dpRing.push_back((0.299792458 * Bfield * R - p0) / p0);
         if (kasa(allHits, R)) dpAll.push_back((0.299792458 * Bfield * R - p0) / p0);
         // Kasa but CHARGE-WEIGHTED (each hit weighted by its charge), to test whether
         // the framework's charge-weighting is what degrades the radius.
         if (allHitsQ.size() >= 5) {
            double Sw = 0, Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
            for (auto &hq : allHitsQ) { double x = hq[0], y = hq[1], w = hq[2], z = x * x + y * y;
               Sw += w; Sx += w * x; Sy += w * y; Sxx += w * x * x; Syy += w * y * y; Sxy += w * x * y; Sxz += w * x * z; Syz += w * y * z; Sz += w * z; }
            double a11 = Sxx, a12 = Sxy, a13 = Sx, a22 = Syy, a23 = Sy, a33 = Sw, b1 = -Sxz, b2 = -Syz, b3 = -Sz;
            double det = a11 * (a22 * a33 - a23 * a23) - a12 * (a12 * a33 - a23 * a13) + a13 * (a12 * a23 - a22 * a13);
            if (std::abs(det) > 1e-9) {
               double D = (b1 * (a22 * a33 - a23 * a23) - a12 * (b2 * a33 - a23 * b3) + a13 * (b2 * a23 - a22 * b3)) / det;
               double E = (a11 * (b2 * a33 - a23 * b3) - b1 * (a12 * a33 - a23 * a13) + a13 * (a12 * b3 - b2 * a13)) / det;
               double F = (a11 * (a22 * b3 - b2 * a23) - a12 * (a12 * b3 - b2 * a13) + b1 * (a12 * a23 - a22 * a13)) / det;
               double cx = -D / 2, cy = -E / 2, arg = cx * cx + cy * cy - F;
               if (arg > 0) { double Rt = std::sqrt(arg); if (Rt > 0 && Rt < 1e5) dpTaubin.push_back((0.299792458 * Bfield * Rt - p0) / p0); }
            }
         }
      }
   }
   printf("\n===== ring-centroid circle fit on %s =====\n", digiFile.Data());
   printf("  per-RING charge centroid : sigma_p/p = %.1f%%  (%zu tracks, ~%d ring-hits each)\n",
          100 * iqr(dpRing), dpRing.size(), kNring);
   printf("  all raw pad hits (Kasa)  : sigma=%.1f%%  MEDIAN=%+.1f%%  (%zu)\n", 100*iqr(dpAll), 100*med(dpAll), dpAll.size());
   printf("  Kasa CHARGE-WEIGHTED (all hits): sigma_p/p = %.1f%%  (%zu)\n", 100 * iqr(dpTaubin), dpTaubin.size());
   printf("  PRA GetGeoRadius (frmwk) : sigma_p/p = %.1f%%  (%zu tracks)  median %+.1f%%\n", 100 * iqr(dpPRA), dpPRA.size(), 100 * med(dpPRA));
   printf("  (full UKF/genfit pipeline ~19%%; MC-truth floor ~2%%)\n");
}
