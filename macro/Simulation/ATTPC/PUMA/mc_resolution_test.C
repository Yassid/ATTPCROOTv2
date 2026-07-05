/// @file mc_resolution_test.C
/// @brief Momentum-resolution FLOOR: fit a circle directly to the TRUTH MC points
///        of each primary pion (sigma_x ~ 0, no digitization/PSA/clustering noise),
///        same annular lever arm as the digi. If this is ~few %, the ~19% digi
///        resolution is set by the point resolution sigma_x (pads+diffusion), not
///        the lever arm; if it is also ~19%, it is the geometry/scattering floor.
///        "every5" downsamples MC points to 1-in-5 (the user's suggestion).
/// Run: root -b -q mc_resolution_test.C
///      root -b -q 'mc_resolution_test.C(5)'   // use every 5th MC point
double iqr(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void mc_resolution_test(int stride = 1, double Bfield = 4.0, double p0 = 374.9,
                        TString simFile = "./data/attpcsim.root")
{
   gSystem->Load("libAtReconstruction.so");
   TFile f(simFile); TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *mcP = new TClonesArray("AtMCPoint"); t->SetBranchAddress("AtTpcPoint", &mcP);

   std::vector<double> dp;           // (p_fit - p0)/p0 from the MC circle fit
   std::vector<double> nptsV, Lv, sagResid;
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      // gather MC points per primary trackID (0,1)
      std::map<int, std::vector<std::pair<double, double>>> pts;
      for (int k = 0; k < mcP->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcP->At(k);
         int tid = mp->GetTrackID();
         if (tid != 0 && tid != 1) continue; // primaries only
         pts[tid].emplace_back(mp->GetX() * 10, mp->GetY() * 10); // mm
      }
      for (auto &kv : pts) {
         auto &v = kv.second;
         // optional downsample to every 'stride'th point
         std::vector<std::pair<double, double>> pp;
         for (size_t i = 0; i < v.size(); i += stride) pp.push_back(v[i]);
         int N = pp.size();
         if (N < 5) continue;
         // Kasa algebraic circle fit: minimize sum (x^2+y^2 + D x + E y + F)^2
         double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
         for (auto &q : pp) { double x = q.first, y = q.second, z = x * x + y * y;
            Sx += x; Sy += y; Sxx += x * x; Syy += y * y; Sxy += x * y; Sxz += x * z; Syz += y * z; Sz += z; }
         // normal equations for D,E,F (n=N)
         double a11 = Sxx, a12 = Sxy, a13 = Sx, a21 = Sxy, a22 = Syy, a23 = Sy, a31 = Sx, a32 = Sy, a33 = N;
         double b1 = -Sxz, b2 = -Syz, b3 = -Sz;
         // solve 3x3 (Cramer)
         double det = a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
         if (std::abs(det) < 1e-9) continue;
         double D = (b1 * (a22 * a33 - a23 * a32) - a12 * (b2 * a33 - a23 * b3) + a13 * (b2 * a32 - a22 * b3)) / det;
         double E = (a11 * (b2 * a33 - a23 * b3) - b1 * (a21 * a33 - a23 * a31) + a13 * (a21 * b3 - b2 * a31)) / det;
         double F = (a11 * (a22 * b3 - b2 * a32) - a12 * (a21 * b3 - b2 * a31) + b1 * (a21 * a32 - a22 * a31)) / det;
         double cx = -D / 2, cy = -E / 2, R = std::sqrt(cx * cx + cy * cy - F);
         if (!(R > 0) || R > 1e5) continue;
         double p = 0.299792458 * Bfield * R; // MeV/c (R in mm)
         dp.push_back((p - p0) / p0);
         nptsV.push_back(N);
         // chord (lever arm) of the MC points
         double x0 = pp.front().first, y0 = pp.front().second, x1 = pp.back().first, y1 = pp.back().second;
         Lv.push_back(std::hypot(x1 - x0, y1 - y0));
      }
   }
   printf("\n===== MC-truth circle-fit resolution (stride=%d, %zu primaries) =====\n", stride, dp.size());
   printf("  lever arm L (MC chord) : %.1f mm     hits used/track : %.0f\n", med(Lv), med(nptsV));
   printf("  sigma_p/p (MC points)  : %.1f%%  (median %+.1f%%)\n", 100 * iqr(dp), 100 * med(dp));
   printf("  --> digi sigma_p/p was ~19%% (UKF).  If MC << digi, sigma_x (pads+diffusion) is the limit;\n");
   printf("      if MC ~ digi, the lever arm / multiple scattering is the floor.\n");
}
