/// @file resistive_centroid_test.C
/// @brief Taubin-VALIDATED test of resistive per-ring charge centroiding. For each
///        track, group hits by ring and take the CHARGE-WEIGHTED centroid per ring
///        (one sub-pad point per ring), then fit with the FRAMEWORK's unbiased
///        Taubin+GN fit (NOT Kasa). Reports sigma AND median so bias is visible.
///        Compares raw-hit Taubin vs ring-centroid Taubin, on baseline (no PRF) and
///        PRF (charge-sharing) digi files.
/// Run: root -b -q 'resistive_centroid_test.C("./data/output_digi_prf0.5.root")'
double iqr(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void resistive_centroid_test(TString digiFile = "./data/output_digi_prf0.5.root", double Bfield = 4.0, double p0 = 374.9)
{
   gSystem->Load("libAtReconstruction.so");
   const double kRin = 62.9, kRout = 121.1; const int kNring = 16;
   std::vector<double> redge(kNring + 1); redge[0] = kRin;
   double dA = (kRout * kRout - kRin * kRin) / kNring;
   for (int i = 0; i < kNring; ++i) redge[i + 1] = std::sqrt(redge[i] * redge[i] + dA);
   auto ringOf = [&](double r) { for (int i = 0; i < kNring; ++i) if (r >= redge[i] && r < redge[i + 1]) return i; return (r < kRin) ? -1 : kNring - 1; };

   TFile f(digiFile); TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pe = new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent", &pe);

   // framework unbiased Taubin (+GN) fit on a set of positions; returns radius or -1
   auto taubin = [](const std::vector<std::pair<double, double>> &pts) -> double {
      if (pts.size() < 5) return -1;
      std::vector<AtHit> hitStore; hitStore.reserve(pts.size());
      std::vector<const AtHit *> hp; hp.reserve(pts.size());
      for (auto &q : pts) hitStore.emplace_back(0, ROOT::Math::XYZPoint(q.first, q.second, 0.0), 1.0);
      for (auto &h : hitStore) hp.push_back(&h);
      AtPatterns::AtPatternCircle2D c;
      c.AtPattern::FitPattern(hp, -1.0); // qThreshold=-1 -> UNWEIGHTED Taubin+GN
      double R = c.GetRadius();
      return (R > 0 && R < 1e5) ? R : -1;
   };

   std::vector<double> dpRaw, dpRing;
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      if (!pe->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pe->At(0))->GetTrackCand()) {
         std::vector<double> sx(kNring, 0), sy(kNring, 0), sq(kNring, 0);
         std::vector<std::pair<double, double>> raw;
         for (auto &h : tr.GetHitArray()) {
            auto p = h->GetPosition(); double q = std::max(1e-6, (double)h->GetCharge());
            double r = std::hypot(p.X(), p.Y()); int ri = ringOf(r); if (ri < 0) continue;
            sx[ri] += q * p.X(); sy[ri] += q * p.Y(); sq[ri] += q;
            raw.emplace_back(p.X(), p.Y());
         }
         std::vector<std::pair<double, double>> ringPts;
         for (int i = 0; i < kNring; ++i) if (sq[i] > 0) ringPts.emplace_back(sx[i] / sq[i], sy[i] / sq[i]);
         double Rr = taubin(raw), Rc = taubin(ringPts);
         if (Rr > 0) dpRaw.push_back((0.299792458 * Bfield * Rr - p0) / p0);
         if (Rc > 0) dpRing.push_back((0.299792458 * Bfield * Rc - p0) / p0);
      }
   }
   printf("\n===== Taubin-validated resistive centroiding: %s =====\n", digiFile.Data());
   printf("  raw hits, Taubin fit        : sigma=%5.1f%%  median=%+6.1f%%  (%zu)\n", 100 * iqr(dpRaw), 100 * med(dpRaw), dpRaw.size());
   printf("  per-RING centroid, Taubin   : sigma=%5.1f%%  median=%+6.1f%%  (%zu)\n", 100 * iqr(dpRing), 100 * med(dpRing), dpRing.size());
   printf("  (MC-truth floor ~2%%; honest baseline no-PRF ~19%%)\n");
}
