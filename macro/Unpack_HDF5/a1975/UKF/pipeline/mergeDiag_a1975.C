/// @file mergeDiag_a1975.C
/// @brief Print the geometric discriminators (centre dist, radius frac, endpoint gap,
/// chord cos, and a CIRCLE-RESIDUAL: how well track B's clusters lie on track A's PRA
/// circle) for the track pairs in selected events, to find what separates a genuine
/// fragmented track (should merge) from two tracks / a kink (should not).
///   root -l -b -q 'mergeDiag_a1975.C("run_0106", "{429,669,1001,2847}", "/mnt/f/a1975/reco/")'

void mergeDiag_a1975(TString runName = "run_0106", TString evList = "429,669,1001,2847",
                     TString recoDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   auto *f = TFile::Open(recoDir + runName + "_reco.root");
   auto *t = (TTree *)f->Get("cbmsim");
   auto *peArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &peArr);

   std::set<Long64_t> evs;
   for (auto *o : *evList.Tokenize(","))
      evs.insert(((TObjString *)o)->GetString().Atoll());

   auto axisEnds = [](const std::vector<AtHitCluster> &cl) {
      int iv = 0, ifar = 0; double rmin = 1e18, rmax = -1;
      for (int i = 0; i < (int)cl.size(); ++i) {
         auto p = cl[i].GetPosition(); double r = std::hypot(p.X(), p.Y());
         if (r < rmin) { rmin = r; iv = i; } if (r > rmax) { rmax = r; ifar = i; }
      }
      return std::pair{iv, ifar};
   };
   auto d3 = [](const AtHitCluster &a, const AtHitCluster &b) {
      auto pa = a.GetPosition(), pb = b.GetPosition();
      return std::sqrt(std::pow(pa.X()-pb.X(),2)+std::pow(pa.Y()-pb.Y(),2)+std::pow(pa.Z()-pb.Z(),2));
   };

   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      if (evs.find(i) == evs.end()) continue;
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      auto &trk = pe->GetTrackCand();
      std::cout << "\n=== ev " << i << " : " << trk.size() << " tracks ===\n";
      for (size_t a = 0; a < trk.size(); ++a)
         for (size_t b = a + 1; b < trk.size(); ++b) {
            auto &A = trk[a], &B = trk[b];
            auto *clA = A.GetHitClusterArray(); auto *clB = B.GetHitClusterArray();
            if (clA->empty() || clB->empty()) continue;
            auto cA = A.GetGeoCenter(); auto cB = B.GetGeoCenter();
            double cd = std::hypot(cA.first-cB.first, cA.second-cB.second);
            double rA = A.GetGeoRadius(), rB = B.GetGeoRadius();
            double rfrac = std::max(rA,rB) > 0 ? std::abs(rA-rB)/std::max(rA,rB) : 0;
            auto [ivA, ifA] = axisEnds(*clA); auto [ivB, ifB] = axisEnds(*clB);
            double gap = std::min({d3((*clA)[ifA],(*clB)[ivB]), d3((*clA)[ifA],(*clB)[ifB]), d3((*clA)[ivA],(*clB)[ifB])});
            double vtxgap = d3((*clA)[ivA],(*clB)[ivB]);
            auto dA = (*clA)[ifA].GetPosition()-(*clA)[ivA].GetPosition();
            auto dB = (*clB)[ifB].GetPosition()-(*clB)[ivB].GetPosition();
            double cosang = (dA.R()>0&&dB.R()>0)? dA.Unit().Dot(dB.Unit()) : 0;
            // circle residual: RMS distance of B's clusters from A's PRA circle (centre cA, radius rA)
            double res = 0; int nb = 0;
            for (auto &c : *clB) { auto p = c.GetPosition();
               res += std::pow(std::hypot(p.X()-cA.first, p.Y()-cA.second) - rA, 2); nb++; }
            res = nb ? std::sqrt(res/nb) : -1;
            std::cout << "  pair(" << a << "," << b << "): centreDist=" << cd << "mm radiusFrac=" << rfrac
                      << " gap=" << gap << "mm vtxGap=" << vtxgap << "mm chordCos=" << cosang
                      << " B-on-A-circle RMS=" << res << "mm  (rA=" << rA << " rB=" << rB << ")\n";
         }
   }
}
