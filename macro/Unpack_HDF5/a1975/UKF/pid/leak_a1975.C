/// @file leak_a1975.C
/// @brief Quantify deuteron leakage into the proton PID gate, numerically
/// (no images). For the IC+quality clean sample it builds, per theta bin, the
/// median brho of proton-gated and deuteron-gated tracks (the two kinematic
/// lines), and counts how many proton-gate tracks sit closer to the deuteron
/// line than the proton line (= leakage).
void leak_a1975(TString cacheTag = "combined", TString protonGate = "proton_pid.json",
                TString deuteronGate = "deuteron_pid.json", Double_t icMin = 950, Double_t icMax = 1350,
                Int_t minClusters = 15, Double_t maxVertexR = 40, Double_t polarMin = 10, Double_t polarMax = 170)
{
   gSystem->Load("libAtTools.so");
   auto pP = AtTools::AtParticleID::LoadJSON(protonGate.Data());
   auto pD = AtTools::AtParticleID::LoadJSON(deuteronGate.Data());

   TFile *fo = TFile::Open(cacheTag + "_pidobs.root");
   TNtuple *nt = (TNtuple *)fo->Get("pidobs");
   float dedx, brho, ncl, ic, polar, vtxr, vtxz;
   nt->SetBranchAddress("dedx", &dedx);
   nt->SetBranchAddress("brho", &brho);
   nt->SetBranchAddress("ncl", &ncl);
   nt->SetBranchAddress("ic", &ic);
   nt->SetBranchAddress("polar", &polar);
   nt->SetBranchAddress("vtxr", &vtxr);
   nt->SetBranchAddress("vtxz", &vtxz);

   const int NB = 18; // 10-degree theta bins
   std::vector<std::vector<double>> pBr(NB), dBr(NB);
   for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
      nt->GetEntry(i);
      if (ncl < minClusters || vtxr > maxVertexR || polar < polarMin || polar > polarMax)
         continue;
      if (ic < icMin || ic > icMax)
         continue;
      int b = (int)(polar / 10.0);
      if (b < 0 || b >= NB)
         continue;
      if (pP.IsInside(dedx, brho))
         pBr[b].push_back(brho);
      if (pD.IsInside(dedx, brho))
         dBr[b].push_back(brho);
   }
   auto med = [](std::vector<double> &v) {
      if (v.empty())
         return -1.0;
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };
   printf("theta[deg]  Np   medBrho_p   Nd   medBrho_d   leak%%(p-gate near d-line)\n");
   long totP = 0, totLeak = 0;
   for (int b = 0; b < NB; ++b) {
      double mp = med(pBr[b]), mdv = med(dBr[b]);
      long leak = 0;
      if (mp > 0 && mdv > 0) {
         double mid = 0.5 * (mp + mdv); // boundary between the two kinematic lines
         for (double x : pBr[b])
            if (x > mid)
               leak++; // proton-gate track sitting on the deuteron side
      }
      totP += pBr[b].size();
      totLeak += leak;
      if (!pBr[b].empty() || !dBr[b].empty())
         printf("  %3d-%-3d  %4zu   %7.3f   %4zu   %7.3f   %5.1f\n", b * 10, b * 10 + 10, pBr[b].size(), mp,
                dBr[b].size(), mdv, pBr[b].empty() ? 0.0 : 100.0 * leak / pBr[b].size());
   }
   printf("\nproton-gate total %ld, on deuteron side %ld -> overall leakage %.1f%%\n", totP, totLeak,
          totP ? 100.0 * totLeak / totP : 0);
}
