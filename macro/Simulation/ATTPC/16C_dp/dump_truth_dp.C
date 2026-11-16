/// @file dump_truth_dp.C  Generated proton of every event of a (d,p) simulation as text "ev tvz tth tke" (the keys.txt format of
/// sim_rebuild/trig/dump_keys.C; tvz = MCTrack start z x 10, same frame as sim_dump_dp.C's tvz). 2026-11-16, for the unbound efficiency.
///   root -l -b -q 'dump_truth_dp.C("p_ex39_sim.root","keys_ex39.txt")'
void dump_truth_dp(const char *sim, const char *out)
{
   TFile f(sim); auto *ts = (TTree *)f.Get("cbmsim"); TClonesArray *mc = nullptr;
   ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc);
   FILE *o = fopen(out, "w"); long n = 0;
   for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
      ts->GetEntry(i);
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (m && m->GetPdgCode() == 2212) { const double p = m->GetP() * 1000.; fprintf(o, "%lld %.4f %.4f %.4f\n", i, m->GetStartZ() * 10, std::acos(m->GetPz() / m->GetP()) * 57.29578, std::sqrt(p * p + 938.272 * 938.272) - 938.272); ++n; break; }
      }
   }
   fclose(o); printf("truth: %ld of %lld events -> %s\n", n, ts->GetEntries(), out);
}
