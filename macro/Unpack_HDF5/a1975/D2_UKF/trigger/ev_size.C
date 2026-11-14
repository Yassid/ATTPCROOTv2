// total AtEventH hits for listed (run, entry) -- environment of real (d,p) candidates
void ev_size(const char *in, const char *out) {
   gErrorIgnoreLevel = kError;
   std::ifstream f(in); std::map<int, std::vector<std::pair<long,int>>> L; int r, e, nh; while (f >> r >> e >> nh) L[r].push_back({e, nh});
   FILE *o = fopen(out, "w");
   for (auto &[run, v] : L) {
      TFile fd(Form("/mnt/f/a1975/reco_d2_dv1104/run_%04d_multifit_reco.root", run)); auto *td = (TTree *)fd.Get("cbmsim"); if (!td) continue;
      TClonesArray *hd = nullptr; td->SetBranchStatus("*", 0); td->SetBranchStatus("AtEventH*", 1); td->SetBranchAddress("AtEventH", &hd);
      for (auto &[ev, n] : v) { if (ev >= td->GetEntries()) continue; td->GetEntry(ev); auto *E = (AtEvent *)hd->At(0); fprintf(o, "%d %ld %d %zu\n", run, ev, n, E ? E->GetHits().size() : 0); }
   }
   fclose(o);
}
