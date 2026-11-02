/// @file order_score_pra_dp.C
/// @brief Scores AtPRA::SetRemoveStrayHits end to end (2026-11-02): pattern recognition re-run with the option (rerun_pra_dp.C), fitted by order_fit_dp.C
///        mode 0, against production clusters (mode 0) and the offline guard test (mode 9), same tracks matched by (event, trackID).
///        what = "d": sim deuterons, "t": sim tritons (GOOD |ln(Brho/true)| < 0.2 over all tracks); "data": deuteron band of run 0026, first 8000 events
///        (GOOD |Ex| < 2 MeV AND |ln(Brho_fit/bA)| < 0.5 over all tracks, the repaired referee of order_score5_dp.C).
///   root -l -b -q 'order_score_pra_dp.C("data")'
namespace ordpra {
struct R { int ok{0}; double ke{-1}, th{-1}; };
double exEl(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ed = ke + md, pd = std::sqrt(Ed * Ed - md * md);
   const double E4 = Eb + md - Ed, px = -pd * std::sin(th), pz = pb - pd * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - mC;
}
/// (event, trackID) of each row of order_fit_dp.C run on the production reco file
void rowMap(TString reco, TString sim, int pdg, Long64_t maxEv, std::vector<std::pair<Long64_t, int>> &rows)
{
   TFile f(reco);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr, *mc = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   TFile *fs = sim.Length() ? TFile::Open(sim) : nullptr;
   TTree *ts = fs ? (TTree *)fs->Get("cbmsim") : nullptr;
   if (ts) { ts->SetBranchStatus("*", 0); ts->SetBranchStatus("MCTrack*", 1); ts->SetBranchAddress("MCTrack", &mc); }
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      auto &c = ((AtPatternEvent *)pa->At(0))->GetTrackCand();
      if (ts) {
         int nl = 0;
         for (auto &tr : c) nl += tr.GetHitArray().size() >= 30;
         if (nl != 1) continue;
         ts->GetEntry(ev);
         bool has = false;
         for (int j = 0; j < mc->GetEntriesFast(); j++) { auto *m = (AtMCTrack *)mc->At(j); if (m && m->GetPdgCode() == pdg) has = true; }
         if (!has) continue;
      }
      for (auto &tr : c) if (tr.GetHitArray().size() >= 30) rows.push_back({ev, tr.GetTrackID()});
   }
   if (fs) fs->Close();
}
}
void order_score_pra_dp(TString what, TString O = "/mnt/f/a1975/order_test/", TString Q = "/mnt/f/a1975/pra_test/", TString pre = "f_")
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   const bool isData = what == "data";
   std::vector<TString> tags, recos, sims;
   int pdg = 0; double mass = 1875.613;
   if (isData) { tags.push_back("data0"); tags.push_back("data1"); }
   else if (what == "d") { pdg = 1000010020; for (int s = 5001; s <= 5004; ++s) { tags.push_back(Form("d%d", s)); recos.push_back(Form("/mnt/f/a1975_C16_dd_sim/el_s%d_reco.root", s)); sims.push_back(Form("/mnt/f/a1975_C16_dd_sim/el_s%d_sim.root", s)); } }
   else { pdg = 1000010030; mass = 2808.921; for (const char *s : {"gs", "ex2"}) { tags.push_back(Form("t%s", s)); recos.push_back(Form("/mnt/f/a1975_C16_dt_sim/%s_s3001_reco.root", s)); sims.push_back(Form("/mnt/f/a1975_C16_dt_sim/%s_s3001_sim.root", s)); } }
   TFile *fp = isData ? TFile::Open("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root") : nullptr;
   TTree *w = isData ? (TTree *)fp->Get("w") : nullptr;
   const int NS = 8;
   const double bl[NS + 1] = {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.80};
   double dm[NS] = {0}, md = 1, bA = 0, s1g = 0, sp3 = 0;
   int v3 = 0, ok1 = 0;
   if (isData) {
      TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0";
      w->Draw("s1_g", sel, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      md = x[x.size() / 2];
      for (int k = 0; k < NS; ++k) {
         TH1D h("h", "", 600, -4, 2);
         w->Draw(Form("log(s1_g/%g*bA)>>h", md), sel + Form("&&bA>%g&&bA<%g", bl[k], bl[k + 1]), "goff");
         h.Smooth(2);
         dm[k] = h.GetBinCenter(h.GetMaximumBin());
      }
      w->SetBranchAddress("bA", &bA); w->SetBranchAddress("s1_g", &s1g); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("ok1", &ok1);
   }
   const int NM = 3;
   const char *nm[NM] = {"production", "offline guard (mode 9)", "AtPRA::SetRemoveStrayHits"};
   long nAll[NS] = {0}, nFit[NS][NM] = {{0}}, nGood[NS][NM] = {{0}}, lost = 0;
   Long64_t rowBase = isData ? 36860 : 0;
   std::vector<std::pair<Long64_t, int>> rowsData;
   if (isData) ordpra::rowMap("/mnt/f/a1975/reco_d2_dv1104/run_0026_multifit_reco.root", "", 0, 8000, rowsData);
   Long64_t iData = 0;
   for (size_t it = 0; it < tags.size(); ++it) {
      // the new fits, by (event, trackID)
      std::map<std::pair<Long64_t, int>, ordpra::R> res;
      {
         TFile fn(Q + pre + tags[it] + ".root");
         auto *n = (TTree *)fn.Get("o");
         Long64_t e; int t2, ok; double ke, th;
         n->SetBranchAddress("ev", &e); n->SetBranchAddress("tid", &t2); n->SetBranchAddress("ok0", &ok); n->SetBranchAddress("ke0", &ke); n->SetBranchAddress("th0", &th);
         for (Long64_t i = 0; i < n->GetEntries(); ++i) { n->GetEntry(i); res[{e, t2}] = {ok, ke, th}; }
      }
      std::vector<std::pair<Long64_t, int>> rows;
      if (!isData) ordpra::rowMap(recos[it], sims[it], pdg, 1000000000, rows);
      TFile fa(O + tags[it] + ".root"), fg(O + "g_" + tags[it] + ".root");
      auto *a = (TTree *)fa.Get("o"); auto *g = (TTree *)fg.Get("o");
      int ok0, ok9, skip; double ke0, th0, ke9, th9, tbr;
      a->SetBranchAddress("ok0", &ok0); a->SetBranchAddress("ke0", &ke0); a->SetBranchAddress("th0", &th0); a->SetBranchAddress("tbr", &tbr); a->SetBranchAddress("skip", &skip);
      g->SetBranchAddress("ok9", &ok9); g->SetBranchAddress("ke9", &ke9); g->SetBranchAddress("th9", &th9);
      for (Long64_t i = 0; i < a->GetEntries(); ++i, ++iData) {
         a->GetEntry(i); g->GetEntry(i);
         int k = -1;
         double ref = -1;
         if (isData) {
            w->GetEntry(iData + rowBase);
            if (!v3 || !ok1 || sp3 >= 6 || s1g <= 0 || bA <= 0) continue;
            for (int j = 0; j < NS; ++j) if (bA > bl[j] && bA < bl[j + 1]) k = j;
            if (k < 0 || std::abs(std::log(s1g / md * bA) - dm[k]) > 0.12) continue;
            ref = bA;
         } else {
            if (skip) continue;
            for (int j = 0; j < NS; ++j) if (tbr > bl[j] && tbr < bl[j + 1]) k = j;
            if (k < 0) continue;
            ref = tbr;
         }
         auto key = isData ? rowsData[iData] : rows[i];
         ordpra::R n;
         auto itn = res.find(key);
         if (itn != res.end()) n = itn->second; else ++lost;
         ++nAll[k];
         const int oks[NM] = {ok0, ok9, n.ok};
         const double kes[NM] = {ke0, ke9, n.ke}, ths[NM] = {th0, th9, n.th};
         for (int m = 0; m < NM; ++m) {
            if (!oks[m]) continue;
            ++nFit[k][m];
            const double b = std::sqrt(kes[m] * kes[m] + 2 * kes[m] * mass) / 299.792458;
            const double lb = std::log(b / ref);
            nGood[k][m] += isData ? (std::abs(ordpra::exEl(kes[m], ths[m], 180.0)) < 2 && std::abs(lb) < 0.5) : (std::abs(lb) < 0.2);
         }
      }
   }
   printf("%s: tracks no longer reaching 30 hits or otherwise absent from the new reco: %ld\n", what.Data(), lost);
   long tot = 0, tF[NM] = {0}, tG[NM] = {0};
   for (int k = 0; k < NS; ++k) {
      if (nAll[k] < 15) continue;
      tot += nAll[k];
      printf(" Brho %.2f-%.2f  N %ld:", bl[k], bl[k + 1], nAll[k]);
      for (int m = 0; m < NM; ++m) { tF[m] += nFit[k][m]; tG[m] += nGood[k][m]; printf("  GOOD %5.1f %%", 100. * nGood[k][m] / nAll[k]); }
      printf("\n");
   }
   printf(" ALL %ld\n", tot);
   for (int m = 0; m < NM; ++m) printf("   %-30s fit %5.1f %%  GOOD %5.1f %%\n", nm[m], 100. * tF[m] / tot, 100. * tG[m] / tot);
}
