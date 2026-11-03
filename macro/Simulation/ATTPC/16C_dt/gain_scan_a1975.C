/// @file gain_scan_a1975.C
/// @brief Chooses the simulation gain at the DATA reco settings (thr 40, HDBSCAN mcs 0, ms 3), per species, by comparing sim and data in bins of Spyral Brho:
///        median hits per track, median charge per hit, median sqrt(dE/dx) of Spyral, and how many truth-matched sim tracks fall inside the data's PID gate.
///        Data: tracks inside the gate, first 8000 events of runs 0016 and 0026. Sim: the track with most hits among those whose Spyral Brho is within 20 % of truth.
///   root -l -b -q 'gain_scan_a1975.C()'
namespace gs {
struct Row { std::vector<double> h, q, d; };
const int NB = 6;
const double bl[NB + 1] = {0.10, 0.20, 0.30, 0.40, 0.55, 0.75, 1.05};
int bin(double b) { for (int k = 0; k < NB; ++k) if (b >= bl[k] && b < bl[k + 1]) return k; return -1; }
double med(std::vector<double> v) { if (v.empty()) return 0; std::sort(v.begin(), v.end()); return v[v.size() / 2]; }
void addTrack(AtTrack &tr, const AtTools::AtSpyralResult &r, Row *row)
{
   const int k = bin(r.brho);
   if (k < 0) return;
   std::vector<double> q;
   for (auto &h : tr.GetHitArray()) q.push_back(h->GetCharge());
   row[k].h.push_back(tr.GetHitArray().size()); row[k].q.push_back(med(q)); row[k].d.push_back(r.sqrtdEdx);
}
void data(const AtTools::AtParticleID &gate, Row *row)
{
   AtTools::AtSpyralPID spy; spy.SetBField(2.85);
   for (const char *run : {"0016", "0026"}) {
      TFile f(Form("/mnt/f/a1975/reco_d2_dv1104/run_%s_multifit_reco.root", run));
      auto *t = (TTree *)f.Get("cbmsim");
      TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
      for (Long64_t e = 0; e < 8000 && e < t->GetEntries(); ++e) {
         if (t->GetEntry(e) <= 0 || !pa->GetEntries()) continue;
         for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
            auto r = spy.Estimate(tr);
            if (r.valid && gate.IsInside(r.sqrtdEdx, r.brho)) addTrack(tr, r, row);
         }
      }
   }
}
void sim(TString reco, TString mcf, int pdg, double massMeV, const AtTools::AtParticleID &gate, Row *row, long *nT, long *nM, long *nG)
{
   AtTools::AtSpyralPID spy; spy.SetBField(2.85);
   TFile f(reco);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("AtPatternEvent*", 1); t->SetBranchAddress("AtPatternEvent", &pa);
   TFile fm(mcf);
   auto *tm = (TTree *)fm.Get("cbmsim");
   TClonesArray *mc = nullptr;
   tm->SetBranchStatus("*", 0); tm->SetBranchStatus("MCTrack*", 1); tm->SetBranchAddress("MCTrack", &mc);
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e); tm->GetEntry(e);
      double tb = -1;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (m && m->GetPdgCode() == pdg) { tb = m->GetP() * 1000. / 299.792458; break; }
      }
      if (tb < 0) continue;
      const int kt = bin(tb);
      if (kt < 0) continue;
      ++nT[kt];
      if (!pa->GetEntries()) continue;
      AtTrack *best = nullptr; AtTools::AtSpyralResult rb; size_t nb = 0;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto r = spy.Estimate(tr);
         if (r.valid && std::abs(std::log(r.brho / tb)) < 0.2 && tr.GetHitArray().size() > nb) { best = &tr; rb = r; nb = tr.GetHitArray().size(); }
      }
      if (!best) continue;
      ++nM[kt];
      nG[kt] += gate.IsInside(rb.sqrtdEdx, rb.brho);
      addTrack(*best, rb, row);
   }
}
void report(const char *sp, TString mcf, int pdg, double mass, const char *gateFile, std::vector<int> gains, TString O)
{
   auto gate = AtTools::AtParticleID::LoadJSON(gateFile);
   Row D[NB];
   data(gate, D);
   printf("\n=== %s, gate %s\n", sp, gateFile);
   printf("Brho bin   : "); for (int k = 0; k < NB; ++k) printf("  %.2f-%.2f", bl[k], bl[k + 1]); printf("\n");
   printf("data N     : "); for (int k = 0; k < NB; ++k) printf("  %9zu", D[k].h.size()); printf("\n");
   printf("data hits  : "); for (int k = 0; k < NB; ++k) printf("  %9.0f", med(D[k].h)); printf("\n");
   printf("data q/hit : "); for (int k = 0; k < NB; ++k) printf("  %9.0f", med(D[k].q)); printf("\n");
   printf("data sqdEdx: "); for (int k = 0; k < NB; ++k) printf("  %9.2f", med(D[k].d)); printf("\n");
   for (int g : gains) {
      Row S[NB]; long nT[NB] = {0}, nM[NB] = {0}, nG[NB] = {0};
      sim(O + Form("%s_g%d.root", sp, g), mcf, pdg, mass, gate, S, nT, nM, nG);
      printf("g %6d N : ", g); for (int k = 0; k < NB; ++k) printf("  %4ld/%4ld", nM[k], nT[k]); printf("\n");
      printf("  hits s/d : "); for (int k = 0; k < NB; ++k) printf("  %9.2f", S[k].h.size() >= 10 && med(D[k].h) > 0 ? med(S[k].h) / med(D[k].h) : -1.); printf("\n");
      printf("  q/h  s/d : "); for (int k = 0; k < NB; ++k) printf("  %9.2f", S[k].q.size() >= 10 && med(D[k].q) > 0 ? med(S[k].q) / med(D[k].q) : -1.); printf("\n");
      printf("  sqdE s/d : "); for (int k = 0; k < NB; ++k) printf("  %9.2f", S[k].d.size() >= 10 && med(D[k].d) > 0 ? med(S[k].d) / med(D[k].d) : -1.); printf("\n");
      printf("  in gate %%: "); for (int k = 0; k < NB; ++k) printf("  %9.1f", nM[k] >= 10 ? 100. * nG[k] / nM[k] : -1.); printf("\n");
   }
}
}
void gain_scan_a1975(TString O = "/mnt/f/a1975/sim_rebuild/gainscan/")
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   TString P = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1975/D2_UKF/pid/";
   std::vector<int> G = {35000, 50000, 75000, 100000, 150000};
   gs::report("t", "/mnt/f/a1975_C16_dt_sim/gs_s3001_sim.root", 1000010030, 2808.921, (P + "triton_d2_dv1104.json").Data(), G, O);
   gs::report("p", "/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root", 2212, 938.272, (P + "proton_d2_dv1104.json").Data(), G, O);
   gs::report("d", "/mnt/f/a1975_C16_dd_sim/el_s5001_sim.root", 1000010020, 1875.613, (P + "deuteron_band_d2.json").Data(), G, O);
}
