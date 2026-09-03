/// @file dump_ic_traces_C15d.C
/// @brief Dump example ion-chamber traces (FRIB generic channel 0) as JSON for the trace viewer.
///
///   root -b -q 'pid/dump_ic_traces_C15d.C("run_0017")'
///
/// WHY. icsum_C15d.C reduces each trace to two numbers -- icmax and npulse -- and every beam
/// selection downstream rests on them. The reduction is not obvious from the numbers alone: the
/// AMPLITUDE is taken over a fixed [1050,1250] while the PULSE COUNT is taken over [800,1500], so
/// a pulse that arrives outside the narrower window is counted but measured as baseline. This dump
/// keeps the raw samples so that reduction can be watched happening on real pulses.
///
/// The sample is deliberately STRATIFIED, not the first N events: the Z ladder spans a factor 7 in
/// amplitude and the interesting failures are ~8 % of events, so a sequential sample would be all
/// carbon and oxygen and would show none of them.

static int CountPulsesJ(const std::vector<Double_t> &adc, double thr, int lo, int hi)
{
   int n = 0, b = lo, N = (int)adc.size();
   if (hi > N) hi = N;
   while (b < hi) {
      if (adc[b] > thr) { ++n; while (b < hi && adc[b] > thr * 0.5) ++b; }
      else ++b;
   }
   return n;
}

void dump_ic_traces_C15d(TString run = "run_0017", TString inDir = "/home/yassid/C15d_ic/tmp/",
                         TString outJson = "", Int_t perClass = 45,
                         // the sample kept per trace: the full 2048 would be mostly flat baseline
                         Int_t tbLo = 700, Int_t tbHi = 1600,
                         // the two windows whose disagreement is the thing being shown
                         Int_t icTbLo = 1050, Int_t icTbHi = 1250, Int_t pkTbLo = 800, Int_t pkTbHi = 1500,
                         Double_t peakThr = 200)
{
   gSystem->Load("libAtReconstruction.so");
   if (outJson.IsNull())
      outJson = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/ic_traces.json";

   TString fin = inDir + run + "_FRIB.root";
   if (gSystem->AccessPathName(fin)) { printf("MISSING %s\n", fin.Data()); return; }
   TFile *f = TFile::Open(fin);
   TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
   if (!t) { printf("no cbmsim in %s\n", fin.Data()); return; }
   TClonesArray *ra = nullptr;
   t->SetBranchAddress("AtRawEvent", &ra);

   // classes: the Z ladder anchored on carbon 1175, plus the two failure modes
   struct Cls { const char *name; double lo, hi; };
   const Cls CLS[] = {
      {"Z=3 Li ~295",   200,  400},
      {"Z=5 B ~775",    650,  900},
      {"Z=6 C ~1175",  1050, 1300},
      {"Z=7 N ~1595",  1480, 1700},
      {"Z=8 O ~2055",  1900, 2200},
   };
   const int NC = sizeof(CLS) / sizeof(CLS[0]);
   std::vector<int> got(NC, 0);
   int gotOutside = 0, gotMulti = 0;

   std::vector<TString> rows;
   const Long64_t N = t->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!ra || ra->GetEntries() == 0) continue;
      auto *r = (AtRawEvent *)ra->At(0);
      if (!r || r->GetGenTraces().empty()) continue;
      const auto &adc = r->GetGenTraces()[0]->GetADC();
      const int n = (int)adc.size();
      if (n < pkTbHi) continue;

      double fixMax = -1e9;
      for (int b = icTbLo; b < icTbHi && b < n; ++b) fixMax = std::max(fixMax, (double)adc[b]);
      double trueMax = -1e9; int trueTb = -1;
      for (int b = pkTbLo; b < pkTbHi && b < n; ++b)
         if (adc[b] > trueMax) { trueMax = adc[b]; trueTb = b; }
      const int npul = CountPulsesJ(adc, peakThr, pkTbLo, pkTbHi);
      const bool outside = (trueTb < icTbLo || trueTb >= icTbHi);

      // decide whether this trace is wanted: fill the classes, and keep the failure modes
      int cls = -1;
      for (int c = 0; c < NC; ++c)
         if (trueMax >= CLS[c].lo && trueMax <= CLS[c].hi) { cls = c; break; }
      bool want = false;
      const char *tag = "";
      if (outside && trueMax > 150 && gotOutside < perClass) { want = true; ++gotOutside; tag = "max OUTSIDE the fixed window"; }
      else if (npul > 1 && gotMulti < perClass)              { want = true; ++gotMulti;   tag = "multi-pulse"; }
      else if (cls >= 0 && got[cls] < perClass)              { want = true; ++got[cls];   tag = CLS[cls].name; }
      if (!want) continue;

      TString adcs = "";
      for (int b = tbLo; b < tbHi && b < n; ++b)
         adcs += TString::Format("%s%d", b > tbLo ? "," : "", (int)std::lround(adc[b]));
      rows.push_back(TString::Format(
         "{\"ev\":%lld,\"cls\":\"%s\",\"fixMax\":%.0f,\"trueMax\":%.0f,\"trueTb\":%d,"
         "\"np\":%d,\"outside\":%d,\"adc\":[%s]}",
         i, tag, fixMax, trueMax, trueTb, npul, outside ? 1 : 0, adcs.Data()));

      bool done = (gotOutside >= perClass && gotMulti >= perClass);
      for (int c = 0; c < NC; ++c) done = done && (got[c] >= perClass);
      if (done) break;
   }

   FILE *o = fopen(outJson.Data(), "w");
   fprintf(o, "{\"run\":\"%s\",\"tbLo\":%d,\"tbHi\":%d,\"icTbLo\":%d,\"icTbHi\":%d,"
              "\"pkTbLo\":%d,\"pkTbHi\":%d,\"thr\":%g,\"traces\":[",
           run.Data(), tbLo, tbHi, icTbLo, icTbHi, pkTbLo, pkTbHi, peakThr);
   for (size_t k = 0; k < rows.size(); ++k) fprintf(o, "%s%s", k ? "," : "", rows[k].Data());
   fprintf(o, "]}");
   fclose(o);

   printf("wrote %s : %zu traces\n", outJson.Data(), rows.size());
   for (int c = 0; c < NC; ++c) printf("  %-14s %3d\n", CLS[c].name, got[c]);
   printf("  %-14s %3d\n  %-14s %3d\n", "outside", gotOutside, "multi-pulse", gotMulti);
}
