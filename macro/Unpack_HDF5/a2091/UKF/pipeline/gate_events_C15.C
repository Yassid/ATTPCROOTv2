/// @file gate_events_C15.C
/// @brief Build a small FIT INPUT holding only IC-gated (15C) events AND, within them,
///        only proton-gated tracks. Reads <run>_slim.root (AtPatternEvent) + <run>_FRIB.root
///        (IC by entry idx). Keeps events with IC max in [icLo,icHi]; inside each, keeps only
///        tracks whose Spyral PID falls in the proton polygon. Writes <outDir>/<run>_reco.root
///        WITH FairRoot metadata (BranchList/FileHeader/...) so the fitters read it.
///
///   root -b -q 'gate_events_C15.C("run_0138", ...)'
static TCutG *LoadCut(const char *path, const char *name)
{
   std::ifstream in(path); if (!in) return nullptr;
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto pos = s.find('[', s.find("vertices")); if (pos == std::string::npos) return nullptr;
   std::vector<double> nums; const char *p = s.c_str() + pos, *end = s.c_str() + s.size(); int depth = 0;
   while (p < end) { if (*p=='[') depth++; if (*p==']'){depth--; if(depth<=0){++p;break;}}
      char *np=nullptr; double v=strtod(p,&np); if(np!=p){nums.push_back(v);p=np;} else ++p; }
   auto *cut = new TCutG(name, nums.size()/2);
   for (size_t i=0,k=0;i+1<nums.size();i+=2,++k) cut->SetPoint(k,nums[i],nums[i+1]);
   return cut;
}

// count pulses (peaks above thr, hysteresis at thr/2) in adc[tbLo,tbHi]
static int CountPulses(const std::vector<Double_t> &adc, double thr, int tbLo, int tbHi)
{
   int n = 0, b = tbLo, N = (int)adc.size(); if (tbHi > N) tbHi = N;
   while (b < tbHi) { if (adc[b] > thr) { ++n; while (b < tbHi && adc[b] > thr * 0.5) ++b; } else ++b; }
   return n;
}

/// icLo/icHi default to the 15C beam peak chosen in pid/ic_15C.json (the dominant 1142 ADC peak,
/// edges in the valleys), NOT the old 950-1350 placeholder which cut through the 1372 contaminant
/// peak. icTbLo/icTbHi/peakThr/pkTbLo/pkTbHi are kept only for reference: the IC observable is now
/// read pre-computed from <run>_ic.root, which was produced with exactly these values.
void gate_events_C15(TString run, TString inDir = "/home/yassid/a2091_C15_reco/",
                      TString outDir = "/home/yassid/a2091_C15_fit/in/", TString refReco = "", double icLo = 979.5,
                      double icHi = 1278.8, Int_t icTbLo = 1050, Int_t icTbHi = 1250, double bField = 2.85,
                      TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_HDF5/a2091/"
                                           "UKF/pid/proton_15C.json",
                      double thMin = 90.0, double peakThr = 200, Int_t pkTbLo = 800, Int_t pkTbHi = 1500,
                      TString icDir = "/home/yassid/a2091_C15_ic/")
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   // IC now comes from the compact <run>_ic.root summary (entry, icmax, npulse) instead of the
   // full <run>_FRIB.root: the raw FRIB traces are ~36 kB/event (~54 GB over the run set) and only
   // these two numbers were ever used, so they are no longer persisted. See pipeline/icsum_C15.C.
   TString ff = icDir + run + "_ic.root";
   // The slim AtPatternEvent cache was never built for a2091; fall back to the full reco, which
   // carries the same AtPatternEvent branch.
   TString rf = inDir + run + "_slim.root";
   if (gSystem->AccessPathName(rf)) rf = inDir + run + "_reco.root";
   if (gSystem->AccessPathName(ff)) { printf("SKIP %s (no IC summary %s)\n", run.Data(), ff.Data()); return; }
   if (gSystem->AccessPathName(rf)) { printf("SKIP %s (no reco %s)\n", run.Data(), rf.Data()); return; }
   gSystem->mkdir(outDir.Data(), kTRUE);
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");
   if (!pcut) { printf("ERR: no proton gate\n"); return; }

   // IC max per entry index -> keep flags. `entry` in the summary IS the FRIB tree entry index,
   // i.e. the same index the reco tree is walked with below.
   TFile *fF = TFile::Open(ff);
   TTree *tF = fF ? (TTree *)fF->Get("ic") : nullptr;
   if (!tF || tF->GetEntries() == 0) { printf("SKIP %s (empty IC summary)\n", run.Data()); if (fF) fF->Close(); return; }
   Int_t icEntry, icNpulse; Float_t icMax;
   tF->SetBranchAddress("entry", &icEntry);
   tF->SetBranchAddress("icmax", &icMax);
   tF->SetBranchAddress("npulse", &icNpulse);
   std::vector<char> icKeep(tF->GetEntries(), 0);
   Long64_t nIC = 0;
   for (Long64_t i = 0; i < tF->GetEntries(); i++) {
      tF->GetEntry(i);
      // IC gate: amplitude in [icLo,icHi] AND exactly ONE pulse (reject pile-up)
      if (icEntry >= 0 && icEntry < (Int_t)icKeep.size() && icMax >= icLo && icMax <= icHi && icNpulse == 1) {
         icKeep[icEntry] = 1;
         ++nIC;
      }
   }
   fF->Close();
   printf("%s : IC gate [%.0f,%.0f] 1-pulse -> %lld / %lld events\n", run.Data(), icLo, icHi, nIC,
          (Long64_t)icKeep.size());

   TFile *fin = TFile::Open(rf);
   TTree *t = (TTree *)fin->Get("cbmsim");
   TClonesArray *pe = nullptr; t->SetBranchAddress("AtPatternEvent", &pe);
   TString of = outDir + run + "_reco.root";
   TFile *fout = new TFile(of, "RECREATE", "", 1);
   TTree *nt = t->CloneTree(0);
   Long64_t nEvt = 0, nTrk = 0;
   for (Long64_t i = 0; i < t->GetEntries(); i++) {
      if (i >= (Long64_t)icKeep.size() || !icKeep[i]) continue; // IC gate on event
      t->GetEntry(i);
      if (pe->GetEntries() == 0) continue;
      auto *p = (AtPatternEvent *)pe->At(0); if (!p) continue;
      std::vector<AtTrack> protons; // proton gate + backward (theta_lab>thMin) on tracks
      for (auto &trk : p->GetTrackCand()) {
         AtTrack &tr = const_cast<AtTrack &>(trk);
         auto r = spy.Estimate(tr);
         if (r.valid && pcut->IsInside(r.sqrtdEdx, r.brho) && r.polar * TMath::RadToDeg() > thMin)
            protons.push_back(tr);
      }
      if (protons.empty()) continue;
      nTrk += protons.size();
      p->SetTrackCand(std::move(protons));
      nt->Fill();
      ++nEvt;
   }
   nt->Write();

   TList bl; for (auto b : *nt->GetListOfBranches()) bl.Add(new TObjString(b->GetName()));
   fout->cd(); bl.Write("BranchList", TObject::kSingleKey);
   TList etb; etb.Write("TimeBasedBranchList", TObject::kSingleKey);
   if (refReco.Length()) { TFile *fr = TFile::Open(refReco, "READ"); if (fr && !fr->IsZombie()) {
      if (auto *fh = fr->Get("FileHeader")) { fout->cd(); fh->Write("FileHeader"); }
      if (auto *cb = fr->Get("cbmout")) { fout->cd(); cb->Write("cbmout", TObject::kSingleKey); } fr->Close(); } }
   fout->Write("", TObject::kOverwrite);
   printf("%s : %lld gated-proton events, %lld tracks -> %s\n", run.Data(), nEvt, nTrk, of.Data());
   fout->Close(); fin->Close();
}
