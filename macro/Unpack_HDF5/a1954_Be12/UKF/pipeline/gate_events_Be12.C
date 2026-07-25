/// @file gate_events_Be12.C
/// @brief Build a small FIT INPUT holding only IC-gated (12Be) events AND, within them,
///        only proton-gated tracks. Reads <run>_slim.root (AtPatternEvent) + <run>_FRIB.root
///        (IC by entry idx). Keeps events with IC max in [icLo,icHi]; inside each, keeps only
///        tracks whose Spyral PID falls in the proton polygon. Writes <outDir>/<run>_reco.root
///        WITH FairRoot metadata (BranchList/FileHeader/...) so the fitters read it.
///
///   root -b -q 'gate_events_Be12.C("run_0142", ...)'
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

void gate_events_Be12(TString run, TString inDir = "/home/yassid/a1954_Be12_reco_hdb_slim/",
                      TString outDir = "/home/yassid/a1954_Be12_fit/in/", TString refReco = "", double icLo = 500,
                      double icHi = 900, Int_t icTbLo = 1050, Int_t icTbHi = 1250, double bField = 2.85,
                      TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954_Be12/"
                                           "UKF/pid/proton_12Be.json",
                      double thMin = 90.0, double peakThr = 200, Int_t pkTbLo = 800, Int_t pkTbHi = 1500)
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   TString ff = inDir + run + "_FRIB.root", rf = inDir + run + "_slim.root";
   if (gSystem->AccessPathName(ff) || gSystem->AccessPathName(rf)) { printf("SKIP %s (missing)\n", run.Data()); return; }
   gSystem->mkdir(outDir.Data(), kTRUE);
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");
   if (!pcut) { printf("ERR: no proton gate\n"); return; }

   // IC max per entry index -> keep flags
   TFile *fF = TFile::Open(ff);
   TTree *tF = fF ? (TTree *)fF->Get("cbmsim") : nullptr;
   if (!tF || tF->GetEntries() == 0) { printf("SKIP %s (no FRIB)\n", run.Data()); if (fF) fF->Close(); return; }
   TClonesArray *ra = nullptr; tF->SetBranchAddress("AtRawEvent", &ra);
   std::vector<char> icKeep(tF->GetEntries(), 0);
   for (Long64_t i = 0; i < tF->GetEntries(); i++) {
      tF->GetEntry(i); if (ra->GetEntries() == 0) continue;
      auto *r = (AtRawEvent *)ra->At(0); if (!r || r->GetGenTraces().empty()) continue;
      auto &adc = r->GetGenTraces()[0]->GetADC(); double mx = -1e9;
      for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); b++) mx = std::max(mx, (double)adc[b]);
      // IC gate: amplitude in [icLo,icHi] AND exactly ONE pulse (reject pile-up)
      if (i < (Long64_t)icKeep.size() && mx >= icLo && mx <= icHi && CountPulses(adc, peakThr, pkTbLo, pkTbHi) == 1)
         icKeep[i] = 1;
   }
   fF->Close();

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
