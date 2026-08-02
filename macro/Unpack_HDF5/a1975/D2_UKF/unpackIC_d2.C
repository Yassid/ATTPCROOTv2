/// @file unpackIC_d2.C
/// @brief Extract ONLY the ion-chamber amplitude per event from the FRIB DAQ stream of an
///        a1975 D2-target run, into a small ntuple -- no 1.4 GB AtRawEvent file per run.
///
/// The D2 (d,p)/(d,t) reconstruction chain has never had a beam gate: reco_d2/ holds no
/// _FRIB.root and every analysis so far used "no IC gate (unavailable for D2)". That is
/// wrong -- the frib/evt group IS present in the D2 h5 files. On run_0016 the IC spectrum
/// shows TWO beam components, ~1150 (the 16C the old C16_dt_anaFit.C gated with 900-1300)
/// and ~2060, in roughly equal numbers. Half the reactions in the D2 data are therefore
/// induced by the wrong beam particle, and they cannot close 16C two-body kinematics.
///
/// The FRIB stream is event-number ordered and lines up 1:1 with the GET stream: on
/// run_0016 the FRIB and reco timestamps differ by a constant offset (17485284 ticks) at
/// every entry, so entry i here is entry i of <run>_reco.root / <run>_genfitter_*.root.
/// The timestamp is stored anyway so the alignment can be re-checked per run (ic_align.C).
/// NOTE the clock WRAPS AT 2^32 partway through a long run -- around entry 40000 of
/// run_0016 the raw difference jumps to 4312452581 = 17485284 + 2^32. Compare the offset
/// modulo 2^32 or a perfectly aligned run looks broken.
///
///   root -b -q 'unpackIC_d2.C+("run_0016", -1, "/mnt/f/a1975/ic_d2/")'
///
/// Output <run>_IC.root, TTree "ic" with: evt, icmax, icsum, ts, mult (~1 MB/run).

/// Run INTERPRETED (no ACLiC "+"): AtHDFUnpacker.h pulls in H5Cpp.h, whose include path is
/// not on the ACLiC search list, so compiling this macro fails. The dictionary from
/// libAtReconstruction is all cling needs.

void unpackIC_d2(TString fileName = "run_0016", Long64_t nEvents = -1, TString outDir = "/mnt/f/a1975/ic_d2/",
                 Int_t icTbLo = 1000, Int_t icTbHi = 1350, Int_t traceIdx = 0)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = "/mnt/f/a1975/h5/" + fileName + ".h5";
   TString mapDir = dir + "/scripts/ANL2023.xml";
   TString digiParFile = dir + "/parameters/ATTPC.a1975_deuterium.par";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   if (gSystem->AccessPathName(inputFile)) {
      printf("\033[1;31mERROR: %s not found\033[0m\n", inputFile.Data());
      return;
   }
   gSystem->mkdir(outDir, kTRUE);

   // The unpacker is driven directly (no FairRunAna output file) so nothing large is written.
   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(outDir + fileName + "_unpackIC_dummy.root");
   run->SetGeomFile(dir + "/geometry/ATTPC_H1bar_geomanager.root");
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);
   rtdb->getContainer("AtDigiPar");

   auto map = std::make_shared<AtTpcMap>();
   map->ParseXMLMap(mapDir.Data());
   map->GeneratePadPlane();

   auto unpacker = std::make_unique<AtFRIBHDFUnpacker>(map);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(1);
   unpacker->SetBaseLineSubtraction(true);
   unpacker->Init();

   Long64_t N = unpacker->GetNumEvents();
   if (nEvents > 0 && nEvents < N)
      N = nEvents;
   printf("\033[1;32m%s: FRIB stream %lld events, extracting IC for %lld\033[0m\n", fileName.Data(),
          unpacker->GetNumEvents(), N);

   TFile out(outDir + fileName + "_IC.root", "RECREATE");
   int evt = 0, mult = 0;
   float icmax = 0, icsum = 0;
   ULong64_t ts = 0;
   TTree *tr = new TTree("ic", "ion chamber per event");
   tr->Branch("evt", &evt);
   tr->Branch("icmax", &icmax);
   tr->Branch("icsum", &icsum);
   tr->Branch("ts", &ts);
   tr->Branch("mult", &mult);

   // Mirror AtUnpackTask::Exec: the event lives in a TClonesArray slot and the stream is
   // polled with IsLastEvent() -- calling FillRawEvent past the end is what made a plain
   // stack AtRawEvent + fixed trip count die mid-run.
   TClonesArray arr("AtRawEvent", 1);
   for (Long64_t i = 0; i < N; ++i) {
      arr.Clear("C");
      auto *rawp = dynamic_cast<AtRawEvent *>(arr.ConstructedAt(0));
      unpacker->FillRawEvent(*rawp);
      AtRawEvent &raw = *rawp;
      evt = raw.GetEventID();
      ts = raw.GetTimestamp(0);
      auto &gt = raw.GetGenTraces();
      mult = gt.size();
      icmax = -1;
      icsum = 0;
      if ((int)gt.size() > traceIdx && gt[traceIdx]) {
         auto &adc = gt[traceIdx]->GetADC();
         for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b) {
            icmax = std::max(icmax, (float)adc[b]);
            icsum += adc[b];
         }
      }
      tr->Fill();
      if (i % 5000 == 0)
         printf("  %lld / %lld\n", i, N);
      if (unpacker->IsLastEvent()) {
         printf("  reached last event at %lld\n", i);
         break;
      }
   }
   out.cd();
   tr->Write();
   out.Close();
   gSystem->Unlink(outDir + fileName + "_unpackIC_dummy.root");
   timer.Stop();
   printf("\033[1;32m%s: wrote %lld IC entries -> %s%s_IC.root (%.0f s)\033[0m\n", fileName.Data(), N, outDir.Data(),
          fileName.Data(), timer.RealTime());
}
