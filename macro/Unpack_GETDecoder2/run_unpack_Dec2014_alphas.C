// Unpacker for the Dec 2014 alpha data (Cris_OneD / Dec2014_alphas).
//
// Adapted from run_unpack_alpha.C, which predates the AT -> At class rename and
// the AtPSAtask refactor (PSA options now live on the AtPSA object, not the task).
//
//   root -l 'run_unpack_Dec2014_alphas.C("runfiles/NSCL/Dec2014_alphas/alpha_run_0080.txt", "alpha_run_0080.root", 300, "Lookup20141208.xml")'

// NB: the libraries are loaded by ./rootlogon.C, which ROOT runs before this macro is
// parsed. Do not add #includes or gSystem->Load calls here -- see rootlogon.C for why.

#define cRED "\033[1;31m"
#define cNORMAL "\033[0m"

bool check_file(const std::string &name)
{
   std::ifstream f(name.c_str());
   return f.good();
}

void run_unpack_Dec2014_alphas(TString dataFile = "runfiles/NSCL/Dec2014_alphas/alpha_run_0080.txt",
                               TString outputFile = "alpha_run_0080.root", Int_t nEvents = 300,
                               TString scriptfile = "Lookup20141208.xml",
                               TString parameterFile = "ATTPC.alpha.par")
{
   if (!check_file(dataFile.Data())) {
      std::cout << cRED << " Run file " << dataFile.Data() << " not found! Terminating..." << cNORMAL << std::endl;
      return;
   }

   TStopwatch timer;
   timer.Start();

   TString dir = getenv("VMCWORKDIR");
   TString scriptdir = dir + "/scripts/" + scriptfile;
   TString dataDir = dir + "/macro/data/";
   TString geomDir = dir + "/geometry/";
   gSystem->Setenv("GEOMPATH", geomDir.Data());

   TString loggerFile = dataDir + "ATTPCLog.log";
   TString digiParFile = dir + "/parameters/" + parameterFile;
   TString geoManFile = dir + "/geometry/ATTPC_v1.2.root";

   FairLogger *fLogger = FairLogger::GetLogger();
   fLogger->SetLogFileName(loggerFile.Data());
   fLogger->SetLogToScreen(kTRUE);
   fLogger->SetLogToFile(kTRUE);
   fLogger->SetLogVerbosityLevel("LOW");

   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(outputFile);
   run->SetGeomFile(geoManFile);

   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setSecondInput(parIo1);

   // The runfiles are .txt lists of per-CoBo .graw files -> separated data.
   Bool_t fUseSeparatedData = dataFile.EndsWith(".txt");

   AtDecoder2Task *fDecoderTask = new AtDecoder2Task();
   fDecoderTask->SetUseSeparatedData(fUseSeparatedData);
   if (fUseSeparatedData)
      fDecoderTask->SetPseudoTopologyFrame(kTRUE);
   fDecoderTask->SetPersistence(kFALSE);
   fDecoderTask->SetMap(scriptdir.Data());
   fDecoderTask->SetMapOpt(0); // ATTPC : 0 - Prototype : 1
   // SetNumCobo is set below, from the CoBos actually present in the run file.
   fDecoderTask->SetEventID(0);

   if (!fUseSeparatedData) {
      fDecoderTask->AddData(dataFile);
      fDecoderTask->SetNumCobo(1);
   } else {
      // AtCore2 allocates decoders 0 .. fNumCobo-1, so the CoBo *number* is not the
      // decoder index: this data has no CoBo 5, so CoBo 6->5, 7->6, 8->7, 9->8.
      // Map each distinct CoBo number to a sequential index in order of appearance.
      // (run_unpack_alpha.C bumped the index on every filename mismatch instead, which
      // misassigns the second chunk of a CoBo once a gap has been passed, and runs off
      // the end of the decoder array on the CoBo9 chunk.)
      std::vector<Int_t> coboSeen;
      std::vector<TString> paths;
      std::vector<Int_t> indices;

      std::ifstream listFile(dataFile.Data());
      TString line;
      while (line.ReadLine(listFile)) {
         Int_t pos = line.Index("CoBo");
         if (pos == kNPOS)
            continue;
         TString num;
         for (Int_t i = pos + 4; i < line.Length() && line[i] >= '0' && line[i] <= '9'; ++i)
            num += line[i];
         if (num.Length() == 0)
            continue;

         Int_t cobo = num.Atoi();
         Int_t idx = -1;
         for (size_t k = 0; k < coboSeen.size(); ++k)
            if (coboSeen[k] == cobo) {
               idx = (Int_t)k;
               break;
            }
         if (idx < 0) {
            coboSeen.push_back(cobo);
            idx = (Int_t)coboSeen.size() - 1;
         }
         paths.push_back(line);
         indices.push_back(idx);
      }

      fDecoderTask->SetNumCobo((Int_t)coboSeen.size());
      for (size_t i = 0; i < paths.size(); ++i)
         fDecoderTask->AddData(paths[i], indices[i]);

      std::cout << " CoBos found: " << coboSeen.size() << " -> ";
      for (size_t k = 0; k < coboSeen.size(); ++k)
         std::cout << "CoBo" << coboSeen[k] << "=" << k << " ";
      std::cout << std::endl;
   }

   run->AddTask(fDecoderTask);

   // PSA: options moved from the task onto the AtPSA object on this branch.
   // AtPSASimple2 is the ATTPC analyzer (old SetPSAMode(1)); AtPSAProto is the prototype one.
   AtPSASimple2 *psa = new AtPSASimple2();
   psa->SetThreshold(20);
   psa->SetMaxFinder();          // use either the max finder or the peak finder, not both
   psa->SetBaseCorrection(kTRUE);
   psa->SetTimeCorrection(kFALSE);

   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(kTRUE);
   run->AddTask(psaTask);

   run->Init();
   run->Run(0, nEvents);

   timer.Stop();
   std::cout << std::endl;
   std::cout << "Macro finished succesfully." << std::endl;
   std::cout << "- Input      : " << dataFile << std::endl;
   std::cout << "- Lookup     : " << scriptfile << std::endl;
   std::cout << "- Output file: " << outputFile << std::endl;
   std::cout << "- Real time  : " << timer.RealTime() << " s" << std::endl;
   std::cout << "- CPU time   : " << timer.CpuTime() << " s" << std::endl;
}
