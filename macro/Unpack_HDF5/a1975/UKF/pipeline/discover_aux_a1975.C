/// @file discover_aux_a1975.C
/// @brief Discover the IC / auxiliary channels in the a1975 HDF5 by registering
/// every electronics address NOT in the pad map as an aux pad, then ranking by
/// occupancy. The ion chamber fires on essentially every beam event, so it shows
/// up as a high-occupancy unmapped channel (amplitude ~950-1350 per the Spyral
/// C16 config gate). Reveals the (cobo,asad,aget,channel) to add via AddAuxPad.
///
/// Run: root -b -q 'discover_aux_a1975.C("run_0116", 500)'

void discover_aux_a1975(TString fileName = "run_0116", Long64_t nEvents = 500, Int_t maxCobo = 11)
{
   gSystem->Load("libAtReconstruction.so");
   TString dir = getenv("VMCWORKDIR");
   TString inputFile = "/mnt/f/a1975/h5/" + fileName + ".h5";
   TString mapXml = dir + "/scripts/ANL2023.xml";
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());

   // --- Parse the pad map XML to collect mapped (cobo,asad,aget,channel) keys ---
   auto key = [](int co, int as, int ag, int ch) { return ((co * 100 + as) * 100 + ag) * 100 + ch; };
   std::set<long> padSet;
   {
      std::ifstream in(mapXml.Data());
      std::string line;
      int co = -1, as = -1, ag = -1, ch = -1;
      auto grab = [](const std::string &l, const char *tag, int &out) {
         auto a = l.find(std::string("<") + tag + ">");
         if (a == std::string::npos)
            return;
         a = l.find('>', a) + 1;
         auto b = l.find('<', a);
         out = std::atoi(l.substr(a, b - a).c_str());
      };
      while (std::getline(in, line)) {
         grab(line, "CoboID", co);
         grab(line, "AsadID", as);
         grab(line, "AgetID", ag);
         grab(line, "ChannelID", ch);
         if (line.find("</ANL2023>") != std::string::npos && co >= 0) {
            padSet.insert(key(co, as, ag, ch));
            co = as = ag = ch = -1;
         }
      }
      printf("Parsed %zu mapped pad channels\n", padSet.size());
   }

   FairRunAna *run = new FairRunAna();
   run->SetOutputFile(fileName + "_auxscan.root");
   run->SetGeomFile(dir + "/geometry/ATTPC_H1bar_geomanager.root");
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open((dir + "/parameters/ATTPC.a1954.par").Data(), "in");
   rtdb->setFirstInput(parIo1);
   rtdb->getContainer("AtDigiPar");

   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   fAtMapPtr->ParseXMLMap(mapXml.Data());
   fAtMapPtr->GeneratePadPlane();

   // Register every UNMAPPED address as aux.
   int nAux = 0;
   for (int co = 0; co <= maxCobo; ++co)
      for (int as = 0; as < 4; ++as)
         for (int ag = 0; ag < 4; ++ag)
            for (int ch = 0; ch < 68; ++ch)
               if (padSet.find(key(co, as, ag, ch)) == padSet.end()) {
                  fAtMapPtr->AddAuxPad({co, as, ag, ch}, Form("c%d_a%d_g%d_ch%02d", co, as, ag, ch));
                  ++nAux;
               }
   printf("Registered %d unmapped addresses as aux\n", nAux);

   auto unpacker = std::make_unique<AtHDFUnpacker>(fAtMapPtr);
   unpacker->SetInputFileName(inputFile.Data());
   unpacker->SetNumberTimestamps(2);
   unpacker->SetBaseLineSubtraction(true);
   auto unpackTask = new AtUnpackTask(std::move(unpacker));
   unpackTask->SetPersistence(true);
   run->AddTask(unpackTask);
   run->Init();
   auto N = unpackTask->GetNumEvents();
   if (nEvents > 0 && nEvents < N)
      N = nEvents;
   run->Run(0, N);

   delete run;
   TFile *f = TFile::Open(fileName + "_auxscan.root");
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *ev = nullptr;
   t->SetBranchAddress("AtRawEvent", &ev);
   std::map<std::string, std::pair<int, double>> seen;
   Long64_t nev = t->GetEntries();
   for (Long64_t i = 0; i < nev; ++i) {
      t->GetEntry(i);
      if (ev->GetEntries() == 0)
         continue;
      auto *r = (AtRawEvent *)ev->At(0);
      if (!r)
         continue;
      for (auto &ap : r->GetAuxPads()) {
         const auto &adc = ap.second.GetADC();
         double mx = 0;
         for (auto v : adc)
            mx = std::max(mx, (double)v);
         if (mx > 80) {
            auto &e = seen[ap.first];
            e.first++;
            e.second = std::max(e.second, mx);
         }
      }
   }
   std::vector<std::pair<std::string, std::pair<int, double>>> v(seen.begin(), seen.end());
   std::sort(v.begin(), v.end(), [](auto &a, auto &b) { return a.second.first > b.second.first; });
   printf("\n\033[1;32m=== Unmapped channels with signal, by occupancy (%lld events) ===\033[0m\n", nev);
   for (int i = 0; i < (int)v.size() && i < 20; ++i)
      printf("  %-16s  fired %5d/%lld (%3.0f%%)  maxADC %.0f\n", v[i].first.c_str(), v[i].second.first, nev,
             100.0 * v[i].second.first / nev, v[i].second.second);
}
