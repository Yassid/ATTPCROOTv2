/// @file run_digi_C14.C
/// @brief Digitize C14_pp_sim.C output with the SAME parameters used to reconstruct the
/// a1954 14C data, so simulation and experiment can be compared one-to-one.
///
/// Uses parameters/ATTPC.a1954_C14.par directly (identical key set to the e20009 sim par
/// the sibling macros use, verified by diff) -- that is the whole point: B = 2.85 T,
/// H2 600 torr, SamplingRate 6.0, DriftVelocity 1.30, TBEntrance 480.
///
/// The drift velocity is the parameter currently under question (the dv scan on data
/// prefers 1.37-1.49, not the 1.30 in the par). Pass a different par here to digitize the
/// same truth at a different dv -- that is the clean way to see what dv does to the
/// reconstructed Ex, with the truth held fixed.
///
///   root -l 'run_digi_C14.C()'
///   root -l 'run_digi_C14.C("ATTPC.a1954_C14_dv145.par", "./data/output_digi_dv145.root")'

bool reduceFunc(AtRawEvent *evt);

void run_digi_C14(TString paramFile = "ATTPC.a1954_C14.par", TString outputFile = "./data/output_digi.root",
                  Int_t nEvents = 0)
{
   TString inOutDir = "./data/";
   TString scriptfile = "Lookup20150611.xml";
   TString dir = getenv("VMCWORKDIR");
   TString mcFile = inOutDir + "attpcsim.root";

   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapParFile.Data());
   mapping->GeneratePadPlane();

   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kTRUE);

   auto psa = std::make_unique<AtPSAMax>();
   psa->SetThreshold(0);

   AtPSAtask *psaTask = new AtPSAtask(std::move(psa));
   psaTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);

   fRun->Init();

   timer.Start();
   if (nEvents > 0)
      fRun->Run(0, nEvents);
   else
      fRun->Run();
   timer.Stop();

   std::cout << std::endl;
   std::cout << "Digitization finished. par = " << paramFile << std::endl;
   std::cout << "Output file is " << outputFile << std::endl;
   std::cout << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl << std::endl;
}

bool reduceFunc(AtRawEvent *evt)
{
   return (evt->GetNumPads() > 0) && evt->IsGood();
}
