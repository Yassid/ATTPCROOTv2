/// @file pidPass_Ar46.C
/// @brief Compute and persist the Spyral PID observables (brho, sqrt(dE/dx) per track) from
/// <tag>_reco.root into <tag>_pid.root. NO GATE, NO SELECTION, NO FIT.
///
/// This is the last stage the automated chain runs. It exists so the PID plane can be drawn and
/// the deuteron gate placed BY HAND before anything is fitted -- AtPIDTask only computes the two
/// observables per pattern track and writes them out; it decides nothing.
///
/// The a1975 version of this (pipeline/pidPass_a1975.C) hardcodes ATTPC.a1954.par. This one takes
/// the 46Ar par so the drift velocity and gas description match the digitisation that produced
/// the input -- AtPIDTask reads the par through the runtime db like every other task.
///
/// AtSpyralPID works off the pattern track and the hits' charge, so no fit is needed to get here.
/// On a1975 the fit was run purely because AtGenfitter was what carried AtPIDTask along, at about
/// an hour per sample against a minute for this.
///
/// B FIELD: 2.85 T, AtPIDTask's own default and the value in the par. The rigidity scale of the
/// plane follows from it, so if the field changes, the gate must be redrawn.
///
///   root -b -q 'pidPass_Ar46.C("gs_s3001","/mnt/f/ar46_3hed_OLD_2.85T_placeholder/","/mnt/f/ar46_3hed_OLD_2.85T_placeholder/")'

void pidPass_Ar46(TString fileName = "gs_s3001", TString inDir = "./data/", TString outDir = "./data/",
                  Double_t bField = 2.85, TString paramFile = "ATTPC.46Ar_3Hed_sim.par")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("WARNING");

   TString dir = getenv("VMCWORKDIR");
   TString inputFile = inDir + fileName + "_reco.root";
   TString outputFile = outDir + fileName + "_pid.root";
   TString digiParFile = dir + "/parameters/" + paramFile;
   if (gSystem->AccessPathName(inputFile.Data())) {
      std::cout << "\033[1;31mERROR: " << inputFile << " not found.\033[0m\n";
      return;
   }

   FairRunAna *run = new FairRunAna();
   run->SetSource(new FairFileSource(inputFile));
   run->SetOutputFile(outputFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();
   FairParAsciiFileIo *parIo = new FairParAsciiFileIo();
   parIo->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo);

   AtPIDTask *pid = new AtPIDTask();
   pid->SetBField(bField);
   pid->SetInputBranch("AtPatternEvent");
   pid->SetOutputBranch("AtPIDEvent");
   pid->SetPersistence(kTRUE);
   run->AddTask(pid);

   TStopwatch t;
   t.Start();
   run->Init();
   run->Run(0, 0);
   t.Stop();
   std::cout << "\033[1;32mpid pass done.\033[0m -> " << outputFile << "  (Real " << t.RealTime() << " s)\n";
}
