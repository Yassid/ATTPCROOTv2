void testUnpacker()
{
   TString dir = gSystem->Getenv("VMCWORKDIR");
   TString mapFile = "e12014_pad_mapping.xml";
   TString mapDir = dir + "/scripts/" + mapFile;
   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   AtRawEvent rawEvent;

   TString dataFile = "/home/software_dev/fair_install/ATTPCROOTv2/macro/tests/MAGNEX/run_0372.txt";

   std::unique_ptr<AtV2745Unpacker> unpacker = std::make_unique<AtV2745Unpacker>(fAtMapPtr, 10000);
   unpacker->SetInputFileName(dataFile.Data());
   unpacker->Init();
}
