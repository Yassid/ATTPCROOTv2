void testSOLARISUnpacker()
{
   TString dir = gSystem->Getenv("VMCWORKDIR");
   TString mapFile = "e12014_pad_mapping.xml";
   TString mapDir = dir + "/scripts/" + mapFile;
   auto fAtMapPtr = std::make_shared<AtTpcMap>();
   AtRawEvent rawEvent;

   TString dataFile = dir + "/macro/tests/MAGNEX/run_0372.txt";

   // Dummy event
   AtRawEvent event;

   std::unique_ptr<AtSOLARISUnpacker> unpacker = std::make_unique<AtSOLARISUnpacker>(
      fAtMapPtr, 10000, "/home/yassid/Desktop/NUMEN_data/tracker_and_sic_006_01_22642_000.sol");
   unpacker->Init();
   unpacker->FillRawEvent(event);
}
