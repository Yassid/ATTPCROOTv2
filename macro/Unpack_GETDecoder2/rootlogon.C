// ROOT executes ./rootlogon.C before any macro on the command line.
//
// This branch's .rootmap autoloading does not work, so the dictionaries have to be
// loaded here, BEFORE the macro is parsed. Loading them from inside a macro function
// is too late -- cling parses the whole function first, so every class comes out as
// "unknown type name". Including the headers instead makes cling re-parse FairRoot's
// Base dictionary payload out of order ("field has incomplete type 'FairRunInfo'").
{
   gSystem->Load("libXMLParser");

   // FairRoot
   gSystem->Load("libFairTools");
   gSystem->Load("libParBase");
   gSystem->Load("libGeoBase");
   gSystem->Load("libBase");

   // ATTPCROOT
   gSystem->Load("libAtTpcMap");
   gSystem->Load("libAtData");
   gSystem->Load("libAtParameter");
   gSystem->Load("libAtTools");
   gSystem->Load("libAtReconstruction");

   // Event display (PCL-free build). Eve and FairRoot's EventDisplay must come first:
   // AtEventManager derives from TEveEventManager.
   gSystem->Load("libEve");
   gSystem->Load("libEventDisplay");
   gSystem->Load("libAtEventDisplay");
}
