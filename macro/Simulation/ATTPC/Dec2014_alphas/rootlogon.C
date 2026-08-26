// ROOT runs this before any macro in this directory.
//
// The dictionaries must be loaded HERE, before the macro is parsed. Loading them inside a
// macro function is too late (cling parses the whole function first, so every class comes
// out "unknown type name"), and #including the headers instead makes cling re-parse
// FairRoot's Base dictionary payload out of order. Same arrangement as
// macro/Unpack_GETDecoder2/rootlogon.C.
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
   gSystem->Load("libAtTpc");          // AtTpcPoint (MC truth) lives here
   gSystem->Load("libAtDigitization");
}
