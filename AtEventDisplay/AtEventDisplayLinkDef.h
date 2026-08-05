#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// Only the AT-TPC display is built PCL-free; the Proto and S800 variants are excluded.
#pragma link C++ class AtEventManager + ;
#pragma link C++ class AtEventDrawTask + ;

#endif
