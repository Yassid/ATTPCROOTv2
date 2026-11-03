// Dump ATTPCROOT's pad-centre table so it can be compared pad-for-pad with Spyral's padxy.csv.
void padmap(){
  auto m = std::make_shared<AtTpcMap>();
  TString dir = gSystem->Getenv("VMCWORKDIR");
  m->ParseXMLMap((dir + "/scripts/e12014_pad_mapping.xml").Data());
  m->GeneratePadPlane();
  FILE* o = fopen("/tmp/claude-1000/-home-yassid/c1d6c7d6-78db-4f23-a0ea-cf85e98aff36/scratchpad/attpc_padxy.csv","w");
  fprintf(o,"pad,x,y\n");
  int n=0;
  for(int p=0;p<10240;p++){ auto c=m->CalcPadCenter(p);
    if(std::isnan(c.X())) continue; fprintf(o,"%d,%.5f,%.5f\n",p,c.X(),c.Y()); n++; }
  fclose(o); printf("pads written: %d\n", n);
}
