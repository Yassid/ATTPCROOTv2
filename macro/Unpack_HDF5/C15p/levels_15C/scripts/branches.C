void branches(){
  const char* F[2]={"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root",
                    "/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root"};
  for(auto f:F){ TFile* t=TFile::Open(f); if(!t||t->IsZombie()){printf("%s MISSING\n",f);continue;}
    printf("\n%s\n",f); t->ls();
    TTree* d=(TTree*)t->Get("dd");
    if(d){ printf("  tree 'dd', %lld entries, branches:\n   ",d->GetEntries());
      TObjArray* b=d->GetListOfBranches();
      for(int i=0;i<b->GetEntries();i++) printf(" %s",b->At(i)->GetName());
      printf("\n"); }
    t->Close(); }
}
