/// Dump, per track: PRA GeoTheta (which sets backwardSeed and therefore WHICH END is called the
/// vertex), and the charge profile, which says which end the particle actually stopped at.
/// If they disagree, AtGenfitter seeds from the stopping end and reads out the energy the
/// ejectile had AFTER losing most of it -- the exact failure the backward-seed fix was written
/// for, but arising from a wrong TAG rather than a missing fix.
void dump_direction_dp(TString listFile="/mnt/f/a1975/caches/.explorer/refit_list.txt",
                       TString outCsv="/tmp/dir_dp.csv",
                       TString gfDir="/mnt/f/a1975/gf_dp_cateloss/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   std::ifstream in(listFile.Data()); std::string r; Long64_t e; int ti;
   std::vector<std::tuple<TString,Long64_t,int>> jobs;
   while (in >> r >> e >> ti) jobs.push_back({TString("run_")+r.c_str()+"_multifit", e, ti});
   std::ofstream o(outCsv.Data());
   o << "run,entry,tid,ncl,geotheta,zlab_first,zlab_last,q_first20,q_mid,q_last15\n";
   TString open; TFile *f=nullptr; TTree *t=nullptr; TClonesArray *te=nullptr;
   for (auto &j : jobs) {
      TString rt=std::get<0>(j); Long64_t en=std::get<1>(j); int id=std::get<2>(j);
      if (rt!=open) { if(f) f->Close(); f=TFile::Open(gfDir+rt+"_genfitter_p.root");
         if(!f||f->IsZombie()){f=nullptr;continue;} t=(TTree*)f->Get("cbmsim"); te=nullptr;
         t->SetBranchAddress("AtTrackingEvent",&te); open=rt; }
      t->GetEntry(en);
      auto *ev=(AtTrackingEvent*)te->At(0); if(!ev) continue;
      auto trks=ev->GetTrackArray(); AtTrack *src=nullptr;
      for (auto &tr:trks) if (tr.GetTrackID()==id) src=&tr;
      if(!src) continue;
      auto *hc=src->GetHitClusterArray(); int n=hc->size(); if(n<20) continue;
      int a=n/5, b=(3*n)/5, c=(17*n)/20;
      auto meanq=[&](int i0,int i1){ double s=0; int m=0;
         for(int i=i0;i<i1&&i<n;++i){s+=(*hc)[i].GetCharge();++m;} return m?s/m:0.0; };
      o << rt << "," << en << "," << id << "," << n << ","
        << src->GetGeoTheta()*TMath::RadToDeg() << ","
        << 1000.0-(*hc)[0].GetPosition().Z() << "," << 1000.0-(*hc)[n-1].GetPosition().Z() << ","
        << meanq(0,a) << "," << meanq(b,c) << "," << meanq(c,n) << "\n";
   }
   if(f) f->Close();
   printf("  wrote %s\n", outCsv.Data());
}
