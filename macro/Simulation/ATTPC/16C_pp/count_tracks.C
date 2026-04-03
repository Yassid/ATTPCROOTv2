void count_tracks(TString fname = "data/output_digi.root")
{
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      std::cout << "Cannot open " << fname << std::endl;
      return;
   }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) {
      std::cout << "No cbmsim tree" << std::endl;
      return;
   }
   TClonesArray *pa = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pa);
   int n0 = 0, n1 = 0, n2 = 0, nMore = 0;
   for (int i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      if (!pa || pa->GetEntries() == 0) {
         n0++;
         continue;
      }
      int nTr = ((AtPatternEvent *)pa->At(0))->GetTrackCand().size();
      if (nTr == 0) n0++;
      else if (nTr == 1) n1++;
      else if (nTr == 2) n2++;
      else nMore++;
   }
   std::cout << "0trk=" << n0 << " 1trk=" << n1 << " 2trk=" << n2 << " >2trk=" << nMore << std::endl;
}
