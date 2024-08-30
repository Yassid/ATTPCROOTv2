void tracking()
{

   FairRunAna *run = new FairRunAna();

   std::ofstream hitsFile;
   hitsFile.open("hits.txt");

   TFile *file = new TFile("test.root", "READ");
   TTree *tree = (TTree *)file->Get("cbmsim");
   Int_t nEvents = tree->GetEntries();
   std::cout << " Number of events : " << nEvents << std::endl;

   TTreeReader Reader("cbmsim", file);
   TTreeReaderValue<TClonesArray> eventHArray(Reader, "AtEventH");
   TTreeReaderValue<TClonesArray> patternArray(Reader, "AtPatternEvent");

   for (Int_t i = 0; i < nEvents; i++) {

      Reader.Next();

      AtEvent *event = (AtEvent *)eventHArray->At(0);
      AtPatternEvent *patternEvent = (AtPatternEvent *)patternArray->At(0);

      if (event && patternEvent) {

         auto &hitArray = event->GetHits();
         auto &tracks = patternEvent->GetTrackCand();
         std::cout << " Number of hits : " << hitArray.size() << std::endl;
         std::cout << " Number of tracks : " << tracks.size() << std::endl;

         /*for(auto &hit : hitArray) {
             auto pos = hit->GetPosition();
             auto charge = hit->GetCharge();
             auto time = hit->GetTimeStamp();
             hitsFile  << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
         }*/

         for (auto &track : tracks) {
            auto &points = track.GetHitArray();
            std::cout << "Track ID: " << track.GetTrackID() << std::endl;
            for (auto &point : points) {
               auto pos = point->GetPosition();
               auto charge = point->GetCharge();
               auto time = point->GetTimeStamp();
               std::cout << pos.X() << " " << pos.Y() << " " << pos.Z() << "  " << time << "  " << charge << std::endl;
            }
         }
      }
   }

   hitsFile.close();
   file->Close();
}