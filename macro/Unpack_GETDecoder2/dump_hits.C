// Dump reconstructed hits to CSV / total charge to a text file.
//
// This is the entry point of the whole python-side analysis (calib/*.py) and it used to
// exist only as ad-hoc `root -e` one-liners, so the CSVs those scripts read could not be
// regenerated from a fresh checkout. Written out properly here.
//
// Do not add #includes or gSystem->Load calls -- see rootlogon.C, which must load the
// dictionaries BEFORE this file is parsed.
//
// usage (from macro/Unpack_GETDecoder2, so that ./rootlogon.C is picked up):
//
//   root -l -b -q 'dump_hits.C("<in.root>","<out.csv>","hits",-1)'
//   root -l -b -q 'dump_hits.C("<in.root>","qtot.txt","qtot",-1)'
//
// modes
//   "hits"  per-hit CSV:  event,pad,tb,x,y,z,xc,yc,zc,q,qtot
//              x ,y ,z    AtHit::GetPosition()      -- raw pad coordinates
//              xc,yc,zc   AtHit::GetPositionCorr()  -- after the Lorentz/tilt correction
//              q          GetCharge()  peak height of the pulse
//              qtot       GetQHit()    integral of the pulse   <- use this for charge sums
//   "qtot"  one line per event: the event's summed GetQHit(), for trigger_eff.py
//
// nEvents < 0 means all of them.
void dump_hits(TString inFile, TString outFile, TString mode = "hits", Int_t nEvents = -1)
{
   TFile *f = TFile::Open(inFile);
   if (!f || f->IsZombie()) {
      std::cout << "ERROR cannot open " << inFile << std::endl;
      return;
   }
   TTree *tree = (TTree *)f->Get("cbmsim");
   if (!tree) {
      std::cout << "ERROR no cbmsim tree in " << inFile << std::endl;
      return;
   }

   // The hits live in a TClonesArray of AtEvent on branch AtEventH. NOTE: AtEvent is
   // ClassDef(3) on this branch and ClassDef(6) on OpenKF-Claude. Reading a file written
   // by the other branch gives ZERO hits with no error whatsoever -- if every event comes
   // out empty, check which branch produced the file before debugging anything else.
   TClonesArray *eventArray = nullptr;
   tree->SetBranchAddress("AtEventH", &eventArray);

   Long64_t nAll = tree->GetEntries();
   Long64_t nRun = (nEvents < 0 || nEvents > nAll) ? nAll : (Long64_t)nEvents;

   std::ofstream out(outFile.Data());
   if (!out.is_open()) {
      std::cout << "ERROR cannot write " << outFile << std::endl;
      return;
   }

   bool hitMode = (mode == "hits");
   if (hitMode)
      out << "event,pad,tb,x,y,z,xc,yc,zc,q,qtot\n";

   Long64_t nHitTot = 0, nEmpty = 0;
   for (Long64_t i = 0; i < nRun; i++) {
      tree->GetEntry(i);
      if (!eventArray || eventArray->GetEntriesFast() == 0) {
         nEmpty++;
         continue;
      }
      AtEvent *event = (AtEvent *)eventArray->At(0);
      if (!event) {
         nEmpty++;
         continue;
      }

      Int_t nHits = event->GetNumHits();
      if (nHits == 0)
         nEmpty++;

      Double_t qSum = 0.;
      for (Int_t j = 0; j < nHits; j++) {
         AtHit *hit = event->GetHit(j);
         if (!hit)
            continue;
         TVector3 p = hit->GetPosition();
         TVector3 pc = hit->GetPositionCorr();
         qSum += hit->GetQHit();
         nHitTot++;
         if (hitMode)
            out << i << "," << hit->GetHitPadNum() << "," << hit->GetTimeStamp() << ","
                << Form("%.3f,%.3f,%.3f", p.X(), p.Y(), p.Z()) << ","
                << Form("%.3f,%.3f,%.3f", pc.X(), pc.Y(), pc.Z()) << ","
                << Form("%.4g,%.4g", hit->GetCharge(), hit->GetQHit()) << "\n";
      }
      // Two columns: summed charge and hit multiplicity. The multiplicity is what
      // distinguishes a genuine low-charge event from a nearly-empty one, and matching this
      // layout keeps the file readable by the scripts written against the original dumps.
      if (!hitMode)
         out << Form("%.6g %d", qSum, nHits) << "\n";

      if (i % 2000 == 0)
         std::cout << "  " << i << " / " << nRun << std::endl;
   }
   out.close();

   std::cout << "DUMPED events=" << nRun << " hits=" << nHitTot << " empty_events=" << nEmpty << " -> " << outFile
             << std::endl;
   if (nHitTot == 0)
      std::cout << "WARNING zero hits. If the file was produced by the OpenKF-Claude branch "
                   "the AtEvent ClassDef differs and the hits cannot be streamed here."
                << std::endl;
}
