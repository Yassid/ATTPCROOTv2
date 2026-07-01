// Dump corrected hits + PRA-seed assignment for the follower validation set.
// One CSV row per hit: event,x,y,z,q,hitID,praTrk (-1 if PRA left it unassigned).
#include <set>
#include <map>
#include <fstream>
void extract_follower(const char* evlist="/home/yassid/spyral_d2/follower_events.json",
                      const char* in="/mnt/f/a1975/reco_d2/run_0016_reco.root",
                      const char* out="/home/yassid/spyral_d2/follower_hits.csv") {
   std::set<long> want;
   { std::ifstream js(evlist); std::string s((std::istreambuf_iterator<char>(js)), {});
     auto p = s.find("\"events\""); auto lb = s.find('[', p), rb = s.find(']', p);
     std::stringstream ss(s.substr(lb+1, rb-lb-1)); std::string tok;
     while (std::getline(ss, tok, ',')) { try { want.insert(std::stol(tok)); } catch(...){} } }
   printf("want %zu events\n", want.size());

   TFile* f = TFile::Open(in);
   TTree* t = (TTree*)f->Get("cbmsim");
   TClonesArray *corr = nullptr, *pat = nullptr;
   t->SetBranchAddress("AtEventCorrected", &corr);
   t->SetBranchAddress("AtPatternEvent", &pat);

   std::ofstream o(out); o << "event,x,y,z,q,hitID,praTrk\n";
   for (long ev : want) {
      if (ev < 0 || ev >= t->GetEntries()) continue;
      t->GetEntry(ev);
      // map hitID -> PRA track index
      std::map<int,int> hid2trk;
      if (pat && pat->GetEntriesFast() > 0) {
         auto& trks = ((AtPatternEvent*)pat->At(0))->GetTrackCand();
         for (size_t k = 0; k < trks.size(); ++k)
            for (auto& hp : trks[k].GetHitArray()) hid2trk[hp->GetHitID()] = k;
      }
      if (!corr || corr->GetEntriesFast() == 0) continue;
      AtEvent* e = (AtEvent*)corr->At(0);
      for (int i = 0; i < e->GetNumHits(); ++i) {
         auto& h = e->GetHit(i); auto p = h.GetPosition();
         int id = h.GetHitID();
         int tk = hid2trk.count(id) ? hid2trk[id] : -1;
         o << ev << "," << p.X() << "," << p.Y() << "," << p.Z() << ","
           << h.GetCharge() << "," << id << "," << tk << "\n";
      }
   }
   o.close();
   printf("DONE -> %s\n", out);
}
