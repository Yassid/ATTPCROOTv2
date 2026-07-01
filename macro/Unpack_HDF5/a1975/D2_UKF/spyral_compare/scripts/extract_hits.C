// Dump ATTPCROOT corrected hits + PRA per-track hit assignment for a list of events
// to JSON, for the interactive point-cloud viewer.
#include <set>
#include <fstream>
void extract_hits(const char* evlist="/home/yassid/spyral_d2/viewer_events.json",
                  const char* in="/mnt/f/a1975/reco_d2/run_0016_reco.root",
                  const char* out="/home/yassid/spyral_d2/attpc_hits.json") {
   // parse wanted events from the json (simple: read integers after "events")
   std::set<long> want;
   {
      std::ifstream js(evlist); std::string s((std::istreambuf_iterator<char>(js)), {});
      auto p = s.find("\"events\"");
      if (p != std::string::npos) {
         auto lb = s.find('[', p), rb = s.find(']', p);
         std::string arr = s.substr(lb+1, rb-lb-1);
         std::stringstream ss(arr); std::string tok;
         while (std::getline(ss, tok, ',')) { try { want.insert(std::stol(tok)); } catch(...){} }
      }
   }
   printf("want %zu events\n", want.size());

   TFile* f = TFile::Open(in);
   TTree* t = (TTree*)f->Get("cbmsim");
   TClonesArray* corrArr = nullptr; TClonesArray* patArr = nullptr;
   t->SetBranchAddress("AtEventCorrected", &corrArr);
   t->SetBranchAddress("AtPatternEvent", &patArr);

   std::ofstream o(out); o << "{\n";
   bool firstEv = true;
   for (long ev : want) {
      if (ev < 0 || ev >= t->GetEntries()) continue;
      t->GetEntry(ev);
      if (!firstEv) o << ",\n"; firstEv = false;
      o << "\"" << ev << "\": {\n";
      // all corrected hits
      o << "  \"hits\": [";
      if (corrArr && corrArr->GetEntriesFast() > 0) {
         AtEvent* e = (AtEvent*)corrArr->At(0);
         int nh = e->GetNumHits(); bool fh = true;
         for (int i = 0; i < nh; ++i) {
            auto& h = e->GetHit(i);
            auto p = h.GetPosition();
            if (!fh) o << ","; fh = false;
            o << "[" << p.X() << "," << p.Y() << "," << p.Z() << "," << h.GetCharge() << "]";
         }
      }
      o << "],\n";
      // per-track hits
      o << "  \"tracks\": [";
      if (patArr && patArr->GetEntriesFast() > 0) {
         AtPatternEvent* pe = (AtPatternEvent*)patArr->At(0);
         auto& trks = pe->GetTrackCand(); bool ft = true;
         for (auto& tr : trks) {
            if (!ft) o << ","; ft = false;
            o << "{\"id\":" << tr.GetTrackID()
              << ",\"theta\":" << tr.GetGeoTheta()*TMath::RadToDeg()
              << ",\"radius\":" << tr.GetGeoRadius() << ",\"pts\":[";
            auto& hits = tr.GetHitArray(); bool fp = true;
            for (auto& hp : hits) {
               auto p = hp->GetPosition();
               if (!fp) o << ","; fp = false;
               o << "[" << p.X() << "," << p.Y() << "," << p.Z() << "," << hp->GetCharge() << "]";
            }
            o << "]}";
         }
      }
      o << "]\n}";
   }
   o << "\n}\n"; o.close();
   printf("DONE -> %s\n", out);
}
