// Complete events in a run = min over CoBos of (last event ID + 1).
// GetBasicFrame(-1) returns the NEXT frame, so walk to the end like AtGRAWUnpacker does.
void count_events(const char *listfile)
{
   std::ifstream in(listfile);
   std::string line;
   std::map<int, std::vector<std::string>> byCobo;
   while (std::getline(in, line)) {
      if (line.empty()) continue;
      auto p = line.rfind("CoBo");
      if (p == std::string::npos) continue;
      byCobo[atoi(line.c_str() + p + 4)].push_back(line);
   }
   int minEvents = INT_MAX;
   for (auto &kv : byCobo) {
      GETDecoder2 d;
      for (auto &f : kv.second) d.AddData(f.c_str());
      d.SetData(0);
      d.SetPseudoTopologyFrame(15, kFALSE);
      int last = -1;
      while (GETBasicFrame *bf = d.GetBasicFrame(-1)) last = bf->GetEventID();
      if (last + 1 < minEvents) minEvents = last + 1;
   }
   std::cout << "NEVENTS " << (minEvents == INT_MAX ? 0 : minEvents) << std::endl;
}
