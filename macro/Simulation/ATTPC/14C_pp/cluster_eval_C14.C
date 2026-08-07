/// @file cluster_eval_C14.C
/// @brief Score the HDBSCAN clustering of the 14C(p,p') simulation AGAINST MC TRUTH.
///
/// The viewer showed single protons chopped into several PRA tracks. This quantifies that:
/// every digitized hit carries AtHit::MCSimPoint (AtPulseTask::SetSaveMCInfo, on in
/// run_reco_C14.C), so each hit has a true MC trackID and the clustering can be scored the
/// way a clustering is normally scored.
///
/// Per TRUE particle (>= minTruth labelled hits in the event):
///   fragments  = number of clusters holding >= fragFrac of its hits   (1 = perfect)
///   efficiency = largest single cluster's share of its hits           (1 = perfect)
///   lost       = share of its hits that no cluster claimed
/// Per CLUSTER:
///   purity     = share of its labelled hits coming from one true particle
///
///   root -b -q 'cluster_eval_C14.C("./data/sim_reco.root","orig")'
///   root -b -q 'cluster_eval_C14.C("./diagnostics/scan/T3_reco.root","T3")'

void cluster_eval_C14(TString inFile = "./data/sim_reco.root", TString tag = "orig", Long64_t maxEvt = -1,
                      Int_t minTruth = 20, Double_t fragFrac = 0.10)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);

   TFile *f = TFile::Open(inFile);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", inFile.Data());
      return;
   }
   auto *t = (TTree *)f->Get("cbmsim");
   if (!t || !t->GetBranch("AtEventH") || !t->GetBranch("AtPatternEvent")) {
      printf("\033[1;31m%s needs BOTH AtEventH and AtPatternEvent (MC info must be on)\033[0m\n", inFile.Data());
      return;
   }
   TClonesArray *ev = nullptr, *pe = nullptr;
   t->SetBranchAddress("AtEventH", &ev);
   t->SetBranchAddress("AtPatternEvent", &pe);

   TH1F *hFrag = new TH1F("hFrag", Form("%s;clusters holding a true track;true tracks", tag.Data()), 12, 0, 12);
   TH1F *hEff = new TH1F("hEff", Form("%s;best-cluster efficiency;true tracks", tag.Data()), 50, 0, 1.001);
   TH1F *hPur = new TH1F("hPur", Form("%s;cluster purity;clusters", tag.Data()), 50, 0, 1.001);
   TH1F *hLost = new TH1F("hLost", Form("%s;unclustered share of a true track;true tracks", tag.Data()), 50, 0, 1.001);

   Long64_t nE = t->GetEntries();
   if (maxEvt > 0 && maxEvt < nE)
      nE = maxEvt;
   long nTrue = 0, nClus = 0, nEvt = 0, nPerfect = 0, nSplit = 0, nLostAll = 0;
   long nTrueP = 0, nPerfectP = 0, nSplitP = 0, nMerged = 0, nMergedP = 0;
   double sumEff = 0, sumPur = 0, sumLost = 0, sumFrag = 0;
   double sumEffP = 0, sumLostP = 0, sumFragP = 0, sumPurP = 0;
   double hitsLabelled = 0, hitsClustered = 0;

   for (Long64_t i = 0; i < nE; ++i) {
      t->GetEntry(i);
      if (!ev || !ev->GetEntriesFast() || !pe || !pe->GetEntriesFast())
         continue;
      auto *e = (AtEvent *)ev->At(0);
      auto *p = (AtPatternEvent *)pe->At(0);
      if (!e || !p)
         continue;
      ++nEvt;

      // truth: hitID -> MC trackID, the size of each true track, and whether it is the proton.
      // The beam (Z=6,A=14) is a long straight easy track and would flatter every number, so
      // it is scored separately from the recoil proton (Z=1,A=1) that the analysis needs.
      std::map<int, int> truthOf;
      std::map<int, long> truthSize;
      std::map<int, bool> isProton;
      for (int h = 0; h < e->GetNumHits(); ++h) {
         auto &hit = e->GetHit(h);
         const auto &mc = hit.GetMCSimPointArray();
         if (mc.empty())
            continue; // unlabelled (noise / no MC contributor)
         truthOf[hit.GetHitID()] = mc[0].trackID; // first point = dominant contributor
         ++truthSize[mc[0].trackID];
         isProton[mc[0].trackID] = (mc[0].Z == 1 && mc[0].A == 1);
         ++hitsLabelled;
      }

      // clusters: for each PRA track, how its labelled hits distribute over true tracks
      std::vector<std::map<int, long>> clus; // per cluster: truthID -> count
      for (auto &trk : p->GetTrackCand()) {
         AtTrack &tr = const_cast<AtTrack &>(trk);
         std::map<int, long> m;
         for (auto &h : tr.GetHitArray()) { // HitVector holds POINTERS
            if (!h)
               continue;
            auto it = truthOf.find(h->GetHitID());
            if (it == truthOf.end())
               continue;
            ++m[it->second];
            ++hitsClustered;
         }
         if (m.empty())
            continue;
         clus.push_back(m);
         ++nClus;
         long tot = 0, best = 0;
         for (auto &kv : m) {
            tot += kv.second;
            best = std::max(best, kv.second);
         }
         double pur = double(best) / tot;
         hPur->Fill(pur);
         sumPur += pur;
      }

      // score every true track big enough to be findable
      for (auto &kv : truthSize) {
         if (kv.second < minTruth)
            continue;
         ++nTrue;
         long best = 0, claimed = 0;
         int frag = 0;
         double bestPur = 0;   // purity OF THE CLUSTER that best covers this true track
         for (auto &c : clus) {
            auto it = c.find(kv.first);
            if (it == c.end())
               continue;
            claimed += it->second;
            if (it->second > best) {
               best = it->second;
               long tot = 0;
               for (auto &kv2 : c)
                  tot += kv2.second;
               bestPur = double(it->second) / tot;   // <1 means the cluster ALSO holds other particles
            }
            if (double(it->second) / kv.second >= fragFrac)
               ++frag;
         }
         double eff = double(best) / kv.second;
         double lost = 1.0 - double(claimed) / kv.second;
         hFrag->Fill(frag);
         hEff->Fill(eff);
         hLost->Fill(lost);
         sumEff += eff;
         sumLost += lost;
         sumFrag += frag;
         // CLEAN = one cluster, holding nearly all of this track AND nearly nothing else.
         // Requiring only (frag==1 && eff>0.9) rewards a cluster that swallows every particle
         // in the event -- that bug made an all-merging clustering look 95.6 % clean.
         if (frag == 1 && eff > 0.9 && bestPur > 0.9)
            ++nPerfect;
         if (frag >= 2)
            ++nSplit;
         if (claimed == 0)
            ++nLostAll;
         if (bestPur > 0 && bestPur <= 0.9)
            ++nMerged;
         if (isProton[kv.first]) {
            ++nTrueP;
            sumEffP += eff;
            sumFragP += frag;
            sumLostP += lost;
            sumPurP += bestPur;
            if (frag == 1 && eff > 0.9 && bestPur > 0.9)
               ++nPerfectP;
            if (frag >= 2)
               ++nSplitP;
            if (bestPur > 0 && bestPur <= 0.9)
               ++nMergedP;
         }
      }
   }

   printf("\n===== clustering vs MC truth : %s =====\n", tag.Data());
   printf("events %ld   true tracks (>=%d hits) %ld   clusters %ld   clusters/true %.2f\n", nEvt, minTruth, nTrue,
          nClus, nTrue ? double(nClus) / nTrue : 0);
   if (!nTrue) {
      printf("\033[1;31mno labelled truth -- was AtPulseTask::SetSaveMCInfo() on?\033[0m\n");
      return;
   }
   printf("mean fragments/true track : %.2f   (1.00 = never split)\n", sumFrag / nTrue);
   printf("SPLIT   (>=2 fragments)   : %5.1f %%\n", 100.0 * nSplit / nTrue);
   printf("MERGED  (best cluster <90%% pure) : %5.1f %%\n", 100.0 * nMerged / nTrue);
   printf("CLEAN   (1 frag, eff>0.9, pure>0.9) : %5.1f %%\n", 100.0 * nPerfect / nTrue);
   printf("mean best-cluster eff     : %.3f\n", sumEff / nTrue);
   printf("mean unclustered share    : %.3f   (fully lost tracks %.1f %%)\n", sumLost / nTrue,
          100.0 * nLostAll / nTrue);
   printf("mean cluster purity       : %.3f\n", nClus ? sumPur / nClus : 0);
   printf("labelled hits placed in a cluster: %.1f %%\n", hitsLabelled ? 100.0 * hitsClustered / hitsLabelled : 0);
   printf("\n-- RECOIL PROTON ONLY (the track the analysis needs) --\n");
   if (nTrueP) {
      printf("true protons %ld   mean fragments %.2f   SPLIT %.1f %%   CLEAN %.1f %%\n", nTrueP, sumFragP / nTrueP,
             100.0 * nSplitP / nTrueP, 100.0 * nPerfectP / nTrueP);
      printf("MERGED with another particle %.1f %%\n", 100.0 * nMergedP / nTrueP);
      printf("mean best-cluster eff %.3f   purity of that cluster %.3f   unclustered %.3f\n", sumEffP / nTrueP,
             sumPurP / nTrueP, sumLostP / nTrueP);
   } else
      printf("none found\n");

   TCanvas *c = new TCanvas("c", "clustering", 1500, 420);
   c->Divide(4, 1);
   c->cd(1); hFrag->SetFillColor(kOrange + 7); hFrag->Draw();
   c->cd(2); hEff->SetFillColor(kAzure + 1); hEff->Draw();
   c->cd(3); hPur->SetFillColor(kSpring - 6); hPur->Draw();
   c->cd(4); hLost->SetFillColor(kGray + 1); hLost->Draw();
   TString png = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/diagnostics/cluster_" + tag + ".png";
   c->SaveAs(png);
   printf("saved %s\n\n", png.Data());
   f->Close();
}
