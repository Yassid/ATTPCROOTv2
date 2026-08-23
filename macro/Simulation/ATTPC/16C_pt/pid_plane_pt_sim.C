/// @file pid_plane_pt_sim.C
/// @brief The Spyral PID plane (sqrt(dE/dx), Brho) for the SIMULATED 16C(p,t)14C, computed from
///        PATTERN TRACKS exactly as pid_plane_pt.C does for the data. NO FIT REQUIRED.
///
/// This is the plane the acceptance gate must be drawn on, and it has to be the SIMULATION's own.
/// The (p,d) channel measured what happens otherwise: applying the data-plane polygon
/// deuteron_tight.json to simulated deuterons keeps only 79.8%, and the loss is STATE-DEPENDENT
/// (75.7 / 77.2 / 80.5 / 85.9 for gs / ex1 / ex2 / ex3) -- a 10-point spread that would distort
/// the relative strengths it is supposed to correct. An equivalent polygon drawn on the sim plane
/// keeps 99.3%, spread 1.1 points. The 79.8% was never an efficiency; it was the OFFSET BETWEEN
/// THE TWO PLANES. So: draw here, on these points, and do not reuse triton_pt.json.
///
/// TRUTH IS CARRIED ALONGSIDE, which the data version cannot do. Each track gets `istriton`, set
/// by matching the pattern track's hits back to the MC track ID via AtPulseTask::SetSaveMCInfo.
/// Drawing is still yours -- a truth-selected gate would be circular, since the acceptance is
/// meant to measure how well a gate drawn on observables alone finds the tritons. But it lets you
/// SEE where the real ones sit while drawing, and it lets the drawn gate be scored afterwards
/// (purity and efficiency against truth), which is a check the data plane can never offer.
///
/// The five states are kept separate in `state` because their kinematics differ -- though far less
/// than in (d,t): 14C has nothing below 6.09 MeV, and across the 6.09-7.01 cluster the triton lab
/// angle moves by under a degree. So a state dependence here would be a surprise, and that is
/// exactly why it is measured rather than assumed.
///
/// AND THE REASON THIS CHANNEL NEEDS IT MOST: 16C(p,t)14C has no measured signal. The excitation
/// peak walks ~11.5 MeV between theta_lab 8 and 40 deg -- the signature of assigning triton mass
/// to particles that are not tritons -- and the existing gate (pt/data/triton_pt.json in the
/// analysis repo) was drawn by eye on a plane with no truth behind it, on a sample now measured to
/// be 31% oxygen-induced. With an expected (p,t) yield of only 1-10% of (p,d), sitting under
/// (p,p) and (p,d) channels 10-100x more numerous, gate PURITY is the whole ballgame. This is the
/// first time it can be measured.
///
///   root -b -q 'pid_plane_pt_sim.C()'                 // all five states
///   root -b -q 'pid_plane_pt_sim.C("gs_s4001")'       // one
void pid_plane_pt_sim(TString statesCSV = "gs_s4001,ex1_s4001,ex2_s4001,ex3_s4001,ex4_s4001",
                      TString simDir = "/mnt/f/a1975_C16_pt_sim/",
                      TString outCache = "data/pid_plane_pt_sim.root",
                      TString outPng = "data/pid_plane_pt_sim.png", double bField = 2.85,
                      double sqrtMax = 60, double brhoMax = 3.0, Long64_t maxEvt = -1)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   auto *h = new TH2F("hpid", "a1975 H2 SIMULATION, 16C(p,t)14C, Spyral PID;#sqrt{dE/dx};B#rho [T m]", 300, 0, sqrtMax,
                      300, 0, brhoMax);
   auto *ht = new TH2F("hpid_t", "tritons only (MC truth);#sqrt{dE/dx};B#rho [T m]", 300, 0, sqrtMax,
                       300, 0, brhoMax);

   TFile *fo = TFile::Open(outCache, "RECREATE");
   auto *tp = new TTree("pts", "Spyral PID points, simulated pattern tracks");
   float x, y, polar, vz, vr, ncl;
   int istriton, stateIdx;
   // Same branch names as pid_plane_pt.C so the gate drawer and every downstream tool work on
   // either file without a special case. `ic` is absent on purpose: the simulation has no ion
   // chamber, and inventing one would let a beam gate be "applied" to a quantity that is not
   // measured. The (d,t) acceptance run on Spyral hit exactly that -- every simulated estimate
   // carried ic = -1, so the data's 900-1400 selection passed 0%.
   tp->Branch("sqrtdedx", &x);
   tp->Branch("brho", &y);
   tp->Branch("polar", &polar);
   tp->Branch("vertexz", &vz);
   tp->Branch("vertexr", &vr);
   tp->Branch("ncl", &ncl);
   tp->Branch("istriton", &istriton);
   tp->Branch("state", &stateIdx);

   long nTrk = 0, nVal = 0, nTri = 0;
   TObjArray *ta = statesCSV.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      stateIdx = it;
      TString fr = simDir + tg + "_reco.root";
      if (gSystem->AccessPathName(fr)) {
         printf("  skip %s (no %s)\n", tg.Data(), fr.Data());
         continue;
      }
      TFile *fi = TFile::Open(fr);
      TTree *t = fi ? (TTree *)fi->Get("cbmsim") : nullptr;
      if (!t) {
         printf("  skip %s (no tree)\n", tg.Data());
         if (fi) fi->Close();
         continue;
      }
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pe);

      Long64_t N = (maxEvt > 0) ? std::min(maxEvt, t->GetEntries()) : t->GetEntries();
      long v0 = nVal, t0 = nTri;
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (!pe || pe->GetEntries() == 0)
            continue;
         auto *p = (AtPatternEvent *)pe->At(0);
         if (!p)
            continue;
         for (auto &trk : p->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(trk);
            ++nTrk;
            auto r = spy.Estimate(tr);
            if (!r.valid)
               continue;
            ++nVal;
            // MC truth. Identify the triton by its SPECIES (A = 3, Z = 1), not by a track ID:
            // an ID is an assumption about the order the generator stacked the particles, which
            // is exactly the kind of thing that silently changes and then mislabels an entire
            // sample. Majority vote over the track's hits rather than the first hit, so one
            // shared cluster cannot flip the label.
            int nT = 0, nTot = 0;
            for (const auto &hit : tr.GetHitArray()) {
               if (!hit)
                  continue;
               const auto &mcs = hit->GetMCSimPointArray();
               if (mcs.empty())
                  continue;
               ++nTot;
               if (mcs[0].A == 3 && mcs[0].Z == 1)
                  ++nT;
            }
            istriton = (nTot > 0 && nT > nTot / 2) ? 1 : 0;
            if (istriton)
               ++nTri;

            x = r.sqrtdEdx;
            y = r.brho;
            polar = r.polar * TMath::RadToDeg();
            vz = r.vertex.Z();
            vr = std::sqrt(r.vertex.X() * r.vertex.X() + r.vertex.Y() * r.vertex.Y());
            ncl = r.nClusters;
            h->Fill(x, y);
            if (istriton)
               ht->Fill(x, y);
            tp->Fill();
         }
      }
      printf("  %-14s %8lld evt   valid %7ld   truth-triton %7ld\n", tg.Data(), N, nVal - v0, nTri - t0);
      fi->Close();
   }
   printf("=== pattern tracks %ld, valid Spyral %ld (%.1f%%), of which truth tritons %ld (%.1f%%) ===\n",
          nTrk, nVal, nTrk ? 100.0 * nVal / nTrk : 0.0, nTri, nVal ? 100.0 * nTri / nVal : 0.0);

   fo->cd();
   tp->Write();
   h->Write();
   ht->Write();
   fo->Close();
   printf("cache -> %s\n", outCache.Data());

   auto *c = new TCanvas("cpid", "pid", 1400, 700);
   c->Divide(2, 1);
   c->cd(1)->SetLogz();
   h->Draw("colz");
   c->cd(2)->SetLogz();
   ht->Draw("colz");
   c->SaveAs(outPng);
   printf("plot  -> %s   (left: everything, right: truth tritons)\n", outPng.Data());
}
