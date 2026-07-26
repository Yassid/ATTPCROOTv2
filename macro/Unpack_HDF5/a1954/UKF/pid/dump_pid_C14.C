/// @file dump_pid_C14.C
/// @brief Dump AtSpyralPID observables per track from <run>_reco.root to a CSV, for the
///        Spyral-1.0.0 vs ATTPCROOT (Triplclust/HDBSCAN) PID comparison. Optionally
///        applies the 14C IC beam gate (event-ID matched to <run>_FRIB.root).
///
///   root -b -q 'dump_pid_C14.C("run_0147,run_0148","/home/yassid/a1954_C14_reco/","/tmp/pid_tc.csv")'
///   with IC gate: pass fribDir + icLo/icHi (>0) -> only 14C-beam tracks are written.
void dump_pid_C14(TString runsCSV, TString recoDir, TString outCsv, TString fribDir = "", double icLo = -1,
                   double icHi = -1, Int_t icTbLo = 1000, Int_t icTbHi = 1350, double bField = 2.85,
                   bool useTraceIntegral = false)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);
   spy.SetUseTraceIntegral(useTraceIntegral); // true -> integrated charge (Spyral dE/dx scale)
   bool useIC = (fribDir.Length() > 0);        // record IC value if FRIB available
   bool hardGate = (useIC && icLo >= 0 && icHi > 0); // additionally skip out-of-gate tracks

   std::ofstream csv(outCsv.Data());
   csv << "sqrtdedx,brho,dedx,polar,arclen,npts,ic\n";

   TObjArray *runs = runsCSV.Tokenize(",");
   long nTrk = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString rf = recoDir + run + "_reco.root";
      if (gSystem->AccessPathName(rf)) {
         printf("skip %s\n", run.Data());
         continue;
      }
      std::map<int, double> icByID;
      if (useIC) {
         TString ff = fribDir + run + "_FRIB.root";
         if (!gSystem->AccessPathName(ff)) {
            TFile *fF = TFile::Open(ff);
            TTree *tF = (TTree *)fF->Get("cbmsim");
            TClonesArray *ra = nullptr;
            tF->SetBranchAddress("AtRawEvent", &ra);
            for (Long64_t i = 0; i < tF->GetEntries(); ++i) {
               tF->GetEntry(i);
               if (ra->GetEntries() == 0)
                  continue;
               auto *raw = (AtRawEvent *)ra->At(0);
               if (!raw || raw->GetGenTraces().empty())
                  continue;
               auto &adc = raw->GetGenTraces()[0]->GetADC();
               double mx = -1e9;
               for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
                  mx = std::max(mx, (double)adc[b]);
               icByID[(int)raw->GetEventID()] = mx;
            }
            fF->Close();
         }
      }

      TFile *fR = TFile::Open(rf);
      TTree *tR = (TTree *)fR->Get("cbmsim");
      TClonesArray *pe = nullptr, *ev = nullptr;
      tR->SetBranchAddress("AtPatternEvent", &pe);
      tR->SetBranchAddress("AtEventH", &ev);
      for (Long64_t i = 0; i < tR->GetEntries(); ++i) {
         tR->GetEntry(i);
         double ic = -1;
         if (useIC) {
            if (ev->GetEntries() == 0)
               continue;
            int id = (int)((AtEvent *)ev->At(0))->GetEventID();
            auto it = icByID.find(id);
            if (it == icByID.end())
               continue;
            ic = it->second;
            if (hardGate && (ic < icLo || ic > icHi))
               continue;
         }
         if (pe->GetEntries() == 0)
            continue;
         auto *p = (AtPatternEvent *)pe->At(0);
         if (!p)
            continue;
         for (auto &trk : p->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(trk);
            auto r = spy.Estimate(tr);
            if (!r.valid)
               continue;
            csv << r.sqrtdEdx << "," << r.brho << "," << r.dEdx << "," << r.polar * TMath::RadToDeg() << ","
                << r.arclength << "," << r.nPoints << "," << ic << "\n";
            ++nTrk;
         }
      }
      fR->Close();
      printf("dumped %s\n", run.Data());
   }
   csv.close();
   printf("wrote %ld tracks -> %s\n", nTrk, outCsv.Data());
}
