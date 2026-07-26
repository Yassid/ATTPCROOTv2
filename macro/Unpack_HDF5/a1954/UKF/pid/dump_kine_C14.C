/// @file dump_kine_C14.C
/// @brief Dump per-fitted-track kinematics + PID + IC to CSV, for the gated 14C excitation
///        spectrum (done in python with the drawn gates). Works for UKF (_ukf) or GENFIT
///        (_genfit) fit files. IC matched by event ID from <run>_FRIB.root.
///
///   root -b -q 'dump_kine_C14.C("run_0196,run_0197","_ukf","/home/yassid/a1954_C14_reco_hdb/","/home/yassid/a1954_C14_reco/","/tmp/kine_ukf.csv")'
///
/// Columns: ke,theta,chi2ndf,sqrtdedx,brho,vertexz,ic,arclen,npts
void dump_kine_C14(TString runsCSV, TString fitSuffix, TString fitDir, TString fribDir, TString outCsv,
                    Int_t icTbLo = 1000, Int_t icTbHi = 1350, double bField = 2.85)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   std::ofstream csv(outCsv.Data());
   csv << "ke,theta,chi2ndf,sqrtdedx,brho,vertexz,ic,arclen,npts\n";

   TObjArray *runs = runsCSV.Tokenize(",");
   long nTrk = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff = fitDir + run + fitSuffix + ".root";
      if (gSystem->AccessPathName(ff)) {
         printf("skip %s (no %s)\n", run.Data(), ff.Data());
         continue;
      }
      // IC map from FRIB (event ID -> amplitude)
      std::map<int, double> icByID;
      TString frf = fribDir + run + "_FRIB.root";
      if (!gSystem->AccessPathName(frf)) {
         TFile *fF = TFile::Open(frf);
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

      // The fit output does NOT propagate the event ID, but its entries are aligned 1:1 with
      // <run>_reco.root (same FairRoot event order). Read the event ID from _reco entry i.
      TFile *fT = TFile::Open(ff);
      TTree *tT = (TTree *)fT->Get("cbmsim");
      TClonesArray *te = nullptr;
      tT->SetBranchAddress("AtTrackingEvent", &te);
      TFile *fR = TFile::Open(fitDir + run + "_reco.root");
      TTree *tR = (fR ? (TTree *)fR->Get("cbmsim") : nullptr);
      TClonesArray *ev0 = nullptr;
      if (tR)
         tR->SetBranchAddress("AtEventH", &ev0);
      Long64_t N = tR ? std::min(tT->GetEntries(), tR->GetEntries()) : tT->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         tT->GetEntry(i);
         if (te->GetEntries() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         double ic = -1;
         if (tR) {
            tR->GetEntry(i);
            if (ev0->GetEntries() > 0) {
               int id = (int)((AtEvent *)ev0->At(0))->GetEventID();
               auto it = icByID.find(id);
               if (it != icByID.end())
                  ic = it->second;
            }
         }
         // original PRA tracks for PID, keyed by track ID
         std::vector<AtTrack> orig = ev->GetTrackArray();
         std::map<int, AtTrack *> byID;
         for (auto &tr : orig)
            byID[tr.GetTrackID()] = &tr;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double c2n = ndf > 0 ? chi2 / ndf : 1e9;
            double ke = k.kineticEnergy, thDeg = k.theta * TMath::RadToDeg();
            if (ke <= 0 || ke > 1000)
               continue;
            double sqrtdedx = -1, brho = -1, vz = -1, arclen = -1;
            int npts = 0;
            auto pit = byID.find(ft->GetTrackID());
            if (pit != byID.end()) {
               auto r = spy.Estimate(*pit->second);
               if (r.valid) {
                  sqrtdedx = r.sqrtdEdx;
                  brho = r.brho;
                  vz = r.vertex.Z();
                  arclen = r.arclength;
                  npts = r.nPoints;
               }
            }
            csv << ke << "," << thDeg << "," << c2n << "," << sqrtdedx << "," << brho << "," << vz << "," << ic << ","
                << arclen << "," << npts << "\n";
            ++nTrk;
         }
      }
      fT->Close();
      printf("dumped %s\n", run.Data());
   }
   csv.close();
   printf("wrote %ld fitted tracks -> %s\n", nTrk, outCsv.Data());
}
