/// @file pid_plane_dt.C
/// @brief Spyral PID plane (sqrt(dE/dx), Brho) for the a1975 D2 target, computed from the
///        PATTERN TRACKS in <run>_multifit_reco.root. NO FIT REQUIRED.
///
/// WHY THIS EXISTS.  pid_draw_gate_d2.C fills the plane from the AtPIDEvent branch inside the
/// <run>_genfitter_*.root files, i.e. it needs a fit to have run first -- and the fit needs a
/// gate, which is what we are trying to draw.  That circle is avoidable: AtSpyralPID::Estimate
/// works on an AtTrack straight out of AtPatternEvent and never looks at the fit.  This is the
/// same fact that lets a gate be redrawn without refitting (verified on the (p,d) channel:
/// three independent computation paths agreed bit for bit over 5101 tracks).
///
/// WHY THE OLD GATE CANNOT BE REUSED.  triton_d2.json was drawn on the dv 1.136 plane.
/// AtSpyralPID takes its polar angle from a linear regression of rho against z, and Brho goes
/// as 1/|sin(polar)|.  Moving the drift velocity to 1.10424 rescales z by 2.9 %, so polar
/// moves, so Brho moves -- and a polygon drawn on the old plane no longer bounds the same
/// tracks.  A gate belongs to the plane it was drawn on.
///
/// Writes a points cache (TTree "pts") so a gate can later be drawn on exactly these numbers
/// rather than on a redrawn approximation, plus a PNG for a first look.
///
///   root -b -q 'pid/pid_plane_dt.C()'                       // every reco in the dv1104 dir
///   root -b -q 'pid/pid_plane_dt.C("run_0016,run_0017")'    // named runs only

void pid_plane_dt(TString runsCSV = "", TString inDir = "/mnt/f/a1975/reco_d2_dv1104/",
                  TString outCache = "pid/pid_plane_dt_dv1104.root",
                  TString outPng = "pid/plots/pid_plane_dt_dv1104.png", double bField = 2.85,
                  double sqrtMax = 60, double brhoMax = 2.0, Long64_t maxEvtPerRun = -1,
                  TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   // ---- which runs -------------------------------------------------------------------
   std::vector<TString> files;
   if (runsCSV.Length() == 0) {
      void *d = gSystem->OpenDirectory(inDir);
      const char *e;
      while ((e = gSystem->GetDirEntry(d))) {
         TString f(e);
         if (f.EndsWith("_multifit_reco.root"))
            files.push_back(inDir + f);
      }
      gSystem->FreeDirectory(d);
      std::sort(files.begin(), files.end());
   } else {
      TObjArray *a = runsCSV.Tokenize(",");
      for (int i = 0; i < a->GetEntries(); ++i)
         files.push_back(inDir + ((TObjString *)a->At(i))->GetString() + "_multifit_reco.root");
   }
   if (files.empty()) {
      printf("no _multifit_reco.root found in %s\n", inDir.Data());
      return;
   }
   printf("=== pid_plane_dt: %zu run(s) from %s ===\n", files.size(), inDir.Data());

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   auto *h = new TH2F("hpid", "a1975 D2 Spyral PID, dv 1.10424;#sqrt{dE/dx};B#rho [T m]", 300, 0, sqrtMax,
                      300, 0, brhoMax);

   TFile *fo = TFile::Open(outCache, "RECREATE");
   auto *tp = new TTree("pts", "Spyral PID points, pattern tracks, dv 1.10424");
   float x, y, polar, vz, vr, ncl, ic;
   int runNo;
   tp->Branch("sqrtdedx", &x);
   tp->Branch("brho", &y);
   tp->Branch("polar", &polar);
   tp->Branch("vertexz", &vz);
   tp->Branch("vertexr", &vr);
   tp->Branch("ncl", &ncl);
   tp->Branch("run", &runNo);
   // ic is STORED, never cut on here: the D2 beam is a cocktail (16C near 1150, a second
   // component near 2060 in comparable numbers), so the beam gate has to stay a knob rather
   // than being baked into the landscape. ic = -1 means no IC file for that run.
   tp->Branch("ic", &ic);

   long nTrk = 0, nVal = 0;
   for (auto &fn : files) {
      TString base = gSystem->BaseName(fn);
      runNo = TString(base(4, 4)).Atoi();
      TFile *fi = TFile::Open(fn);
      if (!fi || fi->IsZombie()) {
         printf("  skip (cannot open) %s\n", base.Data());
         continue;
      }
      TTree *t = (TTree *)fi->Get("cbmsim");
      if (!t) {
         fi->Close();
         continue;
      }
      // --- per-run ion chamber, matched by entry index -------------------------------
      // RETRY: on drvfs a file that is plainly there fails to open INTERMITTENTLY under
      // parallelism. Taking that at face value silently produces an ungated landscape.
      std::vector<float> icv;
      if (icDir.Length()) {
         TString icf = icDir + TString(base(0, 8)) + "_IC.root";
         TFile *fic = nullptr;
         for (int att = 0; att < 5 && !fic; ++att) {
            if (att) gSystem->Sleep(200);
            if (!gSystem->AccessPathName(icf)) fic = TFile::Open(icf);
            if (fic && fic->IsZombie()) { fic->Close(); fic = nullptr; }
         }
         if (fic) {
            TTree *ti = (TTree *)fic->Get("ic");
            float icm = 0;
            ti->SetBranchAddress("icmax", &icm);
            icv.reserve(ti->GetEntries());
            for (Long64_t j = 0; j < ti->GetEntries(); ++j) { ti->GetEntry(j); icv.push_back(icm); }
            fic->Close();
         } else {
            printf("  *** WARNING %s: NO IC after 5 attempts -- every track gets ic = -1 ***\n", base.Data());
         }
      }

      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pe);
      Long64_t N = (maxEvtPerRun > 0) ? std::min(maxEvtPerRun, t->GetEntries()) : t->GetEntries();
      long v0 = nVal;
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (pe->GetEntries() == 0)
            continue;
         auto *p = (AtPatternEvent *)pe->At(0);
         if (!p)
            continue;
         ic = (i < (Long64_t)icv.size()) ? icv[i] : -1.f;
         for (auto &trk : p->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(trk);
            ++nTrk;
            auto r = spy.Estimate(tr);
            if (!r.valid)
               continue;
            ++nVal;
            x = r.sqrtdEdx;
            y = r.brho;
            polar = r.polar * TMath::RadToDeg();
            vz = r.vertex.Z();
            vr = std::sqrt(r.vertex.X() * r.vertex.X() + r.vertex.Y() * r.vertex.Y());
            ncl = r.nClusters;
            h->Fill(x, y);
            tp->Fill();
         }
      }
      printf("  %-34s %8lld evt  %7ld valid\n", base.Data(), N, nVal - v0);
      fi->Close();
   }
   printf("=== pattern tracks %ld, valid Spyral %ld (%.1f%%) ===\n", nTrk, nVal,
          nTrk ? 100.0 * nVal / nTrk : 0.0);

   fo->cd();
   tp->Write();
   h->Write();
   fo->Close();
   printf("cache -> %s\n", outCache.Data());

   auto *c = new TCanvas("cpid", "pid", 1000, 800);
   c->SetLogz();
   h->Draw("colz");
   c->SaveAs(outPng);
   printf("plot  -> %s\n", outPng.Data());
}
