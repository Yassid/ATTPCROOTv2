/// @file spyral_failcodes.C
/// @brief Why AtSpyralPID rejects the tracks that come back with polar == 0.
///
/// AtSpyralResult::failCode is declared //! -- transient, never streamed -- so the reason a track
/// was rejected does not exist in *_pid.root. It DOES exist at runtime. This re-runs Estimate() on
/// the pattern tracks in <tag>_reco.root, exactly as AtPIDTask::Exec does, and reads the code
/// before it is thrown away.
///
/// WHY polar == 0 SAYS NOTHING ON ITS OWN. polar, radius, direction, nPoints and valid are all
/// assigned in one block at the end of Estimate() (AtSpyralPID.cxx:456-466), so every early exit
/// leaves all of them at their defaults. Only nClusters (line 224) is stamped up front. Reading
/// nPoints == 0 as "too few points" is therefore circular -- it is true of every failure mode.
///
/// The call is kept identical to AtPIDTask: same defaults, field left alone (AtPIDTask's default
/// and AtSpyralPID's are both 2.85, and pidPass_a1975.C passes 2.85 explicitly).
///
///   root -b -q 'spyral_failcodes.C("/mnt/f/a1975_C16_pp_pid","s2001")'

/// @param minPoints  override AtSpyralPID::fMinPoints. <=0 leaves the class default (30), which is
///                   what both consumers use: neither AtPIDTask nor AtGenfitter ever calls
///                   SetMinPoints, so data and simulation are both cut at 30.
void spyral_failcodes(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001", Double_t bField = 2.85,
                      Int_t minPoints = -1)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");

   const char *why[] = {"accepted",
                        "1  too few points (< fMinPoints)",
                        "2  no direction",
                        "3  spline fit failed",
                        "4  circle fit failed",
                        "5  testIndex < 2",
                        "6  Sxx == 0 (degenerate rho-z)",
                        "7  slope == 0",
                        "8  vertex outside beam region",
                        "9  polar/direction inconsistent",
                        "10 no charge accumulated",
                        "11 arclength <= 0"};
   const int NC = 12;
   long code[NC] = {0};
   long preCode = 0; // the two exits at AtSpyralPID.cxx:237/251 that set no code at all
   long nTracks = 0;
   TH1D *hclAll = new TH1D("hclAll", "", 80, 0, 80);
   TH1D *hclBad = new TH1D("hclBad", "", 80, 0, 80);
   TH1D *hnp = new TH1D("hnp", "", 80, 0, 80); // usable points for accepted tracks

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0) {
      spy.SetMinPoints(minPoints);
      printf("  fMinPoints overridden to %d (class default is 30)\n", minPoints);
   }

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fn = dir + "/" + tg + "_reco.root";
      if (gSystem->AccessPathName(fn)) { printf("  skip %s\n", tg.Data()); continue; }
      TFile *f = TFile::Open(fn);
      TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t) { if (f) f->Close(); continue; }
      TClonesArray *pa = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pa);

      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!pa || !pa->GetEntriesFast()) continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;
         for (auto &track : pe->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(track);
            auto res = spy.Estimate(tr);
            ++nTracks;
            hclAll->Fill(res.nClusters);
            if (res.valid) { hnp->Fill(res.nPoints); ++code[0]; }
            else {
               hclBad->Fill(res.nClusters);
               if (res.failCode >= 1 && res.failCode < NC) ++code[res.failCode];
               else ++preCode;
            }
         }
      }
      f->Close();
      printf("  %-8s done\n", tg.Data());
   }
   delete ta;

   printf("\n  pattern tracks processed: %ld\n\n", nTracks);
   printf("  %-38s %8s  %6s\n", "outcome", "tracks", "%");
   printf("  %-38s %8ld  %5.1f\n", why[0], code[0], 100.0 * code[0] / std::max(1L, nTracks));
   for (int c = 1; c < NC; ++c)
      if (code[c])
         printf("  %-38s %8ld  %5.1f\n", why[c], code[c], 100.0 * code[c] / std::max(1L, nTracks));
   if (preCode)
      printf("  %-38s %8ld  %5.1f\n", "(empty hit/cluster array, no code)", preCode,
             100.0 * preCode / std::max(1L, nTracks));

   printf("\n  nClusters: all tracks mean %.1f | rejected mean %.1f\n", hclAll->GetMean(), hclBad->GetMean());
   printf("  usable points on ACCEPTED tracks: mean %.1f, min %.0f  (fMinPoints = 30)\n\n", hnp->GetMean(),
          hnp->GetXaxis()->GetBinLowEdge(hnp->FindFirstBinAbove(0)));
}
