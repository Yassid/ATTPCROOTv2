/// @file brho_theta_gated_C14.C
/// @brief 14C Bρ vs θ_lab with the IC beam gate (14C) AND the proton PID gate applied.
///        Reads the LOCAL slim cache (AtPatternEvent) + local <run>_FRIB.root (IC).
///        Beam-ID: FRIB IC amplitude in [icLo,icHi]. Recoil-ID: AtSpyralPID + proton polygon.
///        FRIB<->reco matched by GetEventID().
///
///   root -b -q 'pid/brho_theta_gated_C14.C("run_0055,run_0058", "/home/yassid/a1954_C14_reco_hdb_slim/")'

// depth-aware reader of "vertices":[[x,y],...] -> TCutG (proton polygon in sqrtdEdx-brho)
static TCutG *LoadCut(const char *path, const char *name)
{
   std::ifstream in(path);
   if (!in)
      return nullptr;
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto pos = s.find('[', s.find("vertices"));
   if (pos == std::string::npos)
      return nullptr;
   std::vector<double> nums;
   const char *p = s.c_str() + pos, *end = s.c_str() + s.size();
   int depth = 0;
   while (p < end) {
      if (*p == '[')
         depth++;
      if (*p == ']') {
         depth--;
         if (depth <= 0) { ++p; break; }
      }
      char *np = nullptr;
      double v = strtod(p, &np);
      if (np != p) { nums.push_back(v); p = np; }
      else ++p;
   }
   auto *cut = new TCutG(name, nums.size() / 2);
   for (size_t i = 0, k = 0; i + 1 < nums.size(); i += 2, ++k)
      cut->SetPoint(k, nums[i], nums[i + 1]);
   cut->SetLineColor(kRed);
   cut->SetLineWidth(2);
   return cut;
}

void brho_theta_gated_C14(TString runsCSV = "run_0055", TString inDir = "/home/yassid/a1954_C14_reco_hdb_slim/",
                           double icLo = 500, double icHi = 900, Int_t icTbLo = 1050, Int_t icTbHi = 1250,
                           double bField = 2.85, double arclenMin = 0, TString outTag = "",
                           TString protonGate =
                              "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954/UKF/pid/"
                              "proton_14C.json")
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");
   if (!pcut) printf("WARN: no proton gate loaded from %s\n", protonGate.Data());

   TH1D *hIC = new TH1D("hIC", "IC amplitude (red = 14C gate);ic amplitude;events", 500, 0, 2500);
   TH2F *hPID = new TH2F("hPID", "14C PID, IC-gated (red = proton gate);#sqrt{dEdx};B#rho [T m]", 400, 0, 45, 400, 0, 1.2);
   TH2F *hBt = new TH2F("hBt", "14C  B#rho vs #theta_{lab}  (IC + proton gated);#theta_{lab} [deg];B#rho [T m]", 180,
                        0, 180, 300, 0, 1.5);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nMatch = 0, nIC = 0, nTrk = 0, nGtrk = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff = inDir + run + "_FRIB.root", rf = inDir + run + "_slim.root";
      if (gSystem->AccessPathName(ff) || gSystem->AccessPathName(rf))
         continue;
      // 1) eventID -> IC amplitude
      std::map<int, double> icByID;
      TFile *fF = TFile::Open(ff);
      TTree *tF = fF ? (TTree *)fF->Get("cbmsim") : nullptr;
      if (!tF || tF->GetEntries() == 0) { // empty/failed FRIB (e.g. run_0011 unpack crashed)
         printf("SKIP %s (no FRIB tree)\n", run.Data());
         if (fF) fF->Close();
         continue;
      }
      TClonesArray *ra = nullptr;
      tF->SetBranchAddress("AtRawEvent", &ra);
      for (Long64_t i = 0; i < tF->GetEntries(); ++i) {
         tF->GetEntry(i);
         if (ra->GetEntries() == 0) continue;
         auto *raw = (AtRawEvent *)ra->At(0);
         if (!raw || raw->GetGenTraces().empty()) continue;
         auto &adc = raw->GetGenTraces()[0]->GetADC();
         double mx = -1e9;
         for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
            mx = std::max(mx, (double)adc[b]);
         icByID[(int)raw->GetEventID()] = mx;
      }
      fF->Close();
      // 2) slim tracks, IC-gate by eventID, proton-gate, fill
      TFile *fR = TFile::Open(rf);
      TTree *tR = fR ? (TTree *)fR->Get("cbmsim") : nullptr;
      if (!tR || tR->GetEntries() == 0) { if (fR) fR->Close(); continue; }
      TClonesArray *pe = nullptr;
      tR->SetBranchAddress("AtPatternEvent", &pe);
      for (Long64_t i = 0; i < tR->GetEntries(); ++i) {
         tR->GetEntry(i);
         if (pe->GetEntries() == 0) continue;
         auto *p = (AtPatternEvent *)pe->At(0);
         if (!p) continue;
         // AtPatternEvent carries no event ID (PRA drops it) and AtEventH was slimmed out.
         // Reco AtEventH id == entry index and FRIB id == its entry index, both complete
         // sequential streams (verified: id==entry, counts match +/-1), so match by entry i.
         auto it = icByID.find((int)i);
         if (it == icByID.end()) continue;
         ++nMatch;
         hIC->Fill(it->second);
         bool icOK = (it->second >= icLo && it->second <= icHi);
         if (icOK) ++nIC;
         if (!icOK) continue;
         for (auto &trk : p->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(trk);
            auto r = spy.Estimate(tr);
            if (!r.valid || r.arclength < arclenMin) continue;
            ++nTrk;
            hPID->Fill(r.sqrtdEdx, r.brho);
            if (pcut && !pcut->IsInside(r.sqrtdEdx, r.brho)) continue;
            ++nGtrk;
            hBt->Fill(r.polar * TMath::RadToDeg(), r.brho);
         }
      }
      fR->Close();
   }
   printf("\n== gated Brho-theta ==  matched=%ld  IC-14C=%ld (%.1f%%)  IC-gated tracks=%ld  +proton=%ld (%.1f%%)\n",
          nMatch, nIC, nMatch ? 100.0 * nIC / nMatch : 0, nTrk, nGtrk, nTrk ? 100.0 * nGtrk / nTrk : 0);

   TCanvas *c = new TCanvas("c", "gated", 1650, 500);
   c->Divide(3, 1);
   c->cd(1); c->cd(1)->SetLogy(); hIC->Draw();
   { TLine *l1 = new TLine(icLo, 1, icLo, hIC->GetMaximum()); TLine *l2 = new TLine(icHi, 1, icHi, hIC->GetMaximum());
     l1->SetLineColor(kRed); l2->SetLineColor(kRed); l1->Draw(); l2->Draw(); }
   c->cd(2); hPID->Draw("colz"); if (pcut) pcut->Draw("L same");
   c->cd(3); c->cd(3)->SetRightMargin(0.13); hBt->Draw("colz");
   TString png = plotDir + "brho_theta_gated_C14" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
