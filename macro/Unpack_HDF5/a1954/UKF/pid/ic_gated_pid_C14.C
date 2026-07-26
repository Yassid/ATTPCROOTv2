/// @file ic_gated_pid_C14.C
/// @brief PID plane (sqrt(dEdx) vs Brho) for a1954 14C, with the IC BEAM GATE applied.
///        14C beam = IC-amplitude peak ~625-750 (NOT the dominant ~1900 contaminant).
///
/// Combines beam-ID (FRIB IC) + recoil-ID (AtSpyralPID on reco tracks). FRIB and reco
/// trees are NOT index-aligned, so events are matched by AtEvent/AtRawEvent GetEventID().
/// Fills the PID plane for IC-gated (14C) vs all events; overlays the 16C+p proton band.
///
///   root -b -q 'ic_gated_pid_C14.C("run_0055,run_0147,run_0148,run_0149,run_0150")'
void ic_gated_pid_C14(TString runsCSV = "run_0055", TString inDir = "/home/yassid/a1954_C14_reco/",
                       double icLo = 625, double icHi = 750, Int_t icTbLo = 1000, Int_t icTbHi = 1350,
                       double bField = 2.85, TString outTag = "",
                       TString refGate =
                          "/home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_HDF5/a1975/UKF/pid/proton_band.json")
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TH2F *hAll = new TH2F("hAll", "PID (all beam);#sqrt{dEdx};B#rho [T m]", 400, 0, 45, 400, 0, 1.0);
   TH2F *hGate = new TH2F("hGate", "PID (14C-beam gated, IC 625-750);#sqrt{dEdx};B#rho [T m]", 400, 0, 45, 400, 0,
                          1.0);
   TH1D *hIC = new TH1D("hIC", "IC amplitude (red = 14C gate);ic_amplitude;events", 500, 0, 2500);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nGateEvt = 0, nAllEvt = 0, nGateTrk = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff = inDir + run + "_FRIB.root", rf = inDir + run + "_reco.root";
      if (gSystem->AccessPathName(ff) || gSystem->AccessPathName(rf)) {
         printf("skip %s (missing FRIB or reco)\n", run.Data());
         continue;
      }
      // 1) build eventID -> IC amplitude from the FRIB tree
      std::map<int, double> icByID;
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

      // 2) loop reco tracks, look up IC by event ID, gate, fill PID
      TFile *fR = TFile::Open(rf);
      TTree *tR = (TTree *)fR->Get("cbmsim");
      TClonesArray *pe = nullptr, *ev = nullptr;
      tR->SetBranchAddress("AtPatternEvent", &pe);
      tR->SetBranchAddress("AtEventH", &ev);
      for (Long64_t i = 0; i < tR->GetEntries(); ++i) {
         tR->GetEntry(i);
         if (ev->GetEntries() == 0)
            continue;
         int id = (int)((AtEvent *)ev->At(0))->GetEventID();
         auto it = icByID.find(id);
         if (it == icByID.end())
            continue;
         double ic = it->second;
         hIC->Fill(ic);
         ++nAllEvt;
         bool inGate = (ic >= icLo && ic <= icHi);
         if (inGate)
            ++nGateEvt;
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
            hAll->Fill(r.sqrtdEdx, r.brho);
            if (inGate) {
               hGate->Fill(r.sqrtdEdx, r.brho);
               ++nGateTrk;
            }
         }
      }
      fR->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nmatched events=%ld  14C-gated (IC %.0f-%.0f)=%ld (%.1f%%)  gated tracks=%ld\n", nAllEvt, icLo, icHi,
          nGateEvt, nAllEvt ? 100.0 * nGateEvt / nAllEvt : 0, nGateTrk);

   TGraph *g = nullptr;
   {
      std::ifstream in(refGate.Data());
      if (in) {
         std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
         auto pos = s.find('[', s.find("vertices"));
         std::vector<double> nums;
         const char *p = s.c_str() + pos, *end = s.c_str() + s.size();
         int depth = 0;
         while (p < end) {
            if (*p == '[')
               depth++;
            if (*p == ']') {
               depth--;
               if (depth <= 0) {
                  ++p;
                  break;
               }
            }
            char *np = nullptr;
            double v = strtod(p, &np);
            if (np != p) {
               nums.push_back(v);
               p = np;
            } else
               ++p;
         }
         g = new TGraph();
         for (size_t i = 0; i + 1 < nums.size(); i += 2)
            g->SetPoint(g->GetN(), nums[i], nums[i + 1]);
         if (g->GetN()) {
            double x0, y0;
            g->GetPoint(0, x0, y0);
            g->SetPoint(g->GetN(), x0, y0);
         }
         g->SetLineColor(kRed);
         g->SetLineWidth(2);
      }
   }

   TCanvas *c = new TCanvas("c", "icgate", 1550, 500);
   c->Divide(3, 1);
   c->cd(1);
   c->cd(1)->SetLogy();
   hIC->Draw();
   {
      TLine *l1 = new TLine(icLo, 1, icLo, hIC->GetMaximum());
      TLine *l2 = new TLine(icHi, 1, icHi, hIC->GetMaximum());
      l1->SetLineColor(kRed);
      l2->SetLineColor(kRed);
      l1->Draw();
      l2->Draw();
   }
   c->cd(2);
   hAll->Draw("colz");
   if (g)
      g->Draw("L same");
   c->cd(3);
   hGate->Draw("colz");
   if (g)
      g->Draw("L same");
   TString png = plotDir + "ic_gated_pid_C14" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
