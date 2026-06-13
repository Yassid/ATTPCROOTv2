/// @file cycleEvents.C
/// @brief INTERACTIVE event-by-event stepper for the AtGenfitter output. Draws the 3
///        projections (XY / ZX / ZY) on a live TCanvas and waits for you at the terminal.
///        2D canvas only — works over X11 (no Eve/OpenGL).
///
/// Run interactively (NOT -b):
///   root -l 'pipeline/cycleEvents.C("/tmp/run_0106_genfitter_gated.root")'
///   root -l 'pipeline/cycleEvents.C("/tmp/run_0106_genfitter_gated.root","pid/deuteron_band.json",false)' // all evts
///
/// Controls (type in the terminal, then Enter):
///   <Enter> next   |   p  previous   |   <number>  jump to that event
///   s  save current as PNG   |   g  toggle gated-only   |   q  quit
///
/// gatedOnly=true (default) steps only events that have a track inside the PID gate.

namespace {
TFile *gF = nullptr;
TTree *gT = nullptr;
TClonesArray *gTrk = nullptr, *gPid = nullptr;
TCanvas *gC = nullptr;
AtTools::AtParticleID gGate;
Double_t gZpad = 1000.0;
}

static void drawOneEvent(Int_t iev)
{
   gT->GetEntry(iev);
   auto *te = (AtTrackingEvent *)gTrk->At(0);
   auto *pe = (gPid && gPid->GetEntriesFast() > 0) ? (AtPIDEvent *)gPid->At(0) : nullptr;
   if (!te) { printf("ev %d: no tracking event\n", iev); return; }
   auto tracks = te->GetTrackArray();
   auto &fits = te->GetFittedTracks();

   std::map<int, bool> inGate;
   if (pe)
      for (auto &r : pe->GetSpyral())
         inGate[r.trackID] = r.valid && gGate.IsInside(r.sqrtdEdx, r.brho);
   int nGated = 0; for (auto &kv : inGate) if (kv.second) nGated++;

   // console table
   printf("\n===== event %d : %zu tracks, %zu fits, %d in gate =====\n", iev, tracks.size(), fits.size(), nGated);
   printf("  %-4s %-6s %-7s %-8s %-8s %-9s %-5s %-5s\n", "trk", "gate", "Rpra", "KE", "theta", "chi2/ndf", "conv", "nhit");
   for (auto &tr : tracks) {
      int id = tr.GetTrackID();
      bool g = inGate.count(id) ? inGate[id] : false;
      const AtFittedTrack *ft = nullptr;
      for (auto &fp : fits) if (fp && fp->GetTrackID() == id) { ft = fp.get(); break; }
      double ke = ft ? ft->GetKinematics(0).kineticEnergy : -1;
      double th = ft ? ft->GetKinematics(0).theta * TMath::RadToDeg() : -1;
      double c2 = -1; int cv = -1;
      if (ft) { auto &m = const_cast<AtFittedTrack *>(ft)->GetTrackMetadata();
                if (m) { c2 = m->GetNdf() > 0 ? m->GetChi2() / m->GetNdf() : -1; cv = m->GetFitConverged(); } }
      printf("  %-4d %-6s %-7.1f %-8.2f %-8.2f %-9.2f %-5d %-5zu\n", id, g ? "DEUT" : "-", tr.GetGeoRadius(), ke, th, c2,
             cv, tr.GetHitArray().size());
   }

   const Color_t pal[] = {kBlue + 1, kGreen + 2, kOrange + 7, kViolet + 1, kCyan + 2, kMagenta + 1, kAzure + 7};
   auto pick = [&](int proj, double X, double Y, double Z, double &a, double &b) {
      if (proj == 0) { a = X; b = Y; } else if (proj == 1) { a = Z; b = X; } else { a = Z; b = Y; }
   };
   const char *axX[3] = {"X [mm]", "z_{lab} [mm]", "z_{lab} [mm]"};
   const char *axY[3] = {"Y [mm]", "X [mm]", "Y [mm]"};

   gC->Clear();
   gC->Divide(3, 1);
   for (int proj = 0; proj < 3; ++proj) {
      gC->cd(proj + 1); gPad->SetGrid();
      double aLo, aHi, bLo, bHi;
      if (proj == 0) { aLo = bLo = -250; aHi = bHi = 250; }
      else { aLo = -50; aHi = gZpad + 50; bLo = -250; bHi = 250; }
      gPad->DrawFrame(aLo, bLo, aHi, bHi, Form("ev %d  %s vs %s;%s;%s", iev, axY[proj], axX[proj], axX[proj], axY[proj]));
      for (size_t i = 0; i < tracks.size(); ++i) {
         auto &tr = tracks[i];
         bool g = inGate.count(tr.GetTrackID()) ? inGate[tr.GetTrackID()] : false;
         auto &hits = tr.GetHitArray();
         if (hits.empty()) continue;
         auto *gr = new TGraph();
         for (auto &h : hits) { const auto &p = h->GetPosition(); double a, b; pick(proj, p.X(), p.Y(), gZpad - p.Z(), a, b); gr->SetPoint(gr->GetN(), a, b); }
         gr->SetMarkerStyle(g ? 20 : 24); gr->SetMarkerSize(g ? 0.55 : 0.4);
         gr->SetMarkerColor(g ? pal[i % 7] : (kGray + 1)); gr->Draw("P SAME");
      }
      for (auto &fp : fits) {
         if (!fp) continue;
         bool g = inGate.count(fp->GetTrackID()) ? inGate[fp->GetTrackID()] : false;
         const auto &kin = fp->GetKinematics(0); const auto &v = fp->GetVertex(0);
         Color_t cc = g ? kRed + 1 : (kGray + 2);
         double a, b; pick(proj, v.X(), v.Y(), v.Z(), a, b);
         auto *m = new TMarker(a, b, 29); m->SetMarkerSize(g ? 2.0 : 1.2); m->SetMarkerColor(cc); m->Draw();
         if (kin.kineticEnergy > 0 && std::isfinite(kin.theta) && std::isfinite(kin.phi)) {
            double L = 60, ux = std::sin(kin.theta) * std::cos(kin.phi), uy = std::sin(kin.theta) * std::sin(kin.phi), uz = std::cos(kin.theta);
            double a1, b1; pick(proj, v.X() + L * ux, v.Y() + L * uy, v.Z() + L * uz, a1, b1);
            auto *l = new TLine(a, b, a1, b1); l->SetLineColor(cc); l->SetLineWidth(g ? 3 : 1); l->SetLineStyle(g ? 1 : 2); l->Draw();
         }
      }
   }
   gC->Modified(); gC->Update();
   gSystem->ProcessEvents();
}

void cycleEvents(TString file = "/tmp/run_0106_genfitter_gated.root", TString gateFile = "pid/deuteron_band.json",
                 Bool_t gatedOnly = true, Double_t zPadPlane = 1000.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gZpad = zPadPlane;
   gGate = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   gF = TFile::Open(file);
   gT = (TTree *)gF->Get("cbmsim");
   if (!gT) { printf("no cbmsim tree in %s\n", file.Data()); return; }
   gT->SetBranchAddress("AtTrackingEvent", &gTrk);
   gT->SetBranchAddress("AtPIDEvent", &gPid);
   Long64_t N = gT->GetEntries();

   // precompute the list of events to step through
   std::vector<Int_t> list;
   for (Long64_t i = 0; i < N; ++i) {
      gT->GetEntry(i);
      auto *te = (AtTrackingEvent *)gTrk->At(0);
      auto *pe = (gPid && gPid->GetEntriesFast() > 0) ? (AtPIDEvent *)gPid->At(0) : nullptr;
      if (!te) continue;
      if (!gatedOnly) { if (!te->GetFittedTracks().empty()) list.push_back(i); continue; }
      bool any = false;
      if (pe) for (auto &r : pe->GetSpyral()) if (r.valid && gGate.IsInside(r.sqrtdEdx, r.brho)) { any = true; break; }
      if (any) list.push_back(i);
   }
   printf("cycleEvents: %zu events to step (%s)\n", list.size(), gatedOnly ? "gated-only" : "all-with-fits");
   if (list.empty()) return;

   gC = new TCanvas("cycle", "cycleEvents", 1800, 600);

   printf("\nControls: <Enter>=next  p=prev  <number>=jump-to-event  s=save  g=toggle-gated  q=quit\n");
   int idx = 0;
   bool gatedNow = gatedOnly;
   drawOneEvent(list[idx]);
   std::string line;
   while (true) {
      printf("[%d/%zu  ev %d] > ", idx + 1, list.size(), list[idx]);
      fflush(stdout);
      if (!std::getline(std::cin, line)) break;
      // trim
      while (!line.empty() && (line.back() == ' ' || line.back() == '\t')) line.pop_back();
      if (line == "q") break;
      else if (line == "p") { idx = (idx - 1 + list.size()) % list.size(); drawOneEvent(list[idx]); }
      else if (line == "s") { TString png = Form("/tmp/cycle_ev%05d.png", list[idx]); gC->SaveAs(png); printf("saved %s\n", png.Data()); }
      else if (line == "g") { printf("(toggle-gated: re-run the macro with gatedOnly=%s to change the step list)\n", gatedNow ? "false" : "true"); }
      else if (!line.empty() && (isdigit(line[0]) || line[0] == '-')) {
         Int_t want = atoi(line.c_str());
         // jump to nearest listed event >= want, else just GetEntry(want) directly
         bool found = false;
         for (size_t k = 0; k < list.size(); ++k) if (list[k] == want) { idx = k; found = true; break; }
         if (found) drawOneEvent(list[idx]);
         else if (want >= 0 && want < N) { printf("(ev %d not in step list; showing it anyway)\n", want); drawOneEvent(want); }
         else printf("out of range 0..%lld\n", N - 1);
      }
      else { idx = (idx + 1) % list.size(); drawOneEvent(list[idx]); } // Enter or anything else = next
   }
   printf("done.\n");
}
