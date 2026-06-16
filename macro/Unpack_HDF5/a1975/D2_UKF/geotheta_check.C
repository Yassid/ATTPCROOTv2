/// @file geotheta_check.C
/// @brief Is the PRA GeoTheta (2D-projection, sign=sign(dirX*dirY), AtPRA.cxx:289-315)
///        reliable in forward/backward? Compares it, per PRA track, to an INDEPENDENT
///        3D geometric theta computed from the clusters in the lab frame
///        (z_lab = ZPadPlane - z_digi), with the vertex = cluster nearest the beam axis
///        (the physical reaction vertex; unambiguous for backward tracks whose far end
///        is off-axis). If the two disagree on the >90 (backward) classification often,
///        the 2D GeoTheta is the fragile root cause both fitters inherit.
///
///   root -l -b -q 'geotheta_check.C("run_0016,run_0017,run_0018,run_0019")'

void geotheta_check(TString runsCSV = "run_0016", TString dir = "/mnt/f/a1975/reco_d2/", double zPad = 1000.0)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TH2F *h = new TH2F("h", "3D-geometric #theta vs PRA GeoTheta;PRA GeoTheta [deg];3D geom #theta_{lab} [deg]", 90, 0,
                      180, 90, 0, 180);
   long n = 0, disagreeBwd = 0, praBwd = 0, geoBwd = 0;

   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = dir + run + "_reco.root";
      if (gSystem->AccessPathName(fn))
         continue;
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *pat = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pat);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (pat->GetEntries() == 0)
            continue;
         auto *pe = (AtPatternEvent *)pat->At(0);
         for (auto &tk : pe->GetTrackCand()) {
            double praTh = tk.GetGeoTheta() * TMath::RadToDeg();
            if (praTh <= 0 || praTh > 180)
               continue;
            auto *cl = tk.GetHitClusterArray();
            if (!cl || cl->size() < 3)
               continue;
            // lab-frame positions; vertex = min beam-axis radius, far = max dist from vertex
            int iv = 0;
            double rmin = 1e18;
            for (size_t k = 0; k < cl->size(); ++k) {
               auto p = cl->at(k).GetPosition();
               double r = std::hypot(p.X(), p.Y());
               if (r < rmin) {
                  rmin = r;
                  iv = k;
               }
            }
            auto pv = cl->at(iv).GetPosition();
            double vx = pv.X(), vy = pv.Y(), vz = zPad - pv.Z();
            int ifar = iv;
            double dmax = -1;
            for (size_t k = 0; k < cl->size(); ++k) {
               auto p = cl->at(k).GetPosition();
               double d = std::hypot(std::hypot(p.X() - vx, p.Y() - vy), (zPad - p.Z()) - vz);
               if (d > dmax) {
                  dmax = d;
                  ifar = k;
               }
            }
            auto pf = cl->at(ifar).GetPosition();
            double dx = pf.X() - vx, dy = pf.Y() - vy, dz = (zPad - pf.Z()) - vz;
            double rr = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (rr < 1e-6)
               continue;
            double geoTh = std::acos(dz / rr) * TMath::RadToDeg();
            h->Fill(praTh, geoTh);
            ++n;
            bool pB = praTh >= 90, gB = geoTh >= 90;
            if (pB)
               ++praBwd;
            if (gB)
               ++geoBwd;
            if (pB != gB)
               ++disagreeBwd;
         }
      }
      f->Close();
      printf("processed %s\n", run.Data());
   }

   printf("\n=============== PRA GeoTheta vs independent 3D geometric theta ===============\n");
   printf("  tracks: %ld\n", n);
   printf("  classified BACKWARD (>=90):  PRA GeoTheta = %ld (%.1f%%)   3D-geom = %ld (%.1f%%)\n", praBwd,
          100.0 * praBwd / std::max(1L, n), geoBwd, 100.0 * geoBwd / std::max(1L, n));
   printf("  fwd/bwd DISAGREEMENT: %ld of %ld tracks (%.1f%%)\n", disagreeBwd, n, 100.0 * disagreeBwd / std::max(1L, n));
   printf("==============================================================================\n");
   printf("  (diagonal = agreement; anti-diagonal / off-diagonal = GeoTheta unreliable)\n");

   TCanvas *c = new TCanvas("c", "geotheta_check", 850, 780);
   c->SetRightMargin(0.13);
   c->SetLogz();
   h->Draw("colz");
   TLine *d = new TLine(0, 0, 180, 180);
   d->SetLineColor(kRed);
   d->SetLineStyle(2);
   d->Draw();
   TLine *lx = new TLine(90, 0, 90, 180);
   lx->SetLineColor(kGray + 2);
   lx->SetLineStyle(3);
   lx->Draw();
   TLine *ly = new TLine(0, 90, 180, 90);
   ly->SetLineColor(kGray + 2);
   ly->SetLineStyle(3);
   ly->Draw();
   c->SaveAs("plots/geotheta_check.png");
   printf("saved plots/geotheta_check.png\n");
}
