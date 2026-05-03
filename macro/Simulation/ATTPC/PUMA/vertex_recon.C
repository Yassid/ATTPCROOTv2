/// @file vertex_recon.C
/// @brief Reconstructed-vertex resolution plot for PUMA UKF fits.
///
/// For each AtFittedTrack uses TrackProperties.initialPositionXtr — the
/// closest point to (0,0) along the back-extrapolated pattern — and compares
/// to the MC vertex (MCTrack[0] start vertex).
///
/// Conventions:
///   - MCTrack vertex is in cm (FairMCTrack convention) — convert to mm.
///   - Drift convention: hits in the reco frame have z_recon = -z_lab. Flip
///     MC z when comparing.
///   - Truth species via the same (x,y) hit-to-AtMCPoint match used in
///     pid_confusion / brho_vs_dedx.
///
/// Output: `data/vertex_recon.png` + `.pdf` with 4 panels (Δx, Δy, Δz, 3D r),
///         coloured by truth species.
///
/// Run: root -b -q 'vertex_recon.C("./data/output_ukf_multi.root", "./data/attpcsim.root", "./data/vertex_recon.png")'

namespace {

constexpr double kMatchTol_mm = 3.0;

struct Species {
   const char *legend;
   int pdg;
   Color_t col;
};

const std::vector<Species> kSpecies = {
   {"K+", 321, kRed + 1},
   {"#pi-", -211, kMagenta + 1},
   {"other", 0, kGray + 2},
};

int speciesIdx(int pdg)
{
   for (size_t i = 0; i < kSpecies.size() - 1; ++i)
      if (kSpecies[i].pdg == pdg) return (int)i;
   return (int)kSpecies.size() - 1;
}

} // namespace

void vertex_recon(TString ukfFile = "./data/output_ukf_multi.root",
                  TString simFile = "./data/attpcsim.root",
                  TString outPng = "./data/vertex_recon.png")
{
   auto *fU = TFile::Open(ukfFile);
   auto *fS = TFile::Open(simFile);
   if (!fU || fU->IsZombie() || !fS || fS->IsZombie()) {
      std::cerr << "Cannot open inputs\n";
      return;
   }
   auto *tU = (TTree *)fU->Get("cbmsim");
   auto *tS = (TTree *)fS->Get("cbmsim");
   if (!tU || !tS) {
      std::cerr << "Missing cbmsim trees\n";
      return;
   }

   TClonesArray *teArr = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tU->SetBranchAddress("AtTrackingEvent", &teArr);
   tS->SetBranchAddress("AtTpcPoint", &mcPts);
   tS->SetBranchAddress("MCTrack", &mcTrks);

   const Long64_t N = std::min(tU->GetEntries(), tS->GetEntries());
   const int nSp = (int)kSpecies.size();

   struct Row { double dx, dy, dz, r; int sIdx; };
   std::vector<Row> rows;

   // For per-event averaged reconstructed vertex:
   struct EvRow { double dx, dy, dz, r; int nTracks; };
   std::vector<EvRow> evRows;

   bool diagDone = false;

   for (Long64_t e = 0; e < N; ++e) {
      tU->GetEntry(e);
      if (!teArr || teArr->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      const auto &fits = te->GetFittedTracks();
      if (fits.empty()) continue;
      auto trackArr = te->GetTrackArray();
      tS->GetEntry(e);

      if (mcTrks->GetEntries() == 0) continue;
      auto *mt0 = (AtMCTrack *)mcTrks->At(0);
      // MC vertex (mm). AtFittedTrack stores `initialPositionXtr` in the lab
      // frame (positive z), same convention as MCTrack — no z-flip needed.
      double vx_mc = mt0->GetStartX() * 10.0;
      double vy_mc = mt0->GetStartY() * 10.0;
      double vz_mc = mt0->GetStartZ() * 10.0;

      // Cache MC points for truth labelling
      const int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC);
      std::vector<int> mcTid(nMC);
      for (int k = 0; k < nMC; ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         mcX[k] = mp->GetX() * 10.0;
         mcY[k] = mp->GetY() * 10.0;
         mcTid[k] = mp->GetTrackID();
      }

      double sumX = 0, sumY = 0, sumZ = 0;
      int nGood = 0;

      for (const auto &ft : fits) {
         const auto &props = ft->GetTrackPropertiesStruct();
         double xR = props.initialPositionXtr.X();
         double yR = props.initialPositionXtr.Y();
         double zR = props.initialPositionXtr.Z();

         // Skip if back-extrapolation didn't run (default-zero or absurd).
         if (std::abs(xR) < 1e-12 && std::abs(yR) < 1e-12 && std::abs(zR) < 1e-12)
            continue;

         double dx = xR - vx_mc;
         double dy = yR - vy_mc;
         double dz = zR - vz_mc;
         double r = std::sqrt(dx * dx + dy * dy + dz * dz);

         // Truth species via (x,y) hit match
         int truthPdg = 0;
         int tid = ft->GetTrackID();
         if (tid >= 0 && tid < (int)trackArr.size()) {
            const auto &hits = trackArr[tid].GetHitArray();
            std::map<int, int> tidVotes;
            const double tol2 = kMatchTol_mm * kMatchTol_mm;
            for (const auto &h : hits) {
               const auto &p = h->GetPosition();
               double bestD2 = tol2;
               int bestTid = -1;
               for (int k = 0; k < nMC; ++k) {
                  double ex = p.X() - mcX[k];
                  double ey = p.Y() - mcY[k];
                  double d2 = ex * ex + ey * ey;
                  if (d2 < bestD2) { bestD2 = d2; bestTid = mcTid[k]; }
               }
               if (bestTid >= 0) tidVotes[bestTid]++;
            }
            if (!tidVotes.empty()) {
               int bestTid = -1, bestVotes = 0;
               for (auto &kv : tidVotes)
                  if (kv.second > bestVotes) { bestVotes = kv.second; bestTid = kv.first; }
               if (bestTid >= 0 && bestTid < mcTrks->GetEntries())
                  truthPdg = ((AtMCTrack *)mcTrks->At(bestTid))->GetPdgCode();
            }
         }

         if (!diagDone) {
            std::cout << "[diag] event " << e
                      << "  MC vtx (reco frame) = (" << vx_mc << ", " << vy_mc << ", " << vz_mc << ") mm\n"
                      << "       fit Xtr        = (" << xR << ", " << yR << ", " << zR << ") mm\n"
                      << "       Δ = (" << dx << ", " << dy << ", " << dz << ")  r=" << r << " mm\n";
            diagDone = true;
         }

         rows.push_back({dx, dy, dz, r, speciesIdx(truthPdg)});
         sumX += xR; sumY += yR; sumZ += zR;
         nGood++;
      }

      if (nGood > 0) {
         double avX = sumX / nGood, avY = sumY / nGood, avZ = sumZ / nGood;
         double dx = avX - vx_mc, dy = avY - vy_mc, dz = avZ - vz_mc;
         evRows.push_back({dx, dy, dz, std::sqrt(dx * dx + dy * dy + dz * dz), nGood});
      }
   }

   if (rows.empty()) {
      std::cerr << "No usable fitted tracks with back-extrapolation.\n";
      return;
   }

   // ---- Resolution histograms split by truth species
   // Auto-range from the data 99th percentile (each axis independently).
   auto pctAbs = [](std::vector<double> v, double q) {
      for (auto &x : v) x = std::abs(x);
      std::sort(v.begin(), v.end());
      return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
   };
   std::vector<double> allX, allY, allZ, allR;
   for (auto &row : rows) { allX.push_back(row.dx); allY.push_back(row.dy);
                            allZ.push_back(row.dz); allR.push_back(row.r); }
   double xmax = std::max(5.0, pctAbs(allX, 0.95) * 1.1);
   double ymax = std::max(5.0, pctAbs(allY, 0.95) * 1.1);
   double zmax = std::max(10.0, pctAbs(allZ, 0.95) * 1.1);
   double rmax = std::max(10.0, pctAbs(allR, 0.95) * 1.1);

   const int nb = 80;
   std::array<TH1F *, 3> hX{}, hY{}, hZ{}, hR{};
   for (int s = 0; s < nSp; ++s) {
      hX[s] = new TH1F(Form("hX_%d", s), Form(";#Deltax = x_{xtr} - x_{MC} [mm];tracks"), nb, -xmax, xmax);
      hY[s] = new TH1F(Form("hY_%d", s), Form(";#Deltay = y_{xtr} - y_{MC} [mm];tracks"), nb, -ymax, ymax);
      hZ[s] = new TH1F(Form("hZ_%d", s), Form(";#Deltaz = z_{xtr} - z_{MC} [mm];tracks"), nb, -zmax, zmax);
      hR[s] = new TH1F(Form("hR_%d", s), Form(";|#Deltar| [mm];tracks"), nb, 0, rmax);
      for (auto *h : {hX[s], hY[s], hZ[s], hR[s]}) {
         h->SetLineColor(kSpecies[s].col);
         h->SetLineWidth(2);
         h->SetStats(0);
      }
   }
   for (auto &row : rows) {
      hX[row.sIdx]->Fill(row.dx);
      hY[row.sIdx]->Fill(row.dy);
      hZ[row.sIdx]->Fill(row.dz);
      hR[row.sIdx]->Fill(row.r);
   }

   auto *c = new TCanvas("c_vtx", "PUMA reco vertex - track-level", 1300, 900);
   c->Divide(2, 2);

   auto drawStack = [&](int pad, std::array<TH1F *, 3> &hs, const char *title) {
      c->cd(pad);
      gPad->SetGrid();
      double maxY = 0;
      for (auto *h : hs) maxY = std::max(maxY, h->GetMaximum());
      bool first = true;
      for (auto *h : hs) {
         h->SetMaximum(maxY * 1.15);
         h->SetTitle(title);
         h->Draw(first ? "HIST" : "HIST SAME");
         first = false;
      }
      auto *leg = new TLegend(0.62, 0.68, 0.88, 0.88);
      leg->SetBorderSize(0); leg->SetFillStyle(0);
      for (size_t s = 0; s < kSpecies.size(); ++s)
         leg->AddEntry(hs[s], Form("%s  (n=%lld)", kSpecies[s].legend, (Long64_t)hs[s]->GetEntries()), "L");
      leg->Draw();
   };
   drawStack(1, hX, "#Deltax (track xtr - MC vertex)");
   drawStack(2, hY, "#Deltay");
   drawStack(3, hZ, "#Deltaz");
   drawStack(4, hR, "3D distance |#Deltar|");

   c->SaveAs(outPng);
   TString outPdf = outPng; outPdf.ReplaceAll(".png", ".pdf");
   c->SaveAs(outPdf);

   // ---- Resolution summary
   auto sigmaIQR = [](std::vector<double> v) {
      std::sort(v.begin(), v.end());
      if (v.size() < 4) return 0.0;
      double q25 = v[v.size() / 4];
      double q75 = v[(3 * v.size()) / 4];
      return (q75 - q25) / 1.349; // robust σ for Gaussian-like
   };
   auto med = [](std::vector<double> v) {
      std::sort(v.begin(), v.end());
      return v.empty() ? 0.0 : v[v.size() / 2];
   };

   std::cout << "\nPer-track xtr-vs-MC vertex residuals (n=" << rows.size() << "):\n";
   std::cout << "  median Δx = " << med(allX) << " mm,  σ_IQR = " << sigmaIQR(allX) << " mm\n";
   std::cout << "  median Δy = " << med(allY) << " mm,  σ_IQR = " << sigmaIQR(allY) << " mm\n";
   std::cout << "  median Δz = " << med(allZ) << " mm,  σ_IQR = " << sigmaIQR(allZ) << " mm\n";
   std::cout << "  median |Δr| = " << med(allR) << " mm\n";

   if (!evRows.empty()) {
      std::vector<double> ex, ey, ez, er;
      for (auto &v : evRows) { ex.push_back(v.dx); ey.push_back(v.dy);
                               ez.push_back(v.dz); er.push_back(v.r); }
      std::cout << "\nPer-event averaged reconstructed vertex (n=" << evRows.size() << " events):\n";
      std::cout << "  median Δx = " << med(ex) << "  σ_IQR = " << sigmaIQR(ex) << " mm\n";
      std::cout << "  median Δy = " << med(ey) << "  σ_IQR = " << sigmaIQR(ey) << " mm\n";
      std::cout << "  median Δz = " << med(ez) << "  σ_IQR = " << sigmaIQR(ez) << " mm\n";
      std::cout << "  median |Δr| = " << med(er) << " mm\n";
   }

   std::cout << "\nSaved " << outPng << " and " << outPdf << "\n";
}
