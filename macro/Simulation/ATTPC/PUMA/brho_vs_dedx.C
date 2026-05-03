/// @file brho_vs_dedx.C
/// @brief Brho vs energy-loss (dE/dx) PID plot for PUMA reconstructed tracks,
///        truth-labelled by majority-vote MC PDG.
///
/// For each AtTrack in AtPatternEvent:
///   - Brho [T·m] = B * GetGeoRadius() / 1000
///   - dE/dx [a.u./mm] = sum of hit charges / track 3D arc length
///   - PDG label: each reco hit is matched to the nearest AtMCPoint within
///     fMatchTol_mm; the dominant MC trackID across the track's hits casts
///     the PDG (looked up in MCTrack).
///
/// Run: root -b -q 'brho_vs_dedx.C("./data/output_digi.root", "./data/attpcsim.root", 4.0, "/tmp/brho_dedx.png")'
///
/// The cbmsim trees in the digi and sim files are entry-aligned by construction.

namespace {

constexpr double kMatchTol_mm = 3.0; // hit ↔ MC-point (x,y) match radius
                                     // (3D match suffers an ~8 mm z systematic
                                     // from the drift parametrisation; (x,y) is
                                     // sub-mm and species-discriminating.)

struct Species {
   const char *name;
   int pdg;
   Color_t col;
   Style_t mk;
};

const std::vector<Species> kSpecies = {
   {"K+", 321, kRed + 1, 20},
   {"K-", -321, kBlue + 1, 21},
   {"#pi+", 211, kGreen + 2, 22},
   {"#pi-", -211, kMagenta + 1, 23},
   {"p", 2212, kOrange + 7, 33},
   {"other", 0, kGray + 2, 24},
};

int speciesIndex(int pdg)
{
   for (size_t i = 0; i < kSpecies.size() - 1; ++i)
      if (kSpecies[i].pdg == pdg)
         return (int)i;
   return (int)kSpecies.size() - 1; // "other"
}

} // anonymous namespace

void brho_vs_dedx(TString digiFile = "./data/output_digi.root",
                  TString simFile = "./data/attpcsim.root",
                  double bField_T = 4.0,
                  TString outPng = "/tmp/brho_dedx.png")
{
   auto *fDigi = TFile::Open(digiFile);
   auto *fSim = TFile::Open(simFile);
   if (!fDigi || fDigi->IsZombie() || !fSim || fSim->IsZombie()) {
      std::cerr << "Cannot open input files\n";
      return;
   }
   auto *tDigi = (TTree *)fDigi->Get("cbmsim");
   auto *tSim = (TTree *)fSim->Get("cbmsim");
   if (!tDigi || !tSim) {
      std::cerr << "Missing cbmsim tree\n";
      return;
   }

   TClonesArray *peArr = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack");
   tDigi->SetBranchAddress("AtPatternEvent", &peArr);
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);
   tSim->SetBranchAddress("MCTrack", &mcTrks);

   const Long64_t N = std::min(tDigi->GetEntries(), tSim->GetEntries());

   struct Row {
      double brho;
      double dedx;
      int specIdx;
   };
   std::vector<Row> rows;
   std::array<int, 6> nBySpec{};

   for (Long64_t e = 0; e < N; ++e) {
      tDigi->GetEntry(e);
      if (!peArr || peArr->GetEntries() == 0)
         continue;
      auto *pe = (AtPatternEvent *)peArr->At(0);
      auto &tracks = pe->GetTrackCand();
      if (tracks.empty())
         continue;
      tSim->GetEntry(e);

      // Cache MC point (x,y) (mm) and track IDs for spatial match.
      const int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC);
      std::vector<int> mcTid(nMC);
      for (int k = 0; k < nMC; ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         mcX[k] = mp->GetX() * 10.0; // cm → mm
         mcY[k] = mp->GetY() * 10.0;
         mcTid[k] = mp->GetTrackID();
      }

      for (auto &tr : tracks) {
         double r_mm = tr.GetGeoRadius();
         if (!(r_mm > 0))
            continue;

         const auto &hits = tr.GetHitArray();
         if (hits.size() < 3)
            continue;

         // Charge sign from rotation sense in (x,y). Order hits by distance
         // from vertex (origin) — particles fan out from the PUMA vertex, so
         // r monotonically increases with time. The sign of (p_mid - p_first)
         // × (p_last - p_first) ẑ then encodes B·q.
         int chargeSign = +1;
         {
            std::vector<std::pair<double, std::pair<double, double>>> byR;
            byR.reserve(hits.size());
            for (const auto &h : hits) {
               const auto &p = h->GetPosition();
               byR.emplace_back(p.X() * p.X() + p.Y() * p.Y(),
                                std::make_pair(p.X(), p.Y()));
            }
            std::sort(byR.begin(), byR.end());
            const auto &p0 = byR.front().second;
            const auto &pm = byR[byR.size() / 2].second;
            const auto &pN = byR.back().second;
            double crossZ = (pm.first - p0.first) * (pN.second - p0.second)
                           - (pm.second - p0.second) * (pN.first - p0.first);
            // Convention calibrated against truth labels for PUMA's
            // B = +4 T ẑ in the (z_recon = -z_lab) frame: K+ → crossZ<0.
            chargeSign = (crossZ > 0) ? -1 : +1;
         }
         double brho = chargeSign * bField_T * r_mm / 1000.0; // T·m, signed

         // Phi-sorted arc length and total charge.
         auto cen = tr.GetGeoCenter();
         std::vector<std::tuple<double, double, double, double, double>> sorted; // phi, x, y, z, charge
         sorted.reserve(hits.size());
         for (const auto &h : hits) {
            const auto &p = h->GetPosition();
            double phi = std::atan2(p.Y() - cen.second, p.X() - cen.first);
            sorted.emplace_back(phi, p.X(), p.Y(), p.Z(), h->GetCharge());
         }
         std::sort(sorted.begin(), sorted.end());

         double totalQ = 0;
         double trackLen = 0;
         for (size_t i = 0; i < sorted.size(); ++i) {
            totalQ += std::get<4>(sorted[i]);
            if (i == 0)
               continue;
            const double dx = std::get<1>(sorted[i]) - std::get<1>(sorted[i - 1]);
            const double dy = std::get<2>(sorted[i]) - std::get<2>(sorted[i - 1]);
            const double dz = std::get<3>(sorted[i]) - std::get<3>(sorted[i - 1]);
            trackLen += std::sqrt(dx * dx + dy * dy + dz * dz);
         }
         if (!(trackLen > 0))
            continue;

         // Match each reco hit to nearest MC point in (x,y) within tol;
         // majority-vote the MC trackID, then look up PDG from MCTrack.
         std::map<int, int> tidVotes;
         const double tol2 = kMatchTol_mm * kMatchTol_mm;
         for (const auto &h : hits) {
            const auto &p = h->GetPosition();
            double bestD2 = tol2;
            int bestTid = -1;
            for (int k = 0; k < nMC; ++k) {
               double dx = p.X() - mcX[k];
               double dy = p.Y() - mcY[k];
               double d2 = dx * dx + dy * dy;
               if (d2 < bestD2) {
                  bestD2 = d2;
                  bestTid = mcTid[k];
               }
            }
            if (bestTid >= 0)
               tidVotes[bestTid]++;
         }

         int pdg = 0;
         if (!tidVotes.empty()) {
            int bestTid = -1;
            int bestVotes = 0;
            for (auto &kv : tidVotes) {
               if (kv.second > bestVotes) {
                  bestVotes = kv.second;
                  bestTid = kv.first;
               }
            }
            if (bestTid >= 0 && bestTid < mcTrks->GetEntries()) {
               auto *mt = (AtMCTrack *)mcTrks->At(bestTid);
               pdg = mt->GetPdgCode();
            }
         }

         int sIdx = speciesIndex(pdg);
         rows.push_back({brho, totalQ / trackLen, sIdx});
         nBySpec[sIdx]++;
      }
   }

   if (rows.empty()) {
      std::cerr << "No tracks with valid Brho/length found.\n";
      return;
   }

   auto pct = [](std::vector<double> v, double q) {
      std::sort(v.begin(), v.end());
      return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
   };
   std::vector<double> allBrhoAbs, allDedx;
   for (auto &r : rows) {
      allBrhoAbs.push_back(std::abs(r.brho));
      allDedx.push_back(r.dedx);
   }
   double brhoMax = std::max(0.05, pct(allBrhoAbs, 0.99) * 1.05);
   double dedxMax = std::max(1.0, pct(allDedx, 0.99) * 1.05);

   // ---- Canvas 1: density background + truth-labelled scatter overlay
   auto *c = new TCanvas("c_brho_dedx", "Brho vs dE/dx (truth-labelled)", 1100, 800);
   c->SetGrid();
   c->SetRightMargin(0.18);
   c->SetLogz();

   const int nbX = 250, nbY = 250;
   auto *hAll = new TH2F("h_brho_dedx",
                         Form("PUMA PID  -  %zu tracks  -  B = %.1f T;sign(q) #upoint Br#rho [T m];dE/dx [a.u. / mm]",
                              rows.size(), bField_T),
                         nbX, -brhoMax, brhoMax, nbY, 0, dedxMax);
   for (auto &r : rows)
      hAll->Fill(r.brho, r.dedx);
   hAll->SetStats(0);
   hAll->Draw("COLZ");

   // Truth-labelled per-species scatter overlay.
   std::array<TGraph *, 6> gSpec{};
   for (size_t i = 0; i < kSpecies.size(); ++i) {
      gSpec[i] = new TGraph();
      gSpec[i]->SetMarkerStyle(kSpecies[i].mk);
      gSpec[i]->SetMarkerSize(0.9);
      gSpec[i]->SetMarkerColor(kSpecies[i].col);
      gSpec[i]->SetLineColor(kSpecies[i].col);
   }
   for (auto &r : rows) {
      auto *g = gSpec[r.specIdx];
      g->SetPoint(g->GetN(), r.brho, r.dedx);
   }
   for (auto *g : gSpec)
      if (g->GetN() > 0)
         g->Draw("P SAME");

   auto *leg = new TLegend(0.62, 0.62, 0.81, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   for (size_t i = 0; i < kSpecies.size(); ++i)
      leg->AddEntry(gSpec[i],
                    Form("%s  (n=%d)", kSpecies[i].name, nBySpec[i]),
                    "P");
   leg->Draw();

   c->SaveAs(outPng);
   TString outPdf = outPng;
   outPdf.ReplaceAll(".png", ".pdf");
   c->SaveAs(outPdf);

   // ---- Console summary
   std::cout << "\nTrack count by species:\n";
   for (size_t i = 0; i < kSpecies.size(); ++i)
      std::cout << "  " << std::setw(6) << kSpecies[i].name << " : " << nBySpec[i] << "\n";

   std::vector<double> allBrhoSgn;
   int nPos = 0, nNeg = 0;
   for (auto &r : rows) {
      allBrhoSgn.push_back(r.brho);
      if (r.brho > 0) nPos++; else nNeg++;
   }
   std::cout << "\nsign(q)·Brho [T m]:  n=" << rows.size()
             << "  q>0: " << nPos << "  q<0: " << nNeg
             << "\n  median=" << pct(allBrhoSgn, 0.5)
             << "  p10=" << pct(allBrhoSgn, 0.1)
             << "  p90=" << pct(allBrhoSgn, 0.9) << "\n";
   std::cout << "dE/dx:       median=" << pct(allDedx, 0.5)
             << "  p10=" << pct(allDedx, 0.1)
             << "  p90=" << pct(allDedx, 0.9) << "\n";
   std::cout << "Saved " << outPng << " and " << outPdf << "\n";
}
