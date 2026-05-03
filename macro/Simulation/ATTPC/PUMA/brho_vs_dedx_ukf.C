/// @file brho_vs_dedx_ukf.C
/// @brief Brho vs dE/dx PID using AtFitterUKFMulti's chosen hypothesis.
///
/// For each AtFittedTrack:
///   - charge       = ParticleInfo.charge (signed)
///   - mass_amu     = ParticleInfo.mass
///   - KE [MeV]     = Kinematics.kineticEnergy
///   - dE/dx        = TrackProperties.estimateDeDx [a.u./mm]
///   - p [GeV/c]    = sqrt(KE·(KE + 2 m)) with KE,m in GeV
///   - Brho [T·m]   = sign(charge) · p / (0.299792458 · |charge|)
///
/// The sign is fit-driven (the multi-hypothesis fitter picks ±K, ±π by
/// reduced-chi² minimisation, so no cross-product rotation trick is needed).
///
/// Run: root -b -q 'brho_vs_dedx_ukf.C("./data/output_ukf_multi_riemann_1k.root", "./data/brho_dedx_ukf.png")'

namespace {

constexpr double kAmuToGeV = 0.9314940954; // 1 u in GeV/c²
constexpr double kCKMNS = 0.299792458;     // c in T·m / (GeV/c)

struct Species {
   const char *legend;        // legend label (ROOT-rendered)
   const char *fitterName;    // string AtFitterUKFMulti stores in idPDG
   Color_t col;
   Style_t mk;
};

const std::vector<Species> kSpecies = {
   {"K+",   "K+",  kRed + 1, 20},
   {"K-",   "K-",  kBlue + 1, 21},
   {"#pi+", "pi+", kGreen + 2, 22},
   {"#pi-", "pi-", kMagenta + 1, 23},
   {"other", "",   kGray + 2, 24},
};

int speciesIndex(const TString &pdgName)
{
   for (size_t i = 0; i < kSpecies.size() - 1; ++i)
      if (pdgName == kSpecies[i].fitterName)
         return (int)i;
   return (int)kSpecies.size() - 1;
}

} // anonymous namespace

void brho_vs_dedx_ukf(TString ukfFile = "./data/output_ukf_multi_riemann_1k.root",
                      TString outPng = "./data/brho_dedx_ukf.png")
{
   auto *f = TFile::Open(ukfFile);
   if (!f || f->IsZombie()) {
      std::cerr << "Cannot open " << ukfFile << "\n";
      return;
   }
   auto *t = (TTree *)f->Get("cbmsim");
   if (!t) {
      std::cerr << "No cbmsim tree\n";
      return;
   }

   TClonesArray *teArr = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &teArr);

   const Long64_t N = t->GetEntries();

   struct Row {
      double brho;
      double dedx;
      int specIdx;
   };
   std::vector<Row> rows;
   std::array<int, 5> nBySpec{};

   bool diagDone = false;

   for (Long64_t e = 0; e < N; ++e) {
      t->GetEntry(e);
      if (!teArr || teArr->GetEntries() == 0)
         continue;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      const auto &fits = te->GetFittedTracks();
      auto trackArr = te->GetTrackArray(); // hits live here; dE/dx computed below

      for (const auto &ft : fits) {
         const auto &pinfo = ft->GetParticleInfo();
         const auto &kin = ft->GetKinematics();

         int q = pinfo.charge;
         double mass_amu = pinfo.mass;
         double KE = kin.kineticEnergy; // MeV by FairRoot/Genfit convention

         // Lookup matching AtTrack by track-ID for dE/dx (multi-fitter wrapper
         // does not populate TrackProperties.estimateDeDx).
         int tid = ft->GetTrackID();
         double dedx = 0;
         int nHits = 0;
         if (tid >= 0 && tid < (int)trackArr.size()) {
            const auto &hits = trackArr[tid].GetHitArray();
            nHits = hits.size();
            if (nHits >= 3) {
               auto cen = trackArr[tid].GetGeoCenter();
               std::vector<std::tuple<double, double, double, double, double>> sorted;
               sorted.reserve(nHits);
               for (const auto &h : hits) {
                  const auto &p = h->GetPosition();
                  double phi = std::atan2(p.Y() - cen.second, p.X() - cen.first);
                  sorted.emplace_back(phi, p.X(), p.Y(), p.Z(), h->GetCharge());
               }
               std::sort(sorted.begin(), sorted.end());
               double totalQ = 0, trackLen = 0;
               for (size_t i = 0; i < sorted.size(); ++i) {
                  totalQ += std::get<4>(sorted[i]);
                  if (i == 0) continue;
                  double dx = std::get<1>(sorted[i]) - std::get<1>(sorted[i - 1]);
                  double dy = std::get<2>(sorted[i]) - std::get<2>(sorted[i - 1]);
                  double dz = std::get<3>(sorted[i]) - std::get<3>(sorted[i - 1]);
                  trackLen += std::sqrt(dx * dx + dy * dy + dz * dz);
               }
               if (trackLen > 0)
                  dedx = totalQ / trackLen;
            }
         }

         if (!diagDone) {
            std::cout << "[diag] first fitted track: pdg=" << pinfo.idPDG
                      << " q=" << q << " mass=" << mass_amu << " amu  KE[MeV]=" << KE
                      << "  tid=" << tid << "  nHits=" << nHits << "  dedx=" << dedx << "\n";
            diagDone = true;
         }

         if (q == 0 || !(KE > 0) || !(mass_amu > 0) || !(dedx > 0))
            continue;
         // PUMA total CM available energy is ~620 MeV; reject UKF-divergent
         // fits whose KE is far above the kinematic budget.
         if (KE > 1000.0 /* MeV */)
            continue;

         double KE_GeV = KE * 1e-3;
         double m_GeV = mass_amu * kAmuToGeV;
         double p = std::sqrt(KE_GeV * (KE_GeV + 2.0 * m_GeV));
         double brho = (q > 0 ? +1.0 : -1.0) * p / (kCKMNS * std::abs(q));

         int sIdx = speciesIndex(pinfo.idPDG);
         rows.push_back({brho, dedx, sIdx});
         nBySpec[sIdx]++;
      }
   }

   if (rows.empty()) {
      std::cerr << "No fitted tracks with valid KE/charge/mass found.\n";
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
   // Cap axes to PUMA physical kinematic limits, even after the KE filter
   // above (a tracks at e.g. KE=900 MeV is still nonsense kinematically and
   // would stretch the bulk of the distribution into a sliver).
   double brhoMax = std::min(3.0, std::max(0.05, pct(allBrhoAbs, 0.95) * 1.05));
   double dedxMax = std::min(200.0, std::max(1.0, pct(allDedx, 0.95) * 1.05));

   auto *c = new TCanvas("c_brho_dedx_ukf", "Brho vs dE/dx (UKF hypothesis)", 1100, 800);
   c->SetGrid();
   c->SetRightMargin(0.18);
   c->SetLogz();

   const int nbX = 250, nbY = 250;
   auto *hAll = new TH2F("h_brho_dedx_ukf",
                         Form("PUMA UKF PID  -  %zu tracks;sign(q) #upoint Br#rho [T m];dE/dx [a.u. / mm]",
                              rows.size()),
                         nbX, -brhoMax, brhoMax, nbY, 0, dedxMax);
   for (auto &r : rows)
      hAll->Fill(r.brho, r.dedx);
   hAll->SetStats(0);
   hAll->Draw("COLZ");

   std::array<TGraph *, 5> gSpec{};
   for (size_t i = 0; i < kSpecies.size(); ++i) {
      gSpec[i] = new TGraph();
      gSpec[i]->SetMarkerStyle(kSpecies[i].mk);
      gSpec[i]->SetMarkerSize(0.9);
      gSpec[i]->SetMarkerColor(kSpecies[i].col);
      gSpec[i]->SetLineColor(kSpecies[i].col);
   }
   for (auto &r : rows)
      gSpec[r.specIdx]->SetPoint(gSpec[r.specIdx]->GetN(), r.brho, r.dedx);
   for (auto *g : gSpec)
      if (g->GetN() > 0)
         g->Draw("P SAME");

   auto *leg = new TLegend(0.62, 0.66, 0.81, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   for (size_t i = 0; i < kSpecies.size(); ++i)
      leg->AddEntry(gSpec[i], Form("%s  (n=%d)", kSpecies[i].legend, nBySpec[i]), "P");
   leg->Draw();

   c->SaveAs(outPng);
   TString outPdf = outPng;
   outPdf.ReplaceAll(".png", ".pdf");
   c->SaveAs(outPdf);

   std::cout << "\nUKF-chosen hypothesis counts:\n";
   for (size_t i = 0; i < kSpecies.size(); ++i)
      std::cout << "  " << std::setw(6) << kSpecies[i].legend << " : " << nBySpec[i] << "\n";

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
