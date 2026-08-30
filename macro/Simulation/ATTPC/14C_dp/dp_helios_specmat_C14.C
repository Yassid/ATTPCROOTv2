/// @file dp_helios_specmat_C14.C
/// @brief Two fit-free representations of a long spiral: the HELIOS one and the SpecMAT one.
///
///   root -b -q 'dp_helios_specmat_C14.C("/mnt/f/a1954_C14dp_hits/hits_reco.root","/mnt/f/a1954_C14dp/sims_b285/gs_s8001_sim.root")'
///
/// Both are built from HIT GEOMETRY, with no helix fit anywhere -- which is the point. A backward
/// (d,p) proton spirals for metres and returns to the beam axis; that is exactly the geometry a
/// solenoidal spectrometer exploits, and it is also why a hit-cloud observable can work where a
/// single-momentum helix fit struggles.
///
/// HELIOS. A particle emitted from the axis in a uniform field returns to the axis after one
/// cyclotron period, displaced along z by
///       z_cyc = 2*pi*p_par / (0.3 * B)          [m, GeV/c, T]
/// and the laboratory energy is linear in that displacement:
///       E_lab = E_cm + (1/2) m V_cm^2 + (m V_cm / T_cyc) * z_cyc
/// so states appear as PARALLEL LINES in (z_cyc, E_lab), separated by their Q-value. Here z_cyc is
/// taken from the hits -- the z-advance per turn -- not from a fitted momentum.
///
/// SpecMAT. The largest chord of the hit cloud, d_max, measures the helix diameter and therefore
/// the transverse momentum, with no fit at all. On a1975 the calibration was d_max/2 = 0.93 R
/// because the spiral TIGHTENS as the particle slows, so the largest chord is not the initial
/// diameter -- the same effect that biased the pattern circle here, seen from a different angle.

#include <algorithm>
#include <vector>

static const double U = 931.49401;
static const double MP = 1.007825 * U, M1 = 14.003242 * U, M2 = 2.0141018 * U, M4 = 15.0105993 * U;

static double hs_q(std::vector<double> v, double p)
{
   if (v.size() < 10) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_helios_specmat_C14(TString recoFile = "/mnt/f/a1954_C14dp_hits/hits_reco.root",
                           TString simFile = "/mnt/f/a1954_C14dp/sims_b285/gs_s8001_sim.root",
                           Double_t bField = 2.85, Double_t thMin = 92.0, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);
   const double ZPAD = 1000.0;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   TFile *fr = TFile::Open(recoFile);
   if (!fr || fr->IsZombie()) { printf("\033[1;31mno reco file: %s\033[0m\n", recoFile.Data()); return; }
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);

   auto *hHel = new TH2D("hHel", "A  HELIOS plane, from hit geometry;z per turn [mm];proton KE truth [MeV]", 100, 0,
                         500, 100, 0, 12);
   auto *hDmax = new TH2D("hDmax", "B  SpecMAT: largest chord vs true radius;R_{true} [mm];d_{max}/2 [mm]", 80, 0, 200,
                          80, 0, 200);
   auto *hRatio = new TH1D("hRatio", "C  d_{max} / 2R_{true};d_{max}/(2R_{true});tracks", 60, 0.5, 1.3);
   auto *hZres = new TH1D("hZres", "D  z per turn: hits vs truth;(z_{hits}-z_{true})/z_{true};tracks", 60, -0.5, 0.5);

   std::vector<double> ratio, zres;
   Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
   long n = 0;
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc) continue;
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1) continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pm = std::sqrt(px * px + py * py + pz * pz);
         if (pm <= 0) continue;
         keT = std::sqrt(pm * pm + MP * MP) - MP;
         thT = std::acos(pz / pm) * TMath::RadToDeg();
         break;
      }
      if (keT <= 0 || thT < thMin) continue;
      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : ev->GetTrackCand())
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = const_cast<AtTrack *>(&t); }
      if (!best || nb < 60) continue;

      std::vector<std::array<double, 3>> P;
      for (const auto &h : best->GetHitArray()) {
         auto qp = h->GetPosition();
         P.push_back({qp.X(), qp.Y(), ZPAD - qp.Z()});
      }
      // truth quantities for this proton
      double pTrue = std::sqrt(keT * keT + 2 * keT * MP) / 1000.0;                       // GeV/c
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField); // mm
      double zCycTrue = 2 * TMath::Pi() * 1000.0 * std::fabs(pTrue * std::cos(thT * TMath::DegToRad())) /
                        (0.299792458 * bField);                                          // mm per turn

      // --- SpecMAT observable: the largest chord of the hit cloud, no fit ---
      // O(n^2) is fine at these multiplicities once the cloud is thinned to its convex extremes;
      // here a stride keeps it cheap without changing the maximum appreciably.
      double d2 = 0;
      size_t step = std::max<size_t>(1, P.size() / 300);
      for (size_t a = 0; a < P.size(); a += step)
         for (size_t b = a + step; b < P.size(); b += step) {
            double dx = P[a][0] - P[b][0], dy = P[a][1] - P[b][1];
            double d = dx * dx + dy * dy;
            if (d > d2) d2 = d;
         }
      double dmax = std::sqrt(d2);

      // --- HELIOS observable: the z advance per turn, measured from the hits ---
      // Sort by z (monotonic along a helix), then find the z span over one full turn of azimuth
      // about the cloud's own centre.
      std::sort(P.begin(), P.end(), [](const std::array<double, 3> &a, const std::array<double, 3> &b) {
         return a[2] > b[2];
      });
      double cx = 0, cy = 0;
      for (auto &p : P) { cx += p[0]; cy += p[1]; }
      cx /= P.size(); cy /= P.size();
      double phiPrev = std::atan2(P[0][1] - cy, P[0][0] - cx), unwrapped = 0, zCycHits = -1;
      for (size_t k = 1; k < P.size(); ++k) {
         double ph = std::atan2(P[k][1] - cy, P[k][0] - cx);
         double d = ph - phiPrev;
         while (d > TMath::Pi()) d -= 2 * TMath::Pi();
         while (d < -TMath::Pi()) d += 2 * TMath::Pi();
         unwrapped += d;
         phiPrev = ph;
         if (std::fabs(unwrapped) >= 2 * TMath::Pi()) { zCycHits = std::fabs(P[0][2] - P[k][2]); break; }
      }
      if (rTrue > 0) { hDmax->Fill(rTrue, dmax / 2); ratio.push_back(dmax / (2 * rTrue)); hRatio->Fill(dmax / (2 * rTrue)); }
      if (zCycHits > 0) {
         hHel->Fill(zCycHits, keT);
         if (zCycTrue > 0) { zres.push_back((zCycHits - zCycTrue) / zCycTrue); hZres->Fill((zCycHits - zCycTrue) / zCycTrue); }
      }
      ++n;
   }

   auto *cv = new TCanvas("hs", "hs", 1500, 1000);
   cv->Divide(2, 2);
   cv->cd(1); gPad->SetLogz(); gPad->SetLeftMargin(0.13); hHel->Draw("colz");
   cv->cd(2); gPad->SetLogz(); gPad->SetLeftMargin(0.13); hDmax->Draw("colz");
   { auto *l = new TLine(0, 0, 200, 200); l->SetLineColor(kRed + 1); l->SetLineWidth(2); l->Draw();
     auto *t = new TLatex(0.18, 0.80, "red: d_{max}/2 = R"); t->SetNDC(); t->SetTextSize(0.040);
     t->SetTextColor(kRed + 1); t->Draw(); }
   cv->cd(3); gPad->SetLeftMargin(0.13); hRatio->SetLineWidth(2); hRatio->SetLineColor(kAzure + 2); hRatio->Draw("hist");
   cv->cd(4); gPad->SetLeftMargin(0.13); hZres->SetLineWidth(2); hZres->SetLineColor(kAzure + 2); hZres->Draw("hist");
   cv->SaveAs(outDir + "dp_helios_specmat.png");

   printf("\n  backward tracks used: %ld\n", n);
   printf("  SpecMAT  d_max/(2 R_true): median %.3f, IQR/1.349 %.3f   (a1975 found 0.93 -- the spiral tightens)\n",
          hs_q(ratio, .5), (hs_q(ratio, .75) - hs_q(ratio, .25)) / 1.349);
   printf("  HELIOS   z-per-turn from hits vs truth: median %+.3f, IQR/1.349 %.3f (fractional)\n",
          hs_q(zres, .5), (hs_q(zres, .75) - hs_q(zres, .25)) / 1.349);
   printf("\n  wrote %sdp_helios_specmat.png\n\n", outDir.Data());
   fs->Close(); fr->Close();
}
