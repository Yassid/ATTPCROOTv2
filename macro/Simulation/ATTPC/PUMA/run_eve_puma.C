/// @file run_eve_puma.C
/// @brief Eve-based interactive 3D event display for PUMA. Shows the annular pad
///        plane, hits + PRA tracks (3D view), the pad pulses, the UKF & GENFIT
///        fitted trajectories, and a per-event table of every fit result.
///
///  *** NEEDS A GRAPHICAL DISPLAY (OpenGL) -- run on the desktop, NOT headless: ***
///        root 'run_eve_puma.C("data/hs_reco375.root")'
///
///  In the sidebar, pick the branches to display:
///    - event branch    -> AtEventH        (PSA hits)
///    - pattern branch   -> AtPatternEvent  (PRA track candidates)
///  The two "Fitted" tabs point at the UKF and GENFIT tracking branches, so
///  their trajectories overlay on the Main 3D view. The "Fit results" tab
///  tabulates, for the current event, every fitted track (curve) of both
///  fitters -- charge, momentum, theta, vertex-z, chi2/ndf, and (when a sim
///  file is given) the generator |p| and residual dp/p, colour-coded by quality
///  (green <5%, orange <15%, red otherwise). Use the < / > navigator to step.
///
///  Truth: pass the entry-aligned sim file as the 2nd argument (default
///  data/hs_sim375.root); "" disables the truth columns.
///
/// NB: like the rest of this macro, the helpers below use NO #includes -- the
/// AT-TPC / ROOT / FairRoot types are resolved by ROOT's rootmap autoloading
/// (needs `source build/config.sh`). Adding explicit #includes pulls the AtData
/// headers off the interpreter path and breaks dictionary parsing.

/// Rest mass [MeV/c^2] implied by the fitted PDG label (PUMA species).
static double PumaFitMass(const TString &pdg)
{
   if (pdg.Contains("pi") || pdg.Contains("Pi"))
      return 139.57039;
   if (pdg.Contains("K") || pdg.Contains("kaon"))
      return 493.677;
   if (pdg.Contains("prot") || pdg == "2212")
      return 938.272;
   return 139.57039; // default: pion
}

/// One generator (truth) primary meson for the current event.
struct PumaTruth {
   int sign;     // charge sign (+1 / -1)
   double p;     // generator |p| [MeV/c]
   double theta; // generator theta [deg]
   bool used;    // already matched to a fitted track this block
};

/// Read the primary pi/K generator tracks for one entry of the (entry-aligned)
/// simulation tree. Momentum in AtMCTrack is GeV/c -> x1000 for MeV/c. Empty if
/// no sim tree was provided.
static std::vector<PumaTruth> PumaReadTruth(TTree *simTree, TClonesArray *mcTrks, long entry)
{
   std::vector<PumaTruth> out;
   if (simTree == nullptr || mcTrks == nullptr || entry < 0 || entry >= simTree->GetEntries())
      return out;
   simTree->GetEntry(entry);
   for (int i = 0; i < mcTrks->GetEntries(); ++i) {
      auto *m = static_cast<AtMCTrack *>(mcTrks->At(i));
      int pdg = m->GetPdgCode();
      if (m->GetMotherId() != -1)
         continue; // primaries only
      if (std::abs(pdg) != 211 && std::abs(pdg) != 321)
         continue; // pi+-, K+- : sign(pdg) == sign(charge)
      double px = m->GetPx(), py = m->GetPy(), pz = m->GetPz();
      double p = 1000.0 * std::sqrt(px * px + py * py + pz * pz);
      double th = (p > 0) ? std::acos(1000.0 * pz / p) * 180.0 / M_PI : 0.0;
      out.push_back({pdg > 0 ? +1 : -1, p, th, false});
   }
   return out;
}

/// Colour for a |dp/p| residual: green good, orange marginal, red poor.
static Color_t PumaResidualColor(double dppAbs)
{
   if (dppAbs < 5)
      return kGreen + 2;
   if (dppAbs < 15)
      return kOrange + 7;
   return kRed + 1;
}

/// Draw-event function for the "Fit results" AtTabMacro tab. Renders, for the
/// current event, one row per fitted track (curve) of BOTH fitters, read
/// straight from FairRootManager (same self-contained approach as AtTabFitted).
/// Columns: track, charge, momentum, theta, vertex-z, chi2/ndf; and -- when the
/// entry-aligned sim tree is supplied -- generator |p| and the residual dp/p
/// (colour-coded) from matching each track to a same-charge primary. simTree /
/// mcTrks may be null, in which case the truth columns are omitted.
static void DrawPumaFitTable(AtTabInfo * /*unused*/, TTree *simTree, TClonesArray *mcTrks)
{
   gPad->Clear();
   long entry = AtViewerManager::Instance()->GetCurrentEntry().Get();
   auto truthAll = PumaReadTruth(simTree, mcTrks, entry);
   const bool hasTruth = !truthAll.empty();

   // Fixed NDC column x-positions (aligned across both fitter blocks).
   const double xTrk = 0.03, xQ = 0.13, xP = 0.19, xTh = 0.35, xVz = 0.47, xChi = 0.59, xTru = 0.73, xDp = 0.88;

   TLatex tl;
   tl.SetNDC();
   tl.SetTextFont(42);

   double y = 0.95;
   tl.SetTextSize(0.05);
   tl.SetTextColor(kBlack);
   tl.DrawLatex(xTrk, y, TString::Format("Fit results   -   Entry %ld", entry));
   y -= 0.085;

   auto drawBlock = [&](const char *branch, Color_t c, const char *label, std::vector<PumaTruth> truth) {
      tl.SetTextColor(c);
      tl.SetTextSize(0.038);
      tl.DrawLatex(xTrk, y, label);
      y -= 0.052;

      tl.SetTextColor(kGray + 3);
      tl.SetTextSize(0.03);
      tl.DrawLatex(xTrk, y, "trk");
      tl.DrawLatex(xQ, y, "q");
      tl.DrawLatex(xP, y, "p [MeV/c]");
      tl.DrawLatex(xTh, y, "#theta [#circ]");
      tl.DrawLatex(xVz, y, "v_{z} [mm]");
      tl.DrawLatex(xChi, y, "#chi^{2}/ndf");
      if (hasTruth) {
         tl.DrawLatex(xTru, y, "p_{gen}");
         tl.DrawLatex(xDp, y, "#Deltap/p");
      }
      y -= 0.048;

      auto *rm = FairRootManager::Instance();
      auto *arr = dynamic_cast<TClonesArray *>(rm->GetObject(branch));
      const AtTrackingEvent *te = (arr && arr->GetEntries() > 0) ? static_cast<AtTrackingEvent *>(arr->At(0)) : nullptr;
      if (te == nullptr || te->GetFittedTracks().empty()) {
         tl.SetTextColor(kGray + 1);
         tl.SetTextSize(0.03);
         tl.DrawLatex(xTrk, y, "(no fitted tracks)");
         y -= 0.06;
         return;
      }

      tl.SetTextSize(0.03);
      for (const auto &ft : te->GetFittedTracks()) {
         const auto &kin = ft->GetKinematics(0);
         const auto &pinfo = ft->GetParticleInfo(0);
         double KE = kin.kineticEnergy;
         double p = (KE > 0) ? std::sqrt(KE * KE + 2 * KE * PumaFitMass(pinfo.idPDG)) : -1;

         int sign = 0;
         if (pinfo.charge > 0 || pinfo.idPDG.Contains("+"))
            sign = +1;
         else if (pinfo.charge < 0 || pinfo.idPDG.Contains("-"))
            sign = -1;

         const auto &props = ft->GetTrackPropertiesStruct();
         double vz = props.initialPositionXtr.Z();
         if (std::abs(props.initialPositionXtr.X()) < 1e-9 && std::abs(props.initialPositionXtr.Y()) < 1e-9 &&
             std::abs(vz) < 1e-9)
            vz = ft->GetVertex(0).Z();

         TString chi = "--";
         auto &meta = ft->GetTrackMetadata();
         if (meta && meta->GetNdf() > 0)
            chi = TString::Format("%.1f/%d", meta->GetChi2(), meta->GetNdf());

         tl.SetTextColor(kBlack);
         tl.DrawLatex(xTrk, y, TString::Format("%d", ft->GetTrackID()));
         tl.DrawLatex(xQ, y, sign > 0 ? "+" : (sign < 0 ? "#minus" : "0"));
         tl.DrawLatex(xP, y, p > 0 ? TString::Format("%.1f", p) : TString("--"));
         tl.DrawLatex(xTh, y, TString::Format("%.1f", kin.theta * 180.0 / M_PI));
         tl.DrawLatex(xVz, y, TString::Format("%.1f", vz));
         tl.DrawLatex(xChi, y, chi);

         // match to a same-charge, still-unused truth primary; report dp/p
         for (auto &t : truth)
            if (!t.used && t.sign == sign && sign != 0) {
               t.used = true;
               double dpp = (t.p > 0) ? 100.0 * (p - t.p) / t.p : 0.0;
               tl.SetTextColor(kGray + 3);
               tl.DrawLatex(xTru, y, TString::Format("%.1f", t.p));
               tl.SetTextColor(PumaResidualColor(std::abs(dpp)));
               tl.DrawLatex(xDp, y, TString::Format("%+.1f%%", dpp));
               break;
            }
         y -= 0.045;
      }
      y -= 0.025;
   };

   drawBlock("AtTrackingEventUKF", kAzure + 2, "UKF   (AtTrackingEventUKF)", truthAll);
   drawBlock("AtTrackingEventGenfit", kRed + 1, "GENFIT   (AtTrackingEventGenfit)", truthAll);
}

/// @param inputFile reco file (hits + PRA + both fitters' tracking events)
/// @param simFile   entry-aligned simulation file with MCTrack; enables the
///                  truth |p| and dp/p columns in the 3D overlay. Pass "" to
///                  disable the truth comparison.
void run_eve_puma(TString inputFile = "data/hs_reco375.root", TString simFile = "data/hs_sim375.root")
{
   TString dir = getenv("VMCWORKDIR");
   TString geoFile = dir + "/geometry/ATTPC_PUMA_geomanager.root";

   auto *fRun = new FairRunAna();
   fRun->SetSource(new FairFileSource(inputFile));
   fRun->SetSink(new FairRootFileSink("eve_puma_out.root"));
   fRun->SetGeomFile(geoFile);

   auto *rtdb = fRun->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open((dir + "/parameters/ATTPC.PUMA_sim.par").Data(), "in");
   rtdb->setFirstInput(parIo);

   // PUMA annular pad plane: 16 equal-area rings x 256 azimuthal pads (R = 62.9-121.1 mm)
   auto fMap = std::make_shared<AtTpcPUMAMap>(62.9, 121.1, 16, 256);
   fMap->GeneratePadPlane();

   auto *eveMan = new AtViewerManager(fMap);

   // Optional generator truth: open the (entry-aligned) sim file so the "Fit
   // results" tab can show p_gen and the residual dp/p per fitted track. Opened
   // before the tabs so the draw function can capture the tree.
   TTree *simTree = nullptr;
   TClonesArray *mcTrks = nullptr;
   if (simFile.Length() > 0) {
      auto *fSim = TFile::Open(simFile);
      if (fSim && !fSim->IsZombie()) {
         simTree = dynamic_cast<TTree *>(fSim->Get("cbmsim"));
         if (simTree) {
            mcTrks = new TClonesArray("AtMCTrack");
            simTree->SetBranchAddress("MCTrack", &mcTrks);
            std::cout << "Truth comparison enabled from " << simFile << " (" << simTree->GetEntries() << " entries)."
                      << std::endl;
         }
      }
      if (simTree == nullptr)
         std::cout << "WARNING: could not open truth sim file '" << simFile
                   << "'; Fit results tab shows fit values only." << std::endl;
   }

   // (1) 3D view: hits + PRA track candidates
   auto tabMain = std::make_unique<AtTabMain>();
   tabMain->SetMultiHit(100);
   eveMan->AddTab(std::move(tabMain));

   // (2) UKF fitted trajectories (blue), smoothed polyline, 4 T
   auto tabUKF = std::make_unique<AtTabFitted>("UKF fit", "AtTrackingEventUKF");
   tabUKF->SetDrawSmoothed(true);
   tabUKF->SetBField(4.0);
   tabUKF->SetTrackColor(kAzure + 2);
   eveMan->AddTab(std::move(tabUKF));

   // (3) GENFIT fitted trajectories (red)
   auto tabGF = std::make_unique<AtTabFitted>("GENFIT fit", "AtTrackingEventGenfit");
   tabGF->SetDrawSmoothed(true);
   tabGF->SetBField(4.0);
   tabGF->SetTrackColor(kRed + 1);
   eveMan->AddTab(std::move(tabGF));

   // (4) pad pulses: raw ADC + PSA-processed ADC
   auto tabPad = std::make_unique<AtTabPad>(2, 1, "Pulses");
   tabPad->DrawRawADC(0, 0);
   tabPad->DrawADC(0, 1);
   eveMan->AddTab(std::move(tabPad));

   // (5) per-event table of every fit result (both fitters, one row per curve),
   // with the truth p_gen and residual dp/p when a sim file was opened.
   auto tabFitRes = std::make_unique<AtTabMacro>(1, 1, "Fit results");
   tabFitRes->SetDrawEventFunction(
      [simTree, mcTrks](AtTabInfo *info) { DrawPumaFitTable(info, simTree, mcTrks); }, 0, 0);
   eveMan->AddTab(std::move(tabFitRes));

   eveMan->Init();

   std::cout << "PUMA event display ready. Select AtEventH / AtPatternEvent in the sidebar." << std::endl;
   std::cout << "Open the 'Fit results' tab for per-track fit values + truth dp/p." << std::endl;
}
