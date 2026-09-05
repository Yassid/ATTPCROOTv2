/// @file kin_Ar46.C
/// @brief Build the browser-explorer kinematics cache for one 46Ar(3He,d)47K configuration.
///
///   root -b -q 'kin_Ar46.C("/mnt/f/ar46_3hed_mx_B285_attpc","/mnt/f/ar46_3hed_OLD_2.85T_placeholder","b285_attpc")'
///
/// Writes plots/kin_<tag>.root holding the TNtuple `pk` the explorer generator expects:
///     ke : theta : vertexz : thcm : ex : chi2ndf : state
/// The first six names and their order are the explorer's contract, NOT a local choice --
/// make_explorer_html.C reads them by name and every other channel's cache uses the same six.
/// `state` is an extra trailing column; the generator ignores columns it does not know, and it
/// lets one page separate the three 47K levels that are merged into this file.
///
/// ALL THREE STATES GO IN ONE FILE, because that is what the experiment sees. They were simulated
/// separately (one seed each, so they never share a random sequence), but a spectrum with one
/// level in it answers no question about whether the levels can be told apart. The `state` column
/// keeps the truth available for colouring without pretending the separation is known.
///
/// EVERY NUMBER COMES FROM Ar46::Collect, the same inversion the tables and the figures use.
/// Do not reimplement the kinematics here. In particular the vertex z written out is ALREADY
/// un-mirrored (this simulation reconstructs drift z backwards, r = -1.000 against truth), so a
/// consumer must not apply the mirror a second time. It is written in MILLIMETRES (Collect works
/// in cm; the explorer template is mm) -- see the conversion comment at the Fill.
///
/// THE CACHE STORES THE RECONSTRUCTED QUANTITIES, NOT A FROZEN Ex. `ex` is written for
/// convenience but the page recomputes it from (ke, theta, vertexz) under whatever beam energy
/// the user dials in -- which is the entire point of the viewer, since E_beam here falls
/// 0.957 MeV/cm and the beam arrives at 598 MeV and leaves at ~502.
#include "ex_core_3Hed.h"

void kin_Ar46(TString fitDir = "/mnt/f/ar46_3hed_mx_B285_attpc",
              TString simDir = "/mnt/f/ar46_3hed_OLD_2.85T_placeholder", TString tag = "b285_attpc",
              TString statesCSV = "gs_s3001:0,360_s3011:0.36,2020_s3021:2.02", TString outDir = "plots",
              Double_t dThetaMax = 10.0, Double_t driftLength = 100.0, Bool_t useXtr = kTRUE)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gSystem->mkdir(outDir, kTRUE);

   TString out = outDir + "/kin_" + tag + ".root";
   TFile *f = new TFile(out, "RECREATE");
   TNtuple *pk = new TNtuple("pk", "46Ar(3He,d)47K kinematics", "ke:theta:vertexz:thcm:ex:chi2ndf:state");

   std::unique_ptr<TObjArray> parts(statesCSV.Tokenize(","));
   long total = 0;
   printf("\n%-14s %-14s %8s %8s   %s\n", "state", "Ex_true", "fits", "kept", "cache");
   for (int i = 0; i < parts->GetEntries(); ++i) {
      TString one = ((TObjString *)parts->At(i))->GetString();
      TString tg = one, exs = "0";
      if (one.Contains(":")) {
         tg = one(0, one.Index(":"));
         exs = one(one.Index(":") + 1, one.Length());
      }
      double exTrue = exs.Atof();
      std::unique_ptr<TObjArray> t(tg.Tokenize(","));
      Ar46::Sample S = Ar46::Collect(fitDir, simDir, t.get(), useXtr, dThetaMax, -1.0, driftLength, tag);
      if (!S.ok) {
         printf("%-14s %-14s %8s %8s   NO DATA\n", tg.Data(), exs.Data(), "-", "-");
         continue;
      }
      // theta_cm is NOT stored by Collect and is NOT recomputed here: the page derives it from
      // (ke, theta) with its own kine2b at the beam energy the user has dialled in, so a value
      // frozen at one beam energy would disagree with everything else on the page. -1 = unset.
      // VERTEX z GOES OUT IN MILLIMETRES. Ar46::Collect works in cm because the beam energy
      // loss is quoted in MeV/cm, but the explorer template labels this column "vertex z [mm]",
      // defaults its window to -100..1100 and draws its axis in mm -- and EVERY other channel's
      // cache (ex_Be12.C, ex_C15.C) writes GetVertex().Z() in mm without dividing. Writing cm
      // here does not error: the -100..1100 window accepts 0..100 happily, so the only symptom
      // is a vertex panel squashed into the left 8% of its axis and any mm-based arithmetic on
      // the page silently off by ten.
      for (size_t j = 0; j < S.ex.size(); ++j)
         pk->Fill(S.ke[j], S.theta[j], 10.0 * S.vz[j], -1.0, S.ex[j], S.chi2[j], exTrue);
      total += S.ex.size();
      printf("%-14s %-14.3f %8ld %8zu   mirror %s\n", tg.Data(), exTrue, S.nFit, S.ex.size(),
             S.mirror ? "APPLIED" : "not applied");
   }
   f->cd();
   pk->Write();
   f->Close();
   printf("\nwrote %s  (%ld tracks over %d states)\n", out.Data(), total, parts->GetEntries());
}
