/// @file export_xsec_C14.C
/// @brief Write the 14C(p,p') elastic and inelastic cross sections to a plain text file to send out.
///
/// One file, self-describing, with the elastic in two angular binnings and the inelastic levels in
/// two, so a recipient can pick whichever suits their comparison without asking for a re-export.
/// Everything is absolute mb/sr on the adopted normalisation; nothing here is a shape.
///
/// The header carries the things a number is meaningless without: the beam energy, the vertex slab
/// the yields were counted in, the optical potential the luminosity was measured against and the
/// window it was measured over, and the fact that the quoted errors are statistical only.
///
///   root -b -q 'export_xsec_C14.C()'
void export_xsec_C14(TString outName = "a1954_14C_pp_cross_sections.txt",
                     Double_t lumi = 72.5,
                     // The per-angle yield files, 10 deg and 20 deg binning. These were hardcoded;
                     // they are arguments so the export can follow a different cache without an
                     // edit. The defaults are the adopted production (cat5_s013).
                     TString tag10 = "pr00", TString tag20 = "pr20")
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   std::ofstream o((here + "/" + outName).Data());
   if (!o) { printf("\033[1;31mcannot write %s\033[0m\n", outName.Data()); return; }

   o << "# =====================================================================================\n"
     << "#  14C(p,p') differential cross sections            AT-TPC experiment a1954 (NSCL e15250)\n"
     << "#  generated " << TDatime().AsSQLString() << " by pp/export_xsec_C14.C\n"
     << "# =====================================================================================\n"
     << "#\n"
     << "#  MEASUREMENT\n"
     << "#    reaction            14C(p,p') in inverse kinematics, recoil proton tracked\n"
     << "#    beam                14C at 159.75 MeV (anchored on the resolved 6.094 MeV level,\n"
     << "#                        checked out of sample on the isolated 8.317 MeV 2+)\n"
     << "#    target / fill gas   H2 at 300 torr, room temperature, 3.308e-5 g/cm3\n"
     << "#    runs                14 (0055-0066, 0068, 0069)\n"
     << "#    vertex slab         z = 10-490 mm, applied to BOTH the yields and the acceptance\n"
     << "#\n"
     << "#  RECONSTRUCTION AND FITTING\n"
     << "#    AtPSAMultiFit -> AtDirDeDxCleaner -> AtTrackFinderHDBSCAN\n"
     << "#    GENFIT with material effects ON and CATIMA supplying the stopping power,\n"
     << "#    geometry ATTPC_H300torr_RT, chi2/ndf < 5\n"
     << "#\n"
     << "#  NORMALISATION\n"
     << "#    L = " << lumi << " counts/mb, measured on the elastic channel against the PEREY\n"
     << "#    optical potential over theta_cm 80-120 deg -- the window chosen because L is\n"
     << "#    constant there (chi2/ndf = 1.5 for a constant fit; 5.7-19.1 for the four other\n"
     << "#    global potentials tried, including 15.8 for Koning-Delaroche).\n"
     << "#    dsigma/dOmega = Y / acceptance / dOmega / L,  dOmega = 2pi(cos_lo - cos_hi) exact.\n"
     << "#\n"
     << "#  UNCERTAINTIES\n"
     << "#    The quoted error is STATISTICAL ONLY. Systematics, not included:\n"
     << "#      optical model on the absolute scale   ~10 %\n"
     << "#      optical model on beta                   8 %\n"
     << "#      6.903 MeV blend, on the 7.012 only    -28 %\n"
     << "#    The 6.903 and 7.012 MeV levels are 0.76 sigma apart and are NOT resolved: what is\n"
     << "#    labelled 7.012 is the 6.903+7.012 blend. The 6.091 1- is included for completeness\n"
     << "#    but no collective calculation reproduces its shape.\n"
     << "#\n"
     << "#  COLUMNS   theta_cm [deg]   dsigma/dOmega [mb/sr]   stat. error [mb/sr]\n"
     << "#\n";

   auto dump = [&](const char *title, const char *note, TGraphErrors *g) {
      if (!g) { o << "\n# " << title << " : NOT AVAILABLE\n"; return; }
      o << "\n# -------------------------------------------------------------------------------\n"
        << "# " << title << "\n";
      if (note && *note) o << "#   " << note << "\n";
      o << "# -------------------------------------------------------------------------------\n";
      for (int i = 0; i < g->GetN(); ++i) {
         if (g->GetY()[i] <= 0) continue;
         o << Form("  %7.2f   %12.5e   %12.5e\n", g->GetX()[i], g->GetY()[i], g->GetEY()[i]);
      }
   };

   // ---- elastic, two binnings. The stored graph is Y/A/dOmega; divide by L for mb/sr ----
   const char *elFile[2] = {"plots/elastic_omp_omp.root", "plots/elastic_omp_omp25.root"};
   const char *elName[2] = {"ELASTIC  14C(p,p)  --  5.0 deg binning",
                            "ELASTIC  14C(p,p)  --  2.5 deg binning"};
   for (int k = 0; k < 2; ++k) {
      TFile *f = TFile::Open(here + "/" + elFile[k]);
      auto *g = f && !f->IsZombie() ? (TGraphErrors *)f->Get("elastic_measured") : nullptr;
      if (g) {
         auto *q = new TGraphErrors();
         for (int i = 0; i < g->GetN(); ++i) {
            q->SetPoint(i, g->GetX()[i], g->GetY()[i] / lumi);
            q->SetPointError(i, 0, g->GetEY()[i] / lumi);
         }
         dump(elName[k], "width-tracking window mu +- 2.5w; no background subtracted", q);
      } else dump(elName[k], "", nullptr);
   }

   // ---- inelastic, two binnings, already in mb/sr ----
   const TString inFileS[2] = {"plots/fit_angles_ps_dist_" + tag10 + ".root",
                               "plots/fit_angles_ps_dist_" + tag20 + ".root"};
   const char *inFile[2] = {inFileS[0].Data(), inFileS[1].Data()};
   const char *inBin[2]  = {"10 deg", "20 deg"};
   const char *lv[5]     = {"lvl0", "lvl1", "lvl2", "lvl3", "lvl4"};
   const char *lvName[5] = {"6.091 MeV 1-", "6.728 MeV 3-", "7.012 MeV 2+ (blend with 6.903)",
                            "7.341 MeV 2-", "8.317 MeV 2+"};
   for (int k = 0; k < 2; ++k) {
      TFile *f = TFile::Open(here + "/" + inFile[k]);
      if (!f || f->IsZombie()) continue;
      for (int i = 0; i < 5; ++i)
         dump(Form("INELASTIC  %s  --  %s binning", lvName[i], inBin[k]),
              i == 2 ? "the 6.903 and 7.012 are 0.76 sigma apart and are not resolved"
                     : (i == 0 ? "no collective calculation reproduces this level's shape" : ""),
              (TGraphErrors *)f->Get(lv[i]));
   }
   o << "\n# end\n";
   o.close();
   printf("\n  wrote %s/%s\n", here.Data(), outName.Data());
   gSystem->Exec(Form("wc -l %s/%s", here.Data(), outName.Data()));
}
