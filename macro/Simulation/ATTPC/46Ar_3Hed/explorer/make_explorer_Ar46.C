/// @file make_explorer_Ar46.C
/// @brief Browser explorer for the 46Ar(3He,d)47K simulation: N configurations, one page.
///
///   root -b -q 'explorer/make_explorer_Ar46.C()'
///
/// THE PAGE ITSELF IS NOT WRITTEN HERE. It is a1954_Be12/UKF/pp/explorer_template.html, the same
/// template every other channel uses, reused verbatim by path. Building a parallel viewer has
/// already cost one session; do not start a second one. This file only supplies the data and the
/// config block the template asks for.
///
/// WHY A 46Ar GENERATOR RATHER THAN make_explorer_html.C. That one is capped at THREE sets
/// (TString sn[3]), and this campaign has four configurations -- 2.85/3.9 T x AT-TPC/2 mm. The
/// TEMPLATE has no such cap: it builds one button per key present in the DATA object
/// ("Build one button per set actually present, not a fixed ukf/genfit pair"), so N sets need
/// only a generator that emits N. Everything else here matches make_explorer_html.C's output
/// byte for byte in shape -- same JSON keys, same CFG fields, same two placeholders.
///
/// ONE EXTRA COLUMN: `state`, the TRUE excitation energy of the level a track came from. The
/// three 47K levels are simulated separately and merged into one cache by kin_Ar46.C, because a
/// spectrum with one level in it cannot answer whether the levels separate. `state` keeps the
/// truth available for colouring and filtering without pretending the assignment is measured.
/// The template ignores JSON keys it does not know, so this is safe to emit unconditionally.
///
/// 46Ar NEEDS A BEAM-ENERGY PATCH THAT THE TEMPLATE DOES NOT HAVE, and it is not cosmetic.
/// The template evaluates kine2b(s.ebeam, ...) with ONE constant beam energy for every event, and
/// its E_beam slider runs 80-260 MeV. 46Ar enters at 598 MeV and loses 0.957 MeV/cm, i.e. 95.7 MeV
/// across the drift, so a constant value would smear every kinematic locus into a wide band and
/// no slider position would ever be right. explorer/patch_ar46.py adds the per-event
/// E_beam(z) = ebeam0 - dEdz * vertexz and widens the slider. Generate, then patch -- the launcher
/// does both and refuses to hand over an unpatched page.
#include "../ar46_masses.h"

#include <fstream>
#include <sstream>

/// Stream one cache as the template's per-set JSON. Column names and the float reads are the
/// contract, not a local choice: kin_Ar46.C fills a TNtuple, which is float-only, and reading a
/// float branch as a double yields garbage rather than an error.
static Long64_t dump_pk_ar46(TTree *t, FILE *o)
{
   float ke = 0, th = 0, c2 = 0, vz = 0, st = -1;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("chi2ndf", &c2);
   const bool hasVz = t->GetBranch("vertexz") != nullptr;
   if (hasVz) t->SetBranchAddress("vertexz", &vz);
   const bool hasSt = t->GetBranch("state") != nullptr;
   if (hasSt) t->SetBranchAddress("state", &st);
   const Long64_t N = t->GetEntries();
   fprintf(o, "{\"ke\":[");
   for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); fprintf(o, "%s%.3f", i ? "," : "", ke); }
   fprintf(o, "],\"th\":[");
   for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); fprintf(o, "%s%.3f", i ? "," : "", th); }
   fprintf(o, "],\"c2\":[");
   for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); fprintf(o, "%s%.2f", i ? "," : "", c2); }
   fprintf(o, "],\"vz\":[");
   for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); fprintf(o, "%s%.1f", i ? "," : "", hasVz ? vz : 0.f); }
   if (hasSt) {
      fprintf(o, "],\"state\":[");
      for (Long64_t i = 0; i < N; ++i) { t->GetEntry(i); fprintf(o, "%s%.3f", i ? "," : "", st); }
   }
   fprintf(o, "]}");
   return N;
}

void make_explorer_Ar46(
   // name=cache pairs, in the order the buttons should appear. The names are what the page shows.
   TString setsCSV = "2.85T AT-TPC=../plots/kin_b285_attpc.root,3.9T AT-TPC=../plots/kin_b39_attpc.root,"
                     "2.85T 2mm=../plots/kin_b285_2mm.root,3.9T 2mm=../plots/kin_b39_2mm.root",
   TString outHtml = "", TString tag = "46Ar(3He,d)47K",
   // Masses in amu, DERIVED FROM THE MeV VALUES IN ex_core_3Hed.h (divided by u = 931.49401) so
   // the page's kine2b and the macros' inversion cannot disagree. These are NUCLEAR masses.
   // MASSES ARE DERIVED, NOT TYPED. Passing <= 0 (the default) computes them from ar46_masses.h.
   // Retyping them by hand is precisely what broke this page once: a beam mass out by 2.66 MeV
   // and a residual out by 1.60 gave Q = +3.521 instead of +7.777, i.e. every excitation energy
   // on the page shifted by -4.26 MeV, while every macro stayed correct.
   double ebeam0 = -1, double mBeamAmu = -1, double mTargAmu = -1, double mEjectAmu = -1,
   double mResidAmu = -1, int beamA = 46,
   // The three levels of the proposal, drawn as kinematic loci.
   TString refExCSV = "0:g.s. 1/2+ 1s1/2,0.36:3/2+ 0d3/2,2.02:7/2- 0f7/2",
   TString tplPath = "")
{
   // AN EMPTY setsCSV IS THE DEFAULT, NOT "no sets". A caller passing "" to mean "use the
   // default" would otherwise get a page with zero configurations that still writes, still opens
   // and still reports success -- which is exactly what the launcher did the first time it ran.
   if (setsCSV.IsNull())
      setsCSV = "2.85T AT-TPC=../plots/kin_b285_attpc.root,3.9T AT-TPC=../plots/kin_b39_attpc.root,"
                "2.85T 2mm=../plots/kin_b285_2mm.root,3.9T 2mm=../plots/kin_b39_2mm.root";
   if (ebeam0    <= 0) ebeam0    = Ar46::kTb0;
   if (mBeamAmu  <= 0) mBeamAmu  = Ar46::kMb / Ar46::kU;
   if (mTargAmu  <= 0) mTargAmu  = Ar46::kMt / Ar46::kU;
   if (mEjectAmu <= 0) mEjectAmu = Ar46::kMe / Ar46::kU;
   if (mResidAmu <= 0) mResidAmu = Ar46::kMR / Ar46::kU;
   // Assert the Q value the PAGE will compute from the amu it is handed, not the one this file
   // knows. A rounding of the amu that shifts Q is exactly the failure being guarded against.
   const double qPage = (mBeamAmu + mTargAmu - mEjectAmu - mResidAmu) * Ar46::kU;
   printf("  Q(g.s.) as the page will compute it: %+.4f MeV   (expected %+.4f)\n", qPage, Ar46::kQgs);
   if (std::fabs(qPage - Ar46::kQgs) > 0.01) {
      std::cerr << "MASS MISMATCH: the page would shift every E_x by " << (qPage - Ar46::kQgs)
                << " MeV -- refusing to write it\n";
      return;
   }
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (tplPath.IsNull())
      tplPath = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954_Be12/UKF/pp/"
                "explorer_template.html";
   if (outHtml.IsNull()) outHtml = TString(gSystem->Getenv("HOME")) + "/ar46_3Hed_explorer.html";

   std::ifstream in(tplPath.Data());
   if (!in) { std::cerr << "cannot read template " << tplPath << "\n"; return; }
   std::stringstream ss; ss << in.rdbuf();
   std::string tpl = ss.str();

   // ---- reference loci ----------------------------------------------------
   TString refJson = "[";
   std::unique_ptr<TObjArray> toks(refExCSV.Tokenize(","));
   for (int i = 0; i < toks->GetEntries(); ++i) {
      TString tk = ((TObjString *)toks->At(i))->GetString();
      TString exs = tk, lab = "";
      if (tk.Contains(":")) { exs = tk(0, tk.First(':')); lab = tk(tk.First(':') + 1, tk.Length()); }
      double exv = exs.Atof();
      TString label = TString::Format("E_x = %g", exv);
      if (lab.Length()) label += " (" + lab + ")";
      refJson += TString::Format("%s{\"ex\":%g,\"label\":\"%s\"}", i ? "," : "", exv, label.Data());
   }
   refJson += "]";

   // ---- open every cache BEFORE writing anything --------------------------
   // A half-written page whose later sets are missing looks like a working page with fewer
   // configurations, so every input is opened and checked first and the run aborts on any failure.
   std::unique_ptr<TObjArray> sets(setsCSV.Tokenize(","));
   std::vector<TString> names;
   std::vector<TFile *> files;
   std::vector<TTree *> trees;
   for (int i = 0; i < sets->GetEntries(); ++i) {
      TString one = ((TObjString *)sets->At(i))->GetString();
      if (!one.Contains("=")) { std::cerr << "set '" << one << "' is not name=cache\n"; return; }
      TString nm = one(0, one.First('=')), fn = one(one.First('=') + 1, one.Length());
      if (!fn.BeginsWith("/")) fn = here + "/" + fn;
      TFile *f = TFile::Open(fn);
      TTree *t = (f && !f->IsZombie()) ? (TTree *)f->Get("pk") : nullptr;
      if (!t) { std::cerr << "cannot read pk from " << fn << " -- run kin_Ar46.C for it first\n"; return; }
      names.push_back(nm); files.push_back(f); trees.push_back(t);
   }

   if (names.empty()) { std::cerr << "no sets to write -- refusing to produce an empty page\n"; return; }

   TString cacheNames;
   for (size_t i = 0; i < names.size(); ++i) cacheNames += (i ? ", " : "") + names[i];

   TString cfg = TString::Format(
      "{\"tag\":\"%s\",\"title\":\"%s kinematics explorer\",\"eyebrow\":\"46Ar(3He,d) proposal . AT-TPC sim\","
      "\"cache\":\"%s\",\"pngName\":\"explorer_ar46\",\"ebeam0\":%.4f,\"beamA\":%d,"
      "\"mBeamAmu\":%.6f,\"mTargAmu\":%.6f,\"mEjectAmu\":%.6f,\"mResidAmu\":%.6f,"
      "\"dEdz\":%.6f,\"driftLength\":%.1f,"
      "\"refEx\":%s}",
      tag.Data(), tag.Data(), cacheNames.Data(), ebeam0, beamA, mBeamAmu, mTargAmu, mEjectAmu, mResidAmu,
      Ar46::kdEdz / 10.0, 1000.0,   // MeV per MILLIMETRE, and the drift in mm: the cache and the
                                    // template are both mm, so dEdz must be too.
      refJson.Data());

   const std::string cfgKey = "/*__CFG__*/ null", datKey = "/*__DATA__*/ null";
   size_t pc = tpl.find(cfgKey);
   if (pc == std::string::npos) { std::cerr << "template has no CFG placeholder\n"; return; }
   tpl.replace(pc, cfgKey.size(), cfg.Data());
   size_t pd = tpl.find(datKey);
   if (pd == std::string::npos) { std::cerr << "template has no DATA placeholder\n"; return; }

   FILE *o = fopen(outHtml.Data(), "w");
   if (!o) { std::cerr << "cannot write " << outHtml << "\n"; return; }
   fwrite(tpl.data(), 1, pd, o);
   fprintf(o, "{");
   std::vector<Long64_t> n(names.size());
   for (size_t i = 0; i < names.size(); ++i) {
      fprintf(o, "%s\"%s\":", i ? "," : "", names[i].Data());
      n[i] = dump_pk_ar46(trees[i], o);
   }
   fprintf(o, "}");
   fwrite(tpl.data() + pd + datKey.size(), 1, tpl.size() - pd - datKey.size(), o);
   fclose(o);
   for (auto *f : files) f->Close();

   FileStat_t st;
   gSystem->GetPathInfo(outHtml, st);
   printf("\nwrote %s  (%.1f MB)\n", outHtml.Data(), st.fSize / 1048576.);
   for (size_t i = 0; i < names.size(); ++i) printf("  %-16s : %lld tracks\n", names[i].Data(), n[i]);
   printf("\nNOT YET USABLE -- run explorer/patch_ar46.py on it for the E_beam(z) term and the\n"
          "598 MeV slider range, or every locus is drawn at one beam energy the data never has.\n");
}
