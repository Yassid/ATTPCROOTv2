/// @file make_explorer_html.C
/// @brief Bake kinematics caches (proton_kin*.root from pp/ex_Be12.C) into a single
/// self-contained HTML file: the same explorer as pp/explore_Be12.C, but running in a
/// browser instead of X11 -- useful when WSLg/X forwarding is not cooperating.
///
/// The page recomputes Ex and theta_cm from (KE, theta_lab) in JavaScript with the very
/// same two-body expressions, so the beam energy, the cuts and every binning stay live.
/// Pass BOTH a UKF and a GENFIT cache to get the in-page fitter switch (and the
/// "overlay the other fitter" comparison); one cache alone still works.
///
///   root -b -q 'pp/make_explorer_html.C()'      // the (p,p') 300 torr UKF cache
///   ./open_explorer.sh pp|pd|pt                 // what you actually want: it wires the
///                                               // caches, masses and level scheme per channel
///
/// PORTED 2026-09-02 from C15d/dd + a2091 (repo HEAD "C15d viewer: run range, IC multiplicity,
/// and a per-run Ex map"): this replaces the 2026-07-27 fork, which had no third set, no run /
/// IC / npulse columns and no in-page KE-correction panel.
///
/// Then open it: pp/open_explorer.sh, or the Windows browser on the written file.

#include <string>

/// stream one cache's (ke, theta, vertex z, chi2/ndf) columns as a compact JSON object
static Long64_t dump_pk(TTree *t, FILE *o)
{
   float ke = 0, th = 0, c2 = 0, vz = 0;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("chi2ndf", &c2);
   // vertexz is absent from very old caches, and was hardcoded to 0 in ex_C15.C before
   // 2026-07-26 -- the page falls back to hiding the z panel when the column is flat.
   const bool hasVz = t->GetBranch("vertexz") != nullptr;
   if (hasVz)
      t->SetBranchAddress("vertexz", &vz);
   // ic / npulse / run are optional and the page hides their controls when the column is absent,
   // rather than showing one that selects nothing. EVERY column here is read as a float on
   // purpose: ex_Be12.C fills a TNtuple, which is float-only -- the C15d generator reads doubles
   // because ITS cache is a TTree of doubles, and reading one as the other yields garbage.
   float ic = -1, npulse = 0, runNo = 0;
   const bool hasIc = t->GetBranch("ic") != nullptr;
   const bool hasNp = t->GetBranch("npulse") != nullptr;
   const bool hasRun = t->GetBranch("run") != nullptr;
   if (hasIc)
      t->SetBranchAddress("ic", &ic);
   if (hasNp)
      t->SetBranchAddress("npulse", &npulse);
   if (hasRun)
      t->SetBranchAddress("run", &runNo);
   const Long64_t N = t->GetEntries();
   fprintf(o, "{\"ke\":[");
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      fprintf(o, "%s%.3f", i ? "," : "", ke);
   }
   fprintf(o, "],\"th\":[");
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      fprintf(o, "%s%.3f", i ? "," : "", th);
   }
   fprintf(o, "],\"c2\":[");
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      fprintf(o, "%s%.2f", i ? "," : "", c2);
   }
   fprintf(o, "],\"vz\":[");
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      fprintf(o, "%s%.1f", i ? "," : "", hasVz ? vz : 0.f);
   }
   if (hasIc) {
      fprintf(o, "],\"ic\":[");
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         fprintf(o, "%s%.0f", i ? "," : "", ic);
      }
   }
   if (hasNp) {
      fprintf(o, "],\"np\":[");
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         fprintf(o, "%s%.0f", i ? "," : "", npulse);
      }
   }
   if (hasRun) {
      fprintf(o, "],\"run\":[");
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         fprintf(o, "%s%.0f", i ? "," : "", runNo);
      }
   }
   fprintf(o, "]}");
   return N;
}

void make_explorer_html(TString cache = "", TString outHtml = "", TString tag = "12Be(p,p')",
                        // 155 MeV is a PLACEHOLDER, not a calibration: it came from scanning
                        // E_beam until the elastic peak sat at zero, in-sample and against a gas
                        // density that was wrong by 2x. The page recomputes Ex live, so move it.
                        double ebeam0 = 155.0, double mBeamAmu = 12.026921, double mTargAmu = 1.007825,
                        double mEjectAmu = 1.007825, double mResidAmu = 12.026921, int beamA = 12,
                        TString refExCSV = "", TString cacheGenfit = "", TString cacheThird = "",
                        TString setNamesCSV = "")
{
   // cacheThird / setNamesCSV are additive and default to empty, so every existing caller gets
   // exactly the old two-set page. setNamesCSV renames the JSON keys, which are what the page
   // shows on its selector buttons -- "ukf"/"genfit" is meaningless once the three sets are
   // genfit-without-material-effects, genfit-with, and the UKF.
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   // default to the GATED caches (IC 15C beam + proton PID gate); fall back to whatever
   // ungated cache exists if the gated pass has not been run
   if (cache.IsNull()) {
      // the 2026-08-25 rebuild (300 torr + CATIMA, IC 500-800) wrote the *800* tags
      cache = here + "/plots/proton_kin_pp800_ukf.root";
      if (gSystem->AccessPathName(cache)) cache = here + "/plots/proton_kin_300_ukf.root";
   }
   // Derive the GENFIT sibling FROM THE UKF CACHE NAME -- never from a fixed path. Hardcoding
   // proton_kin_g_genfit_nomat.root here silently paired the (p,p') GENFIT protons with the
   // (p,d) UKF deuterons: the page offered a "fitter switch" between two different reactions,
   // both interpreted with the (p,d) masses. Deriving the name keeps the pair in one channel.
   // NB: an empty TString is IsNull(), so passing "" explicitly also lands here -- hence the
   // sibling must be channel-correct rather than merely present on disk.
   if (cacheGenfit.IsNull() && cache.Contains("_ukf")) {
      for (const char *suffix : {"_genfit_nomat", "_genfit"}) {
         TString g = cache;
         g.ReplaceAll("_ukf", suffix);
         if (g != cache && !gSystem->AccessPathName(g)) { cacheGenfit = g; break; }
      }
   }
   if (outHtml.IsNull())
      outHtml = TString(gSystem->Getenv("HOME")) + "/a1954_Be12_explorer.html";
   TString tplPath = here + "/explorer_template.html";

   // ---- template ----------------------------------------------------------
   std::ifstream tf(tplPath.Data());
   if (!tf) {
      std::cerr << "cannot read template " << tplPath << "\n";
      return;
   }
   std::string tpl((std::istreambuf_iterator<char>(tf)), std::istreambuf_iterator<char>());
   tf.close();

   // ---- caches ------------------------------------------------------------
   TFile *fU = TFile::Open(cache);
   TTree *tU = (fU && !fU->IsZombie()) ? (TTree *)fU->Get("pk") : nullptr;
   if (!tU) {
      std::cerr << "cannot read ntuple 'pk' from " << cache << "\n";
      return;
   }
   TFile *fG = cacheGenfit.Length() ? TFile::Open(cacheGenfit) : nullptr;
   TTree *tG = (fG && !fG->IsZombie()) ? (TTree *)fG->Get("pk") : nullptr;
   if (cacheGenfit.Length() && !tG)
      std::cerr << "WARNING: genfit cache unusable (" << cacheGenfit << ") -- writing UKF only\n";

   TFile *fT = cacheThird.Length() ? TFile::Open(cacheThird) : nullptr;
   TTree *tT = (fT && !fT->IsZombie()) ? (TTree *)fT->Get("pk") : nullptr;
   if (cacheThird.Length() && !tT)
      std::cerr << "WARNING: third cache unusable (" << cacheThird << ") -- skipping it\n";

   // JSON keys for the three slots. A missing/short CSV falls back to the historical names.
   TString sn[3] = {"ukf", "genfit", "third"};
   if (setNamesCSV.Length()) {
      std::unique_ptr<TObjArray> parts(setNamesCSV.Tokenize(","));
      for (int i = 0; i < 3 && i < parts->GetEntries(); ++i) {
         TString v = ((TObjString *)parts->At(i))->GetString();
         v = v.Strip(TString::kBoth);
         if (v.Length())
            sn[i] = v;
      }
   }

   // ---- reference levels drawn as kinematic loci (per channel) ------------
   //  "Ex:label" pairs; default set chosen from the residual nucleus
   // 15C level scheme, NOT 14C's. The ported default listed 6.09 1-, 6.59 0+_2, 7.34 2+,
   // which are 14C levels -- wrong nucleus, and they would have drawn loci at energies where
   // 15C has no bound structure at all. 15C has exactly ONE bound excited state, the 5/2+ at
   // 0.740 MeV, and is neutron-unbound above Sn = 1.218 MeV; the Sn marker is worth showing
   // because everything above it is continuum rather than a level.
   // 14C set includes the 6.73 3- so the whole 6-7 MeV cluster seen in 15C(p,d)14C is marked,
   // and Sn = 8.18 because the structure above it (the ~9.3 MeV bump) is continuum, not a level.
   if (refExCSV.IsNull())
      refExCSV = (mResidAmu > 11.5) ? "0:g.s.,2.10:2+_1,4.56:,5.7:"                        // 12Be
               : (mResidAmu > 10.5) ? "0:g.s. 1/2+,0.320:1/2-,1.778:5/2+,2.654:3/2-"       // 11Be
                                    : "0:g.s. 0+,3.368:2+_1,5.958:2+_2,6.179:0+_2,"
                                      "6.263:2-,6.812:Sn,7.371:3-,7.542:2+_3,9.27:4-";     // 10Be
   TString refJson = "[";
   TObjArray *toks = refExCSV.Tokenize(",");
   for (int i = 0; i < toks->GetEntries(); ++i) {
      TString tk = ((TObjString *)toks->At(i))->GetString();
      TString exs = tk, lab = "";
      if (tk.Contains(":")) {
         exs = tk(0, tk.First(':'));
         lab = tk(tk.First(':') + 1, tk.Length());
      }
      double exv = exs.Atof();
      TString label = TString::Format("E_x = %g", exv);
      if (lab.Length())
         label += " (" + lab + ")";
      refJson += TString::Format("%s{\"ex\":%g,\"label\":\"%s\"}", i ? "," : "", exv, label.Data());
   }
   refJson += "]";

   // ---- config block ------------------------------------------------------
   TString cacheNames = gSystem->BaseName(cache);
   if (tG)
      cacheNames += TString(" + ") + gSystem->BaseName(cacheGenfit);
   TString cfg = TString::Format(
      "{\"tag\":\"%s\",\"title\":\"%s excitation explorer\",\"eyebrow\":\"a1954 . AT-TPC\","
      "\"cache\":\"%s\",\"pngName\":\"explorer_%s\",\"ebeam0\":%.4f,\"beamA\":%d,"
      "\"mBeamAmu\":%.6f,\"mTargAmu\":%.6f,\"mEjectAmu\":%.6f,\"mResidAmu\":%.6f,"
      "\"refEx\":%s}",
      tag.Data(), tag.Data(), cacheNames.Data(), gSystem->BaseName(outHtml), ebeam0, beamA, mBeamAmu, mTargAmu,
      mEjectAmu, mResidAmu, refJson.Data());

   const std::string cfgKey = "/*__CFG__*/ null", datKey = "/*__DATA__*/ null";
   size_t pc = tpl.find(cfgKey);
   if (pc == std::string::npos) {
      std::cerr << "template has no CFG placeholder\n";
      return;
   }
   tpl.replace(pc, cfgKey.size(), cfg.Data());
   size_t pd = tpl.find(datKey);
   if (pd == std::string::npos) {
      std::cerr << "template has no DATA placeholder\n";
      return;
   }

   // ---- write: head + streamed data + tail --------------------------------
   FILE *o = fopen(outHtml.Data(), "w");
   if (!o) {
      std::cerr << "cannot write " << outHtml << "\n";
      return;
   }
   fwrite(tpl.data(), 1, pd, o);
   fprintf(o, "{\"%s\":", sn[0].Data());
   Long64_t nU = dump_pk(tU, o), nG = 0, nT = 0;
   if (tG) {
      fprintf(o, ",\"%s\":", sn[1].Data());
      nG = dump_pk(tG, o);
   }
   if (tT) {
      fprintf(o, ",\"%s\":", sn[2].Data());
      nT = dump_pk(tT, o);
   }
   fprintf(o, "}");
   fwrite(tpl.data() + pd + datKey.size(), 1, tpl.size() - pd - datKey.size(), o);
   fclose(o);
   fU->Close();
   if (fG)
      fG->Close();
   if (fT)
      fT->Close();

   FileStat_t st;
   gSystem->GetPathInfo(outHtml, st);
   printf("\nwrote %s  (%.1f MB)\n", outHtml.Data(), st.fSize / 1048576.);
   printf("  %-10s : %lld tracks\n", sn[0].Data(), nU);
   if (tG)
      printf("  %-10s : %lld tracks   (fitter switch enabled in the page)\n", sn[1].Data(), nG);
   else
      printf("  %-10s : not supplied -- that button stays disabled\n", sn[1].Data());
   if (tT)
      printf("  %-10s : %lld tracks\n", sn[2].Data(), nT);
   printf("open it with:  pp/open_explorer.sh   (or the Windows browser on %s)\n\n", outHtml.Data());
}
