/// @file make_explorer_html.C
/// @brief Bake kinematics caches (proton_kin*.root from pp/ex_C15.C) into a single
/// self-contained HTML file: the same explorer as pp/explore_C15.C, but running in a
/// browser instead of X11 -- useful when WSLg/X forwarding is not cooperating.
///
/// The page recomputes Ex and theta_cm from (KE, theta_lab) in JavaScript with the very
/// same two-body expressions, so the beam energy, the cuts and every binning stay live.
/// Pass BOTH a UKF and a GENFIT cache to get the in-page fitter switch (and the
/// "overlay the other fitter" comparison); one cache alone still works.
///
///   root -b -q 'pp/make_explorer_html.C()'      // clean195 UKF cache, 161 MeV
///   root -b -q 'pp/make_explorer_html.C("plots/proton_kin_clean195.root","/home/yassid/pp.html",\
///                "15C(p,p'\'')",161,15.0105993,1.007825,1.007825,15.0105993,12,"",\
///                "plots/proton_kin_clean195_genfit.root")'
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
   fprintf(o, "]}");
   return N;
}

void make_explorer_html(TString cache = "", TString outHtml = "", TString tag = "15C(p,p')",
                        double ebeam0 = 195.0, double mBeamAmu = 15.0105993, double mTargAmu = 1.007825,
                        double mEjectAmu = 1.007825, double mResidAmu = 15.0105993, int beamA = 14,
                        TString refExCSV = "", TString cacheGenfit = "")
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (cache.IsNull())
      cache = here + "/plots/proton_kin_clean195.root";
   if (outHtml.IsNull())
      outHtml = TString(gSystem->Getenv("HOME")) + "/a2091_C15_explorer.html";
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

   // ---- reference levels drawn as kinematic loci (per channel) ------------
   //  "Ex:label" pairs; default set chosen from the residual nucleus
   if (refExCSV.IsNull())
      refExCSV = (mResidAmu > 13.5) ? "0:g.s.,6.09:1-,6.59:0+_2,7.34:2+"   // 15C
                                    : "0:g.s.,3.09:1/2+,3.68:3/2-,3.85:5/2+"; // 13C
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
      "{\"tag\":\"%s\",\"title\":\"%s excitation explorer\",\"eyebrow\":\"a2091 . AT-TPC\","
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
   fprintf(o, "{\"ukf\":");
   Long64_t nU = dump_pk(tU, o), nG = 0;
   if (tG) {
      fprintf(o, ",\"genfit\":");
      nG = dump_pk(tG, o);
   }
   fprintf(o, "}");
   fwrite(tpl.data() + pd + datKey.size(), 1, tpl.size() - pd - datKey.size(), o);
   fclose(o);
   fU->Close();
   if (fG)
      fG->Close();

   FileStat_t st;
   gSystem->GetPathInfo(outHtml, st);
   printf("\nwrote %s  (%.1f MB)\n", outHtml.Data(), st.fSize / 1048576.);
   printf("  UKF    : %lld tracks\n", nU);
   if (tG)
      printf("  GENFIT : %lld tracks   (fitter switch enabled in the page)\n", nG);
   else
      printf("  GENFIT : not supplied -- the page's GENFIT button stays disabled\n");
   printf("open it with:  pp/open_explorer.sh   (or the Windows browser on %s)\n\n", outHtml.Data());
}
