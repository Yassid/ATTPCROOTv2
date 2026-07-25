/// @file make_explorer_html.C
/// @brief Bake a kinematics cache (proton_kin*.root from pp/ex_Be12.C) into a single
/// self-contained HTML file: the same explorer as pp/explore_Be12.C, but running in a
/// browser instead of X11 -- useful when WSLg/X forwarding is not cooperating.
///
/// The page recomputes Ex and theta_cm from (KE, theta_lab) in JavaScript with the very
/// same two-body expressions, so the beam energy, the cuts and every binning stay live.
///
///   root -b -q 'pp/make_explorer_html.C()'                       // clean155 cache, 155 MeV
///   root -b -q 'pp/make_explorer_html.C("plots/proton_kin_pd.root","/home/yassid/pd.html",\
///                "12Be(p,d)11Be",155,12.026921,1.007825,2.014102,11.021658,12)'
///
/// Then just open the file (WSL -> Windows browser):   explorer.exe <file>.html

#include <string>

void make_explorer_html(TString cache = "", TString outHtml = "", TString tag = "12Be(p,p')",
                        double ebeam0 = 155.0, double mBeamAmu = 12.026921, double mTargAmu = 1.007825,
                        double mEjectAmu = 1.007825, double mResidAmu = 12.026921, int beamA = 12)
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (cache.IsNull())
      cache = here + "/plots/proton_kin_clean155.root";
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

   // ---- cache -------------------------------------------------------------
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      std::cerr << "cannot open cache " << cache << "\n";
      return;
   }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) {
      std::cerr << "no ntuple 'pk' in " << cache << "\n";
      return;
   }
   float ke = 0, th = 0, c2 = 0;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("chi2ndf", &c2);
   const Long64_t N = t->GetEntries();

   // ---- config block ------------------------------------------------------
   TString cfg = TString::Format(
      "{\"tag\":\"%s\",\"title\":\"%s excitation explorer\",\"eyebrow\":\"a1954 . AT-TPC . UKF\","
      "\"cache\":\"%s\",\"pngName\":\"explorer_%s\",\"ebeam0\":%.4f,\"beamA\":%d,"
      "\"mBeamAmu\":%.6f,\"mTargAmu\":%.6f,\"mEjectAmu\":%.6f,\"mResidAmu\":%.6f}",
      tag.Data(), tag.Data(), gSystem->BaseName(cache), gSystem->BaseName(outHtml), ebeam0, beamA, mBeamAmu, mTargAmu,
      mEjectAmu, mResidAmu);

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
   fprintf(o, "]}");
   fwrite(tpl.data() + pd + datKey.size(), 1, tpl.size() - pd - datKey.size(), o);
   fclose(o);
   f->Close();

   FileStat_t st;
   gSystem->GetPathInfo(outHtml, st);
   printf("\nwrote %s  (%lld tracks, %.1f MB)\n", outHtml.Data(), N, st.fSize / 1048576.);
   printf("open it with:  pp/open_explorer.sh   (or the Windows browser on %s)\n\n", outHtml.Data());
}
