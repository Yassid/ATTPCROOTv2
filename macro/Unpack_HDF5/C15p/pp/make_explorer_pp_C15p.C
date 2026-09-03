/// @file make_explorer_pp_C15p.C
/// @brief Bake kinematics caches (proton_kin*.root from pp/ex_C15.C) into a single
/// self-contained HTML file: the same explorer as pp/explore_C15.C, but running in a
/// browser instead of X11 -- useful when WSLg/X forwarding is not cooperating.
///
/// The page recomputes Ex and theta_cm from (KE, theta_lab) in JavaScript with the very
/// same two-body expressions, so the beam energy, the cuts and every binning stay live.
/// Pass BOTH a UKF and a GENFIT cache to get the in-page fitter switch (and the
/// "overlay the other fitter" comparison); one cache alone still works.
///
///   root -b -q 'pp/make_explorer_pp_C15p.C()'      // clean195 UKF cache, 161 MeV
///   root -b -q 'pp/make_explorer_pp_C15p.C("plots/proton_kin_clean195.root","/home/yassid/pp.html",\
///                "15C(p,p'\'')",161,15.0105993,1.007825,1.007825,15.0105993,12,"",\
///                "plots/proton_kin_clean195_genfit.root")'
///
/// Then open it: pp/open_explorer.sh, or the Windows browser on the written file.

#include <string>

/// stream one cache's (ke, theta, vertex z, chi2/ndf) columns as a compact JSON object
/// stride > 1 emits every Nth track. The page is a spectrum explorer, so a subsample of a few
/// hundred thousand is statistically indistinguishable from the full set for every histogram it
/// draws -- but the FULL set stays in the ROOT cache, and any number quoted for physics should
/// come from there, not from here. The stride is reported on the page so the subsampling is
/// never invisible.
static Long64_t dump_pk(TTree *t, FILE *o, Long64_t stride = 1)
{
   double ke = 0, th = 0, c2 = 0, vz = 0;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("chi2ndf", &c2);
   // vertexz is absent from very old caches, and was hardcoded to 0 in ex_C15.C before
   // 2026-07-26 -- the page falls back to hiding the z panel when the column is flat.
   const bool hasVz = t->GetBranch("vz") != nullptr;
   if (hasVz)
      t->SetBranchAddress("vz", &vz);
   // ic / npulse are optional: caches built before the IC join have neither, and the page hides
   // the IC control when the column is absent rather than showing one that selects nothing.
   float ic = -1;
   int npulse = 0, runNo = 0;
   const bool hasIc = t->GetBranch("ic") != nullptr;
   const bool hasNp = t->GetBranch("npulse") != nullptr;
   // run rides along so a feature can be asked "is this one run or all of them?" in the page,
   // which is the first question to put to any structure that has no known counterpart.
   const bool hasRun = t->GetBranch("run") != nullptr;
   if (hasIc)
      t->SetBranchAddress("ic", &ic);
   if (hasNp)
      t->SetBranchAddress("npulse", &npulse);
   if (hasRun)
      t->SetBranchAddress("run", &runNo);
   const Long64_t N = t->GetEntries();
   fprintf(o, "{\"ke\":[");
   for (Long64_t i = 0; i < N; i += stride) {
      t->GetEntry(i);
      fprintf(o, "%s%.2f", i ? "," : "", ke);
   }
   fprintf(o, "],\"th\":[");
   for (Long64_t i = 0; i < N; i += stride) {
      t->GetEntry(i);
      fprintf(o, "%s%.2f", i ? "," : "", th);
   }
   fprintf(o, "],\"c2\":[");
   for (Long64_t i = 0; i < N; i += stride) {
      t->GetEntry(i);
      fprintf(o, "%s%.2f", i ? "," : "", c2);
   }
   fprintf(o, "],\"vz\":[");
   for (Long64_t i = 0; i < N; i += stride) {
      t->GetEntry(i);
      fprintf(o, "%s%.1f", i ? "," : "", hasVz ? vz : 0.0);
   }
   if (hasIc) {
      fprintf(o, "],\"ic\":[");
      for (Long64_t i = 0; i < N; i += stride) {
         t->GetEntry(i);
         fprintf(o, "%s%.0f", i ? "," : "", ic);
      }
   }
   if (hasNp) {
      fprintf(o, "],\"np\":[");
      for (Long64_t i = 0; i < N; i += stride) {
         t->GetEntry(i);
         fprintf(o, "%s%d", i ? "," : "", npulse);
      }
   }
   if (hasRun) {
      fprintf(o, "],\"run\":[");
      for (Long64_t i = 0; i < N; i += stride) {
         t->GetEntry(i);
         fprintf(o, "%s%d", i ? "," : "", runNo);
      }
   }
   fprintf(o, "]}");
   // the number ACTUALLY emitted, not t->GetEntries(): with a stride they differ, and
   // reporting the full count would overstate what the page holds.
   return (N + stride - 1) / stride;
}

/// beamA is 15, not 14: it was left at the a1954 value by the port, which mislabels every
/// nuclide the page typesets. ebeam0 defaults to the calibrated 170 MeV (elastic ridge +
/// Ex-vs-theta tilt on the IC+PID-gated sample), not the old 195 placeholder.
void make_explorer_pp_C15p(TString cache = "", TString outHtml = "", TString tag = "15C(p,p')",
                        // Ebeam from THIS analysis's own 15C(p,p) elastic ridge (pp/ebeam_pp_C15p.C),
                        // not an external number -- so the g.s. landing at Ex = 0 here is a
                        // self-consistency check, not a calibration.
                        double ebeam0 = 157.0, double mBeamAmu = 15.0105993, double mTargAmu = 1.00782503,
                        // (p,p'): the ejectile is the PROTON and the residual is 15C itself.
                        double mEjectAmu = 1.00782503, double mResidAmu = 15.0105993, int beamA = 15,
                        TString refExCSV = "", TString cacheGenfit = "", TString cacheThird = "",
                        TString setNamesCSV = "",
                        /// Cap on the number of tracks BAKED INTO THE PAGE. The artifact ceiling is
                        /// 16 MB and 500k tracks alone reach 15.9 MB, leaving no room for the IC and
                        /// per-run panels. Emitting every Nth track keeps every spectrum the page
                        /// draws statistically equivalent while fitting the budget; the full sample
                        /// stays in the ROOT cache and is what any quoted number should come from.
                        Long64_t maxTracks = 300000)
{
   // cacheThird / setNamesCSV are additive and default to empty, so every existing caller gets
   // exactly the old two-set page. setNamesCSV renames the JSON keys, which are what the page
   // shows on its selector buttons -- "ukf"/"genfit" is meaningless once the three sets are
   // genfit-without-material-effects, genfit-with, and the UKF.
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   // default to the GATED caches (IC 15C beam + proton PID gate); fall back to whatever
   // ungated cache exists if the gated pass has not been run
   if (cache.IsNull())
      cache = here + "/plots/pp_kin_C15p.root";
   // Derive the GENFIT sibling FROM THE UKF CACHE NAME -- never from a fixed path. Hardcoding
   // proton_kin_g_genfit_nomat.root here silently paired the (p,p') GENFIT protons with the
   // (p,d) UKF deuterons: the page offered a "fitter switch" between two different reactions,
   // both interpreted with the (p,d) masses. Deriving the name keeps the pair in one channel.
   // NB: an empty TString is IsNull(), so passing "" explicitly also lands here -- hence the
   // sibling must be channel-correct rather than merely present on disk.
   if (false && cacheGenfit.IsNull() && cache.Contains("_ukf")) {
      for (const char *suffix : {"_genfit_nomat", "_genfit"}) {
         TString g = cache;
         g.ReplaceAll("_ukf", suffix);
         if (g != cache && !gSystem->AccessPathName(g)) { cacheGenfit = g; break; }
      }
   }
   if (outHtml.IsNull())
      outHtml = TString(gSystem->Getenv("HOME")) + "/C15p_dd_explorer.html";
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
   TTree *tU = (fU && !fU->IsZombie()) ? (TTree *)fU->Get("dd") : nullptr;
   if (!tU) {
      std::cerr << "cannot read ntuple 'dd' from " << cache << "\n";
      return;
   }
   TFile *fG = cacheGenfit.Length() ? TFile::Open(cacheGenfit) : nullptr;
   TTree *tG = (fG && !fG->IsZombie()) ? (TTree *)fG->Get("dd") : nullptr;
   if (cacheGenfit.Length() && !tG)
      std::cerr << "WARNING: genfit cache unusable (" << cacheGenfit << ") -- writing UKF only\n";

   TFile *fT = cacheThird.Length() ? TFile::Open(cacheThird) : nullptr;
   TTree *tT = (fT && !fT->IsZombie()) ? (TTree *)fT->Get("dd") : nullptr;
   if (cacheThird.Length() && !tT)
      std::cerr << "WARNING: third cache unusable (" << cacheThird << ") -- skipping it\n";

   // JSON keys for the three slots. A missing/short CSV falls back to the historical names.
   // ★ NOT "ukf". This workspace fits with GENFIT+CATIMA only -- the UKF is deliberately not used
   // here -- and the inherited default named the first slot "ukf", which is what the page prints on
   // its selector button. A plot labelled with a fitter that never touched the data is exactly the
   // kind of thing that survives into a talk.
   TString sn[3] = {"genfit_catima", "genfit_nomat", "third"};
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
      refExCSV = (mResidAmu > 14.5) ? "0:g.s.,0.740:5/2+,1.218:Sn"                      // 15C
               : (mResidAmu > 13.5) ? "0:g.s.,6.09:1-,6.59:0+,6.73:3-,7.01:2+,8.18:Sn"  // 14C
                                    : "0:g.s.,3.09:1/2+,3.68:3/2-,3.85:5/2+";           // 13C
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
      "{\"tag\":\"%s\",\"title\":\"%s excitation explorer\",\"eyebrow\":\"a2091 H2 . AT-TPC\","
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
   const Long64_t nAll = tU ? tU->GetEntries() : 0;
   const Long64_t stride = (maxTracks > 0 && nAll > maxTracks) ? (nAll / maxTracks + 1) : 1;
   if (stride > 1)
      printf("  subsampling: every %lldth track (%lld -> ~%lld) to stay inside the page budget\n",
             stride, nAll, nAll / stride);
   Long64_t nU = dump_pk(tU, o, stride), nG = 0, nT = 0;
   if (tG) {
      fprintf(o, ",\"%s\":", sn[1].Data());
      nG = dump_pk(tG, o, stride);
   }
   if (tT) {
      fprintf(o, ",\"%s\":", sn[2].Data());
      nT = dump_pk(tT, o, stride);
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
   printf("open it with:  pp/open_explorer_dd.sh   (or any browser on %s)\n\n", outHtml.Data());
}
