/// @file make_pid_explorer_html.C
/// @brief Bake the C15d PID plane into a single self-contained HTML explorer.
///
///   root -b -q 'make_pid_explorer_html.C()'
///   root -b -q 'make_pid_explorer_html.C("/home/yassid/C15d_reco/", "/home/yassid/C15d_pid.html")'
///
/// Same idea as the a2091 kinematics explorer (pp/make_explorer_html.C): one HTML file with the
/// data baked in, so it runs in a browser and needs no X11. Different content, because at this
/// stage there are no fits to plot -- what exists is the PID plane, and what is blocking the
/// analysis is drawing gates on it.
///
/// So this page DRAWS GATES. Click round a band, close the polygon, and it writes spyral_utils
/// Cut2D JSON with Z/A -- the same format AtCut2D reads and draw_gate_C15d.C writes. Save it under
/// gates/ and apply_gate_C15d.C does the exact per-track selection.
///
/// ★ RAW dE/dx IS BAKED, AND THE GAIN TABLE ALONGSIDE IT. The page applies the factor itself, so
/// the matched/raw toggle is live and the underlying measurement is always present. That also
/// means the page cannot silently show a matched plane it cannot reproduce: with no table it says
/// so, in red, because a gate drawn on a raw plane does not transfer (the per-run drift is ~29 %
/// over this run set).
///
/// Values are stored as scaled integers -- sqrt(dE/dx) x100, Brho x1000, vertex R x1 -- which is
/// what keeps 400k tracks to a few MB of JSON instead of a few tens.

#include "gain_C15d.h"

void make_pid_explorer_html(TString inDir = "/home/yassid/C15d_reco/", TString outHtml = "",
                            TString gainTable = "gainmatch_C15d.csv", Int_t runMin = 17,
                            Int_t runMax = 103, Double_t xmax = 60.0, Double_t ymax = 2.0,
                            Int_t bins = 300, TString tag = "15C + d  --  PID plane",
                            /// Keep every Nth track. 1 = all. Only raise it if the page is too
                            /// heavy for the browser; the counts scale but the bands do not move.
                            Int_t stride = 1)
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (outHtml.IsNull())
      outHtml = TString(gSystem->Getenv("HOME")) + "/C15d_pid_explorer.html";
   TString tplPath = here + "/pid_explorer_template.html";

   std::ifstream tf(tplPath.Data());
   if (!tf) {
      std::cerr << "cannot read template " << tplPath << "\n";
      return;
   }
   std::string tpl((std::istreambuf_iterator<char>(tf)), std::istreambuf_iterator<char>());
   tf.close();

   // ---- gain table (resolved next to this macro if given as a bare name) -------------------
   TString gt = gainTable;
   if (gt.Length() && !gt.BeginsWith("/") && gSystem->AccessPathName(gt))
      gt = here + "/" + gainTable;
   auto gain = LoadGainTable_C15d(gt);

   // ---- collect the per-run caches ---------------------------------------------------------
   TChain ch("pid");
   TSystemDirectory dir(inDir, inDir);
   TList *files = dir.GetListOfFiles();
   if (files == nullptr) {
      std::cerr << "cannot list " << inDir << "\n";
      return;
   }
   Int_t nRuns = 0;
   TIter next(files);
   while (auto *o = dynamic_cast<TSystemFile *>(next())) {
      TString n = o->GetName();
      if (o->IsDirectory() || !n.EndsWith("_pid.root"))
         continue;
      TString d = n;
      d.ReplaceAll("run_", "");
      d.ReplaceAll("_pid.root", "");
      const Int_t rn = d.Atoi();
      if (rn < runMin || rn > runMax)
         continue;
      ch.Add(inDir + n);
      ++nRuns;
   }
   if (nRuns == 0) {
      std::cerr << "no *_pid.root in " << inDir << " for runs " << runMin << "-" << runMax << "\n";
      return;
   }

   Int_t run, valid, nClusters;
   Double_t sqrtdEdx, brho, vtxR;
   ch.SetBranchAddress("run", &run);
   ch.SetBranchAddress("valid", &valid);
   ch.SetBranchAddress("nClusters", &nClusters);
   ch.SetBranchAddress("sqrtdEdx", &sqrtdEdx);
   ch.SetBranchAddress("brho", &brho);
   ch.SetBranchAddress("vtxR", &vtxR);

   // ---- data ------------------------------------------------------------------------------
   std::string sx, by, rr, nc, vr;
   sx.reserve(1 << 22);
   by.reserve(1 << 22);
   rr.reserve(1 << 21);
   nc.reserve(1 << 21);
   vr.reserve(1 << 21);
   char buf[64];
   Long64_t kept = 0, seen = 0, inView = 0;
   std::set<int> runsSeen;
   const Long64_t nEnt = ch.GetEntries();
   for (Long64_t i = 0; i < nEnt; ++i) {
      ch.GetEntry(i);
      if (valid != 1)
         continue;
      if (stride > 1 && (seen++ % stride))
         continue;
      const char *sep = kept ? "," : "";
      snprintf(buf, sizeof(buf), "%s%d", sep, (int)llround(sqrtdEdx * 100.0));
      sx += buf;
      snprintf(buf, sizeof(buf), "%s%d", sep, (int)llround(brho * 1000.0));
      by += buf;
      snprintf(buf, sizeof(buf), "%s%d", sep, run);
      rr += buf;
      snprintf(buf, sizeof(buf), "%s%d", sep, nClusters);
      nc += buf;
      snprintf(buf, sizeof(buf), "%s%d", sep, (int)llround(vtxR));
      vr += buf;
      runsSeen.insert(run);
      if (sqrtdEdx >= 0 && sqrtdEdx <= xmax && brho >= 0 && brho <= ymax)
         ++inView;
      ++kept;
   }
   if (kept == 0) {
      std::cerr << "no valid tracks in runs " << runMin << "-" << runMax << "\n";
      return;
   }

   // ---- gain map as JSON, only for runs actually present ------------------------------------
   std::string gj = "{";
   bool first = true;
   Int_t nMissing = 0;
   for (Int_t r = runMin; r <= runMax; ++r) {
      auto it = gain.find(r);
      if (it == gain.end())
         continue;
      snprintf(buf, sizeof(buf), "%s\"%d\":%.6f", first ? "" : ",", r, it->second);
      gj += buf;
      first = false;
   }
   gj += "}";

   std::string cfg = "{";
   snprintf(buf, sizeof(buf), "\"xmax\":%.1f,\"ymax\":%.2f,\"bins\":%d,", xmax, ymax, bins);
   cfg += buf;
   snprintf(buf, sizeof(buf), "\"rmin\":%d,\"rmax\":%d,", runMin, runMax);
   cfg += buf;
   cfg += "\"tag\":\"" + std::string(tag.Data()) + "\",";
   cfg += "\"gainName\":\"" + std::string(gSystem->BaseName(gt.Data())) + "\",";
   snprintf(buf, sizeof(buf), "\"sub\":\"%d runs%s\",", nRuns, stride > 1 ? Form(", 1/%d sampled", stride) : "");
   cfg += buf;
   cfg += "\"gain\":" + gj + "}";

   std::string data = "{\"sx\":[" + sx + "],\"by\":[" + by + "],\"run\":[" + rr + "],\"nc\":[" + nc +
                      "],\"vr\":[" + vr + "]}";

   const std::string kCfg = "__CFG__", kData = "__DATA__";
   size_t p = tpl.find(kCfg);
   if (p == std::string::npos) {
      std::cerr << "template has no __CFG__ placeholder\n";
      return;
   }
   tpl.replace(p, kCfg.size(), cfg);
   p = tpl.find(kData);
   if (p == std::string::npos) {
      std::cerr << "template has no __DATA__ placeholder\n";
      return;
   }
   tpl.replace(p, kData.size(), data);

   std::ofstream out(outHtml.Data());
   if (!out) {
      std::cerr << "cannot write " << outHtml << "\n";
      return;
   }
   out << tpl;
   out.close();

   std::cout << "\033[1;33m=== C15d PID explorer ===\033[0m\n"
             << "  runs      : " << nRuns << "  (" << runMin << "-" << runMax << ")\n"
             << "  tracks    : " << kept << " valid" << (stride > 1 ? Form(" (1/%d sampled)", stride) : "")
             << "\n"
             << "  gain      : " << (gain.empty() ? "NONE -- page will show a RAW plane" : gt.Data()) << "\n"
             << "  \033[1;32mwrote\033[0m " << outHtml << "  ("
             << (double)(tpl.size()) / (1024 * 1024) << " MB)\n"
             << "  open it in a browser; no X11 and no server needed.\n";
   // Compare against the runs actually PRESENT, not every integer in the range: most numbers in
   // 17-103 have no run on disk, and counting those as "missing a gain entry" cried wolf on 12
   // runs that do not exist.
   if (!gain.empty()) {
      std::vector<int> absent;
      for (int r : runsSeen)
         if (gain.find(r) == gain.end())
            absent.push_back(r);
      if (!absent.empty()) {
         std::cout << "\033[1;33m  WARNING: " << absent.size()
                   << " run(s) with data have NO gain entry and show unmatched:";
         for (int r : absent)
            std::cout << " " << r;
         std::cout << "\033[0m\n";
      }
   }
   // The default axes deliberately do not cover the full range: Brho carries a 1/sin(polar)
   // divergence for near-axis tracks, so a few percent sit at arbitrarily high rigidity. Say what
   // fraction is on screen rather than let a clipped tail pass as the whole sample.
   std::cout << "  in default view (x<=" << xmax << ", y<=" << ymax << ") : " << inView << " = "
             << (kept ? 100.0 * inView / kept : 0.) << "% -- the rest is the high-Brho tail,\n"
             << "  reachable with the y-max slider.\n";
}
