/// @file mkpid_C15d.C
/// @brief Merge the per-run gain-matched PID caches into one plane and draw it.
///
///   root -b -q 'mkpid_C15d.C()'                                  // whatever is cached so far
///   root -b -q 'mkpid_C15d.C("/home/yassid/C15d_reco/", "plots/", 60, 0, 4, 200, 0, 2.5)'
///
/// Reads every <run>_pid.root in inDir, fills sqrt(dE/dx) vs Brho, writes
/// plots/pid_C15d.png and plots/pid_C15d.root (the TH2 plus the merged TChain's
/// selection counts). Cheap: the caches are a couple of MB per run, so re-binning
/// and re-cutting is instant.
///
/// GATES OVERLAY: pass a comma-separated list of Cut2D JSON files (as written by
/// draw_gate_C15d.C) and each is drawn on top with its in-gate count.
///
/// GAIN: the caches hold RAW dE/dx, so pass gainTable to see the matched plane. Without it the
/// axis label says so explicitly.
///
/// The `valid` flag is stored, not pre-filtered, so the rejected population stays countable.
/// Default cuts are deliberately loose: quality cuts belong downstream of seeing the plane,
/// not baked into it.

/// Default axis ranges cover the whole populated plane (sqrt_dEdx up to ~85, brho up to 2.5) so
/// nothing is clipped out of view by the plot itself -- a clipped axis looks exactly like a band
/// that ends.
#include "gain_C15d.h"

void mkpid_C15d(TString inDir = "/home/yassid/C15d_reco/", TString outDir = "plots/", Int_t nbx = 340,
                Double_t xlo = 0.0, Double_t xhi = 85.0, Int_t nby = 300, Double_t ylo = 0.0,
                Double_t yhi = 2.5, TString gates = "", Int_t minClusters = 0, Double_t maxVtxR = -1.0,
                /// Set once measure_gain_C15d.C has produced a table AND the caches were built with
                /// gain matching on -- it only controls the axis label, but a mislabelled plane is
                /// how an unmatched plane ends up being trusted as a matched one.
                TString gainTable = "",
                /// Restrict to a run range. Runs >=106 are a HYDROGEN target -- plotting them on
                /// the same plane overlays two different gases and doubles every band.
                Int_t runMin = 17, Int_t runMax = 103)
{
   gSystem->Load("libAtTools.so");
   gSystem->mkdir(outDir, kTRUE);

   TChain ch("pid");
   TSystemDirectory dir(inDir, inDir);
   TList *files = dir.GetListOfFiles();
   if (files == nullptr) {
      std::cout << "\033[1;31mERROR: cannot list " << inDir << "\033[0m\n";
      return;
   }
   Int_t nRuns = 0;
   TIter next(files);
   while (auto *o = dynamic_cast<TSystemFile *>(next())) {
      TString name = o->GetName();
      if (o->IsDirectory() || !name.EndsWith("_pid.root"))
         continue;
      TString d = name;
      d.ReplaceAll("run_", "");
      d.ReplaceAll("_pid.root", "");
      const Int_t rn = d.Atoi();
      if (rn < runMin || rn > runMax)
         continue;
      ch.Add(inDir + name);
      ++nRuns;
   }
   if (nRuns == 0) {
      std::cout << "\033[1;31mERROR: no *_pid.root in " << inDir << " -- run reco_batch.sh first.\033[0m\n";
      return;
   }

   TString sel = "valid==1";
   if (minClusters > 0)
      sel += TString::Format(" && nClusters>=%d", minClusters);
   if (maxVtxR > 0)
      sel += TString::Format(" && vtxR<%g", maxVtxR);

   const Long64_t nAll = ch.GetEntries();
   const Long64_t nSel = ch.GetEntries(sel);
   std::cout << "\033[1;33m=== C15d PID plane ===\033[0m\n"
             << "  runs        : " << nRuns << "  (range " << runMin << "-" << runMax << ")\n"
             << "  tracks      : " << nAll << "\n"
             << "  selected    : " << nSel << "  (" << sel << ")  = " << (nAll ? 100.0 * nSel / nAll : 0.)
             << "%\n";

   // The reco persists RAW dE/dx (unpackReco_C15d.C has gainMatch off by default), so say so on
   // the axis. A plane labelled "gain matched" that is not is the kind of thing that survives
   // into a talk.
   auto gainMap = LoadGainTable_C15d(gainTable);
   TString xlabel = gainMap.empty() ? "#sqrt{dE/dx}  (raw, per-run gain NOT matched)"
                                    : "#sqrt{dE/dx}  (gain matched)";
   auto *h = new TH2D("hpid", ";" + xlabel + ";B#rho  (T m)", nbx, xlo, xhi, nby, ylo, yhi);

   // Explicit loop rather than TChain::Draw: the gain factor is per RUN, so it cannot be written
   // as a single TTree expression over a chain of runs.
   {
      Int_t b_run, b_valid, b_nClusters;
      Double_t b_sqrt, b_brho, b_vtxR;
      ch.SetBranchAddress("run", &b_run);
      ch.SetBranchAddress("valid", &b_valid);
      ch.SetBranchAddress("nClusters", &b_nClusters);
      ch.SetBranchAddress("sqrtdEdx", &b_sqrt);
      ch.SetBranchAddress("brho", &b_brho);
      ch.SetBranchAddress("vtxR", &b_vtxR);
      Long64_t nMissing = 0;
      std::set<int> missingRuns;
      for (Long64_t e = 0; e < ch.GetEntries(); ++e) {
         ch.GetEntry(e);
         if (b_valid != 1)
            continue;
         if (minClusters > 0 && b_nClusters < minClusters)
            continue;
         if (maxVtxR > 0 && b_vtxR >= maxVtxR)
            continue;
         bool missing = false;
         const double f = GainFactor_C15d(gainMap, b_run, missing);
         if (missing) {
            ++nMissing;
            missingRuns.insert(b_run);
         }
         h->Fill(b_sqrt * std::sqrt(f), b_brho);
      }
      if (nMissing > 0) {
         std::cout << "\033[1;31m  WARNING: " << nMissing << " tracks from " << missingRuns.size()
                   << " run(s) are NOT in the gain table and went in UNMATCHED:";
         for (int r : missingRuns)
            std::cout << " " << r;
         std::cout << "\033[0m\n";
      }
      ch.ResetBranchAddresses();
   }

   TCanvas *c = new TCanvas("cpid", "C15d PID", 1100, 850);
   c->SetLogz();
   c->SetRightMargin(0.13);
   h->Draw("colz");

   // Overlay any gates and report how many tracks each one holds.
   std::vector<TGraph *> keep;
   if (gates.Length()) {
      TObjArray *toks = gates.Tokenize(",");
      const int cols[] = {kRed + 1, kBlack, kMagenta + 1, kGreen + 2, kOrange + 7, kCyan + 2};
      for (int i = 0; i < toks->GetEntries(); ++i) {
         TString path = ((TObjString *)toks->At(i))->GetString().Strip(TString::kBoth);
         auto cut = AtTools::AtCut2D::LoadJSON(path.Data());
         if (!cut.IsValid()) {
            std::cout << "\033[1;31m  gate " << path << " : INVALID / not found\033[0m\n";
            continue;
         }
         const auto &v = cut.GetVertices();
         auto *g = new TGraph(static_cast<Int_t>(v.size()) + 1);
         for (size_t k = 0; k < v.size(); ++k)
            g->SetPoint(static_cast<Int_t>(k), v[k].first, v[k].second);
         g->SetPoint(static_cast<Int_t>(v.size()), v[0].first, v[0].second); // close it
         g->SetLineColor(cols[i % 6]);
         g->SetLineWidth(3);
         g->Draw("L same");
         keep.push_back(g);

         // Count in-gate with the polygon test itself, not a bounding box.
         Long64_t in = 0;
         Double_t x, y;
         Int_t valid;
         TChain cc("pid");
         TIter n2(files);
         while (auto *o2 = dynamic_cast<TSystemFile *>(n2())) {
            TString nm = o2->GetName();
            if (!o2->IsDirectory() && nm.EndsWith("_pid.root")) {
               TString d2 = nm;
               d2.ReplaceAll("run_", "");
               d2.ReplaceAll("_pid.root", "");
               const Int_t r2 = d2.Atoi();
               if (r2 >= runMin && r2 <= runMax)
                  cc.Add(inDir + nm);
            }
         }
         cc.SetBranchAddress("sqrtdEdx", &x);
         cc.SetBranchAddress("brho", &y);
         cc.SetBranchAddress("valid", &valid);
         for (Long64_t e = 0; e < cc.GetEntries(); ++e) {
            cc.GetEntry(e);
            if (valid == 1 && cut.IsInside(x, y))
               ++in;
         }
         std::cout << "  gate " << cut.GetName() << " (" << cut.GetXAxis() << " x " << cut.GetYAxis()
                   << ", " << v.size() << " vertices) : " << in << " tracks = "
                   << (nSel ? 100.0 * in / nSel : 0.) << "% of selected\n";
      }
      delete toks;
   }

   c->SaveAs(outDir + "pid_C15d.png");
   TFile fo(outDir + "pid_C15d.root", "RECREATE");
   h->Write();
   fo.Close();
   std::cout << "\033[1;32mwrote\033[0m " << outDir << "pid_C15d.png and .root\n";
}
