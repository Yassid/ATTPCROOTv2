/// @file csv_to_root.C
/// @brief Turn the browser explorer's "save data" CSV into a ROOT file.
///
/// This is the way to save WHAT YOU SEE. pp/save_panels_root.C rebuilds panels from the cache with
/// whatever parameters are passed to it, which is not the same thing: the page's Ebeam, cuts,
/// binning and theta correction are live, and only the page knows their current values. saveCsv()
/// dumps the panels exactly as drawn, so converting that CSV preserves the on-screen state
/// including any correction applied.
///
/// Workflow:
///   1. set the viewer up the way you want it
///   2. click "save data"  -> <stem>_hist.csv lands in the browser's download folder
///   3. root -b -q 'pp/csv_to_root.C("/mnt/c/Users/Yassid/Downloads/xxx_hist.csv","panels.root")'
///
/// Each panel becomes a TH1D or TH2D named after its panel id, with the axis ranges the page used
/// (carried in the '# axes ...' line). The leading '#' comment block is preserved verbatim in a
/// TNamed called "provenance", so the beam energy and fitter that produced the numbers travel with
/// them. A TCanvas 'cPanels' draws them all.
///
/// Requires an explorer built on or after the template change that adds the '# axes' / '# title'
/// lines; an older CSV has no axis ranges and is rejected rather than guessed at.

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

static std::string trim_(std::string s)
{
   size_t a = s.find_first_not_of(" \t\r\n");
   if (a == std::string::npos) return "";
   size_t b = s.find_last_not_of(" \t\r\n");
   return s.substr(a, b - a + 1);
}

/// pull "key=value" out of a comma-separated '# axes ...' line
static double axVal_(const std::string &line, const std::string &key, bool &ok)
{
   ok = false;
   size_t p = line.find(key + "=");
   if (p == std::string::npos) return 0;
   p += key.size() + 1;
   size_t e = line.find_first_of(",", p);
   ok = true;
   return atof(line.substr(p, e == std::string::npos ? std::string::npos : e - p).c_str());
}

void csv_to_root(TString csvPath, TString outFile = "")
{
   gStyle->SetOptStat(0);
   std::ifstream in(csvPath.Data());
   if (!in) { printf("cannot open %s\n", csvPath.Data()); return; }
   if (outFile.IsNull()) {
      outFile = csvPath;
      outFile.ReplaceAll(".csv", "");
      outFile += ".root";
   }

   std::vector<std::string> header;      // the leading '#' block: title, fitter, Ebeam, timestamp
   std::vector<TObject *> objs;
   std::string line, curId, curKind, curTitle;
   bool inPanel = false, sawAxes = false, sawData = false;
   int nb = 0, nx = 0, ny = 0;
   double lo = 0, hi = 0, xlo = 0, xhi = 0, ylo = 0, yhi = 0;
   TH1D *h1 = nullptr; TH2D *h2 = nullptr;
   int nPanel = 0, nRejected = 0;

   auto closePanel = [&]() {
      if (!inPanel) return;
      if (!sawAxes) {
         printf("  panel '%s': no '# axes' line -- SKIPPED (regenerate the explorer to get it)\n",
                curId.c_str());
         ++nRejected;
      } else if (h1) { objs.push_back(h1); ++nPanel; }
      else if (h2)   { objs.push_back(h2); ++nPanel; }
      h1 = nullptr; h2 = nullptr; inPanel = false; sawAxes = false; sawData = false; curTitle.clear();
   };

   while (std::getline(in, line)) {
      std::string s = trim_(line);
      if (s.empty()) continue;

      if (s[0] == '#') {
         if (s.find("# ---- panel ") == 0) {
            closePanel();
            std::string rest = s.substr(13);
            size_t sp = rest.find(" (");
            curId = rest.substr(0, sp);
            curKind = rest.substr(sp + 2, rest.find(')', sp) - sp - 2);
            inPanel = true;
            continue;
         }
         if (inPanel && s.find("# axes") == 0) {
            bool ok = false;
            if (curKind == "hist") {
               nb = (int)axVal_(s, "nb", ok); if (!ok) continue;
               lo = axVal_(s, "lo", ok); hi = axVal_(s, "hi", ok);
            } else {
               nx = (int)axVal_(s, "nx", ok); if (!ok) continue;
               xlo = axVal_(s, "xlo", ok); xhi = axVal_(s, "xhi", ok);
               ny = (int)axVal_(s, "ny", ok); ylo = axVal_(s, "ylo", ok); yhi = axVal_(s, "yhi", ok);
            }
            sawAxes = ok;
            continue;
         }
         if (inPanel && s.find("# title") == 0) { curTitle = trim_(s.substr(7)); continue; }
         if (!inPanel) header.push_back(s.substr(1));
         continue;
      }
      if (!inPanel || !sawAxes) continue;
      if (s.find("bin_centre") == 0 || s.find("x_centre") == 0) continue;   // column header

      if (!sawData) {   // first data row: create the object
         if (curKind == "hist")
            h1 = new TH1D(curId.c_str(), Form("%s;E_{x} [MeV];counts",
                                              curTitle.empty() ? curId.c_str() : curTitle.c_str()),
                          nb, lo, hi);
         else
            h2 = new TH2D(curId.c_str(), (curTitle.empty() ? curId : curTitle).c_str(),
                          nx, xlo, xhi, ny, ylo, yhi);
         sawData = true;
      }
      std::stringstream ss(s);
      std::string a, b, c;
      if (curKind == "hist") {
         std::getline(ss, a, ','); std::getline(ss, b, ',');
         if (h1) h1->SetBinContent(h1->GetXaxis()->FindBin(atof(a.c_str())), atof(b.c_str()));
      } else {
         std::getline(ss, a, ','); std::getline(ss, b, ','); std::getline(ss, c, ',');
         if (h2) h2->SetBinContent(h2->GetXaxis()->FindBin(atof(a.c_str())),
                                   h2->GetYaxis()->FindBin(atof(b.c_str())), atof(c.c_str()));
      }
   }
   closePanel();

   if (objs.empty()) {
      printf("no panels recovered from %s%s\n", csvPath.Data(),
             nRejected ? " (all lacked '# axes' -- rebuild the explorer first)" : "");
      return;
   }

   TFile out(outFile, "RECREATE");
   std::string prov;
   for (auto &l : header) prov += trim_(l) + " | ";
   auto *pv = new TNamed("provenance", prov.c_str());

   int nc = (int)objs.size();
   int cols = nc > 2 ? 2 : nc, rows = (nc + cols - 1)/cols;
   auto *c = new TCanvas("cPanels", "explorer panels", 700*cols, 500*rows);
   c->Divide(cols, rows);
   for (int i = 0; i < nc; ++i) {
      c->cd(i+1);
      if (objs[i]->InheritsFrom("TH2")) { gPad->SetLogz(); objs[i]->Draw("colz"); }
      else objs[i]->Draw("hist");
   }
   for (auto *o : objs) o->Write();
   pv->Write();
   c->Write();
   out.Close();
   printf("wrote %s  (%d panels%s)\n  provenance: %s\n", outFile.Data(), nPanel,
          nRejected ? Form(", %d skipped", nRejected) : "", prov.c_str());
}
