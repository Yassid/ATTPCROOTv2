/// @file plot_padgain_radial.C
/// @brief Radial profile (ring medians, 10 mm in rho) of pad-gain maps from pad_gain_map_dp.C, several maps
///        on one canvas -- e.g. early / middle / late data periods and the simulation (flat reference).
///        maps = "label:file;label:file;..." (files need the x y columns, i.e. maps made after 2026-10-31).
///   root -l -b -q 'plot_padgain_radial.C("early:pid/padgain_d2_early.txt;sim:pid/padgain_sim_dp_gs.txt","plots/padgain_radial.png")'
void plot_padgain_radial(TString maps, TString png, int minHits = 50, int minPads = 20)
{
   gStyle->SetOptStat(0);
   auto *c = new TCanvas("c", "c", 1000, 700);
   gPad->SetGrid();
   auto *fr = gPad->DrawFrame(0, 0.6, 150, 1.2, "Relative pad gain vs distance from the beam axis;#rho [mm];ring median of the pad gain");
   (void)fr;
   auto *lg = new TLegend(0.55, 0.15, 0.89, 0.42);
   const int col[8] = {kBlue + 1, kGreen + 2, kRed + 1, kBlack, kGray + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
   int k = 0;
   for (auto *o : *maps.Tokenize(";")) {
      auto *f = ((TObjString *)o)->GetString().Tokenize(":");
      TString label = ((TObjString *)f->At(0))->GetString(), file = ((TObjString *)f->At(1))->GetString();
      std::map<int, std::vector<double>> ring;
      std::ifstream in(file.Data());
      std::string line;
      while (std::getline(in, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream ss(line);
         int p; double g; long n; double x, y, gl;
         if (!(ss >> p >> g >> n >> x >> y >> gl) || n < minHits) continue;
         ring[int(std::hypot(x, y) / 10)].push_back(g);
      }
      auto *gr = new TGraph;
      for (auto &[r, v] : ring) {
         if ((int)v.size() < minPads || r * 10 + 5 > 148) continue;
         std::sort(v.begin(), v.end());
         gr->AddPoint(r * 10 + 5, v[v.size() / 2]);
      }
      gr->SetLineColor(col[k % 8]); gr->SetMarkerColor(col[k % 8]); gr->SetMarkerStyle(20 + k); gr->SetLineWidth(2);
      gr->Draw("LP");
      lg->AddEntry(gr, label, "lp");
      printf("%-30s %d rings\n", label.Data(), gr->GetN());
      ++k;
   }
   lg->Draw();
   c->SaveAs(png);
}
