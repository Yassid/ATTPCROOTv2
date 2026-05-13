/// @file plot_scan_p.C
/// @brief Plot σ_p/p vs pion momentum from data/scan_p_results.csv.

void plot_scan_p(const char *csv = "data/scan_p_results.csv",
                 const char *outPng = "data/scan_p.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.04, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   std::ifstream f(csv);
   if (!f) { std::cerr << "no " << csv << "\n"; return; }
   std::string line;
   std::getline(f, line); // header
   std::vector<double> p, sig;
   while (std::getline(f, line)) {
      double p_, ke_, nf, nt, pmc, bias, s;
      char c;
      std::stringstream ss(line);
      ss >> p_ >> c >> ke_ >> c >> nf >> c >> nt >> c >> pmc >> c >> bias >> c >> s;
      p.push_back(p_);
      sig.push_back(s * 100.); // %
   }

   auto *g = new TGraph(p.size(), p.data(), sig.data());
   g->SetMarkerStyle(20); g->SetMarkerSize(1.4); g->SetMarkerColor(kBlue + 1);
   g->SetLineColor(kBlue + 1); g->SetLineWidth(2);
   g->SetTitle("HYDRA Prototype: #sigma_{p}/p vs p, B = 2 T, 500 evts/pt;p_{MC} (MeV/c);#sigma_{p}/p (%)");

   auto *c = new TCanvas("c", "", 900, 600);
   c->SetGrid();
   g->Draw("APL");
   g->GetYaxis()->SetRangeUser(0., 10.);
   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
}
