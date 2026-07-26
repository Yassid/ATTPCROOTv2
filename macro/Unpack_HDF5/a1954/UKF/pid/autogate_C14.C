/// @file autogate_C14.C
/// @brief Build a PID gate polygon AUTOMATICALLY from the ridge of the proton band in the
///        (sqrt(dEdx), Brho) plane, so no interactive X11 drawing is needed (WSL has no
///        working ROOT GUI here). Reads pid_C14_data.csv written by dump_pid_C14.C and
///        writes an AtParticleID JSON gate readable by gate_events_C14.C.
///
/// Method: slice in Brho; in each slice take the sqrt(dEdx) projection, find the modal bin
/// (the band ridge), and take the half-maximum width around it. The left/right edges are
/// smoothed over neighbouring slices and widened by `widen` (in units of the local half
/// width) so the gate wraps the full band rather than only its core.
///
///   root -b -q 'pid/autogate_C14.C("proton_14C",1,1)'
///   root -b -q 'pid/autogate_C14.C("deuteron_14C",1,2,0.06,1.0,1.6,15,"","",2.0)'   // shifted band
void autogate_C14(TString name = "proton_14C", int Z = 1, int A = 1, double brLo = 0.06, double brHi = 1.05,
                  double widen = 1.6, int minEntries = 12, TString csvFile = "", TString outJson = "",
                  double ridgeScale = 1.0, double sdMax = 60)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   TString dir = gSystem->DirName(__FILE__);
   if (csvFile.IsNull())
      csvFile = dir + "/pid_C14_data.csv";
   if (outJson.IsNull())
      outJson = dir + "/" + name + ".json";
   if (gSystem->AccessPathName(csvFile)) {
      printf("\033[1;31mERROR: %s not found (run dump_pid_C14.C first)\033[0m\n", csvFile.Data());
      return;
   }

   TTree *t = new TTree("pid", "pid");
   t->ReadFile(csvFile, "", ',');
   Float_t sd, br;
   t->SetBranchAddress("sqrtdedx", &sd);
   t->SetBranchAddress("brho", &br);

   const int nSl = 40; // Brho slices, log-spaced (the band is steep at low Brho)
   std::vector<std::vector<float>> slice(nSl);
   const double lLo = std::log(brLo), lHi = std::log(brHi);
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (br <= brLo || br >= brHi || sd <= 0 || sd > sdMax)
         continue;
      int k = (int)(nSl * (std::log(br) - lLo) / (lHi - lLo));
      if (k >= 0 && k < nSl)
         slice[k].push_back(sd);
   }

   // ridge + half-width per slice
   std::vector<double> brC(nSl, 0), xLo(nSl, 0), xHi(nSl, 0);
   std::vector<int> ok(nSl, 0);
   for (int k = 0; k < nSl; ++k) {
      brC[k] = std::exp(lLo + (lHi - lLo) * (k + 0.5) / nSl);
      if ((int)slice[k].size() < minEntries)
         continue;
      auto &v = slice[k];
      std::sort(v.begin(), v.end());
      // mode via the densest window holding 40% of the slice
      size_t w = std::max<size_t>(3, v.size() * 0.4);
      size_t best = 0;
      double bestW = 1e30;
      for (size_t i = 0; i + w <= v.size(); ++i)
         if (v[i + w - 1] - v[i] < bestW) {
            bestW = v[i + w - 1] - v[i];
            best = i;
         }
      double med = v[best + w / 2];               // ridge position
      double hw = std::max(0.5 * bestW, 0.08 * med); // half width (floor: 8% of the ridge)
      med *= ridgeScale;
      xLo[k] = std::max(0.2, med - widen * hw);
      xHi[k] = med + widen * hw;
      ok[k] = 1;
   }

   // smooth the edges over +-1 slice (only over filled slices)
   std::vector<double> sLo = xLo, sHi = xHi;
   for (int k = 0; k < nSl; ++k) {
      if (!ok[k])
         continue;
      double aL = xLo[k], aH = xHi[k];
      int n = 1;
      for (int d = -1; d <= 1; d += 2)
         if (k + d >= 0 && k + d < nSl && ok[k + d]) {
            aL += xLo[k + d];
            aH += xHi[k + d];
            ++n;
         }
      sLo[k] = aL / n;
      sHi[k] = aH / n;
   }

   // polygon: up the LEFT edge (low Brho -> high Brho), back down the RIGHT edge
   std::vector<std::pair<double, double>> poly;
   for (int k = 0; k < nSl; ++k)
      if (ok[k])
         poly.emplace_back(sLo[k], brC[k]);
   for (int k = nSl - 1; k >= 0; --k)
      if (ok[k])
         poly.emplace_back(sHi[k], brC[k]);
   if (poly.size() < 6) {
      printf("\033[1;31mERROR: only %zu ridge points -- not enough statistics\033[0m\n", poly.size());
      return;
   }
   poly.push_back(poly.front());

   std::ofstream js(outJson.Data());
   js << "{\n  \"name\": \"" << name << "\",\n  \"xaxis\": \"sqrtdEdx\",\n  \"yaxis\": \"brho\",\n";
   js << "  \"Z\": " << Z << ",\n  \"A\": " << A << ",\n  \"vertices\": [\n";
   for (size_t i = 0; i < poly.size(); ++i)
      js << "    [" << poly[i].first << ", " << poly[i].second << "]" << (i + 1 < poly.size() ? "," : "") << "\n";
   js << "  ]\n}\n";
   js.close();
   printf("wrote %s  (%zu vertices, %d filled slices)\n", outJson.Data(), poly.size(),
          (int)std::count(ok.begin(), ok.end(), 1));

   // control plot: band + polygon + fraction inside
   TH2F *h = new TH2F("h", TString::Format("a1954 14C PID with the auto %s gate;#sqrt{dEdx};B#rho [T m]", name.Data()),
                      300, 0, sdMax, 300, 0, 1.2 * brHi);
   auto *g = new TGraph();
   for (auto &p : poly)
      g->SetPoint(g->GetN(), p.first, p.second);
   g->SetLineColor(kRed);
   g->SetLineWidth(2);
   auto *cut = new TCutG("cut", poly.size());
   for (size_t i = 0; i < poly.size(); ++i)
      cut->SetPoint(i, poly[i].first, poly[i].second);
   long nIn = 0, nAll = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      h->Fill(sd, br);
      ++nAll;
      if (cut->IsInside(sd, br))
         ++nIn;
   }
   printf("tracks inside the gate: %ld / %ld (%.1f%%)\n", nIn, nAll, 100.0 * nIn / std::max(1L, nAll));
   TCanvas *c = new TCanvas("c", "autogate", 900, 700);
   c->SetLogz();
   h->Draw("colz");
   g->Draw("L same");
   TString png = dir + "/plots/autogate_" + name + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
