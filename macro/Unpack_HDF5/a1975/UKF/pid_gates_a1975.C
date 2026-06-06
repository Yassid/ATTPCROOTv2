/// @file pid_gates_a1975.C
/// @brief Inspect proton/deuteron separation: high-res IC+quality-gated PID plane,
/// optional AtParticleID gates overlaid, and the theta-vs-brho kinematics of the
/// tracks inside each gate (to expose deuteron leakage into a proton gate).
///
/// Reads the cached observables from pid_ic_a1975 (combined_pidobs.root).
///
/// Run:
///   root -b -q 'pid_gates_a1975.C("combined")'                         // just the bands
///   root -b -q 'pid_gates_a1975.C("combined","proton_pid.json,deuteron_pid.json")'

void pid_gates_a1975(TString cacheTag = "combined", TString gateFiles = "", Double_t icMin = 950, Double_t icMax = 1350,
                     Int_t minClusters = 15, Double_t maxVertexR = 40, Double_t polarMin = 10, Double_t polarMax = 170,
                     Double_t dedxMax = 2000, Double_t brMax = 0.6)
{
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   TFile *fo = TFile::Open(cacheTag + "_pidobs.root");
   if (!fo || fo->IsZombie()) {
      printf("\033[1;31mno cache %s_pidobs.root\033[0m\n", cacheTag.Data());
      return;
   }
   TNtuple *nt = (TNtuple *)fo->Get("pidobs");
   float dedx, brho, ncl, ic, polar, vtxr, vtxz;
   nt->SetBranchAddress("dedx", &dedx);
   nt->SetBranchAddress("brho", &brho);
   nt->SetBranchAddress("ncl", &ncl);
   nt->SetBranchAddress("ic", &ic);
   nt->SetBranchAddress("polar", &polar);
   nt->SetBranchAddress("vtxr", &vtxr);
   nt->SetBranchAddress("vtxz", &vtxz);

   // Load gates
   std::vector<AtTools::AtParticleID> gates;
   std::vector<TString> gnames;
   if (gateFiles.Length() > 0) {
      TObjArray *gt = gateFiles.Tokenize(",");
      for (int i = 0; i < gt->GetEntries(); ++i) {
         TString gf = ((TObjString *)gt->At(i))->GetString().Strip(TString::kBoth);
         auto pid = AtTools::AtParticleID::LoadJSON(gf.Data());
         if (pid.IsValid()) {
            gates.push_back(pid);
            gnames.push_back(pid.GetName().c_str());
         }
      }
      delete gt;
   }

   TH2F *hpid = new TH2F("hpid", "IC+quality PID;dEdx [counts];B#rho [T m]", 500, 0, dedxMax, 500, 0, brMax);
   std::vector<TH2F *> hkin;
   std::vector<long> gcount(gates.size(), 0);
   for (size_t g = 0; g < gates.size(); ++g)
      hkin.push_back(new TH2F(Form("hkin%zu", g), Form("Kinematics gated: %s;#theta [deg];B#rho [T m]", gnames[g].Data()),
                              180, 0, 180, 400, 0, brMax));

   long clean = 0;
   for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
      nt->GetEntry(i);
      if (ncl < minClusters || vtxr > maxVertexR || polar < polarMin || polar > polarMax)
         continue;
      if (ic < icMin || ic > icMax)
         continue;
      ++clean;
      hpid->Fill(dedx, brho);
      for (size_t g = 0; g < gates.size(); ++g)
         if (gates[g].IsInside(dedx, brho)) {
            hkin[g]->Fill(polar, brho);
            gcount[g]++;
         }
   }
   printf("\033[1;32mclean (IC+quality) tracks: %ld\033[0m\n", clean);
   for (size_t g = 0; g < gates.size(); ++g)
      printf("  gate '%s': %ld tracks inside\n", gnames[g].Data(), gcount[g]);

   // Save each panel as its own PNG (<=900px) so they're individually viewable.
   {
      TCanvas *c0 = new TCanvas("c0", "pid", 900, 720);
      c0->SetLogz();
      hpid->Draw("colz");
      for (size_t g = 0; g < gates.size(); ++g) {
         const auto &v = gates[g].GetCut().GetVertices();
         auto *pl = new TPolyLine(v.size() + 1);
         for (size_t k = 0; k < v.size(); ++k)
            pl->SetPoint(k, v[k].first, v[k].second);
         pl->SetPoint(v.size(), v[0].first, v[0].second);
         pl->SetLineColor(g == 0 ? kRed : kMagenta);
         pl->SetLineWidth(3);
         pl->Draw();
      }
      c0->SaveAs(cacheTag + "_pid_plane.png");
   }
   for (size_t g = 0; g < gates.size(); ++g) {
      TCanvas *ck = new TCanvas(Form("ck%zu", g), "kin", 900, 700);
      ck->SetLogz();
      hkin[g]->Draw("colz");
      ck->SaveAs(cacheTag + "_kin_" + gnames[g] + ".png");
   }
   printf("Saved %s_pid_plane.png and per-gate kinematics PNGs\n", cacheTag.Data());
}
