/// @file dt_sim_dists.C
/// @brief What the 16C(d,t)15C generator is putting into the chamber, next to what the data holds.
///
/// Truth only -- nothing here is digitised or reconstructed, so this says what was GENERATED, not
/// what would be measured. The data curve is reconstructed and includes every Ex (peaks plus
/// continuum), while the simulation is a single state, so the two are expected to agree in TREND
/// and not point by point. What would be a real disagreement: a different slope of KE against
/// theta_cm, or tritons landing outside the theta_lab window the analysis cuts on.
///
///   root -l 'dt_sim_dists.C("/mnt/f/a1975_C16_dt_sim/gs_s3001_sim.root",0.0)'

void dt_sim_dists(TString simFile = "./data/dt_conv2.root", Double_t resEx = 0.0, Double_t Ebeam = 184.17,
                  TString dataFile = "/home/yassid/a1975_C16dt_analysis/data/dt_spyral_solver.root")
{
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double u = 931.49401;
   const double m_C16 = 16.0147 * u, m_d = 2.0141017781 * u, m_t = 3.0160492779 * u;
   const double E1 = Ebeam + m_C16, p1 = std::sqrt(E1 * E1 - m_C16 * m_C16), Etot = E1 + m_d;
   const double beta = p1 / Etot, gam = 1.0 / std::sqrt(1 - beta * beta);

   TFile *fs = TFile::Open(simFile);
   TTree *ts = fs && !fs->IsZombie() ? (TTree *)fs->Get("cbmsim") : nullptr;
   if (!ts) { printf("\n  cannot open %s\n\n", simFile.Data()); return; }
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);

   auto *hKEth = new TH2D("hKEth", "TRUTH: triton locus;#theta_{lab} [deg];KE_{t} [MeV]", 70, 0, 70, 80, 0, 40);
   auto *hKEcm = new TH2D("hKEcm", "TRUTH: KE vs #theta_{cm};#theta_{cm} [deg];KE_{t} [MeV]", 70, 0, 70, 80, 0, 40);
   auto *hCM = new TH1D("hCM", "generated #theta_{cm};#theta_{cm} [deg];tritons", 35, 0, 70);
   auto *hZ = new TH1D("hZ", "reaction vertex;z [mm];reactions", 40, 0, 1000);
   auto *hTh = new TH1D("hTh", "#theta_{lab};#theta_{lab} [deg];tritons", 70, 0, 70);
   std::vector<std::vector<double>> Ks(7);

   for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
      ts->GetEntry(i);
      if (!mc) continue;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *tr = (AtMCTrack *)mc->At(k);
         if (!tr || tr->GetPdgCode() != 1000010030 || tr->GetMotherId() != -1) continue;
         double px = tr->GetPx() * 1000, py = tr->GetPy() * 1000, pz = tr->GetPz() * 1000;
         double p = std::sqrt(px * px + py * py + pz * pz);
         if (p <= 0) continue;
         double E3 = std::sqrt(p * p + m_t * m_t), ke = E3 - m_t, th = std::acos(pz / p);
         double t3 = std::atan2(p * std::sin(th), gam * (p * std::cos(th) - beta * E3)) * TMath::RadToDeg();
         double tcm = 180.0 - t3; // residual convention, the one the data analysis uses
         hKEth->Fill(th * TMath::RadToDeg(), ke);
         hKEcm->Fill(tcm, ke);
         hCM->Fill(tcm);
         hTh->Fill(th * TMath::RadToDeg());
         hZ->Fill(tr->GetStartZ() * 10.0);
         int b = (int)(tcm / 10);
         if (b >= 0 && b < 7) Ks[b].push_back(ke);
         break;
      }
   }

   // the data, for the trend comparison
   auto *gData = new TGraph();
   auto *gSim = new TGraph();
   TFile *fd = TFile::Open(dataFile);
   TTree *td = fd && !fd->IsZombie() ? (TTree *)fd->Get("dt") : nullptr;
   if (td) {
      double ke, thl, tcm;
      td->SetBranchAddress("ke", &ke);
      td->SetBranchAddress("thl", &thl);
      td->SetBranchAddress("tcm", &tcm);
      std::vector<std::vector<double>> Kd(7);
      for (Long64_t i = 0; i < td->GetEntries(); ++i) {
         td->GetEntry(i);
         if (!std::isfinite(ke) || thl <= 8 || thl >= 56) continue;
         int b = (int)(tcm / 10);
         if (b >= 0 && b < 7) Kd[b].push_back(ke);
      }
      for (int b = 0; b < 7; ++b) {
         if (Kd[b].size() > 10) {
            std::sort(Kd[b].begin(), Kd[b].end());
            gData->SetPoint(gData->GetN(), 10 * b + 5, Kd[b][Kd[b].size() / 2]);
         }
         if (Ks[b].size() > 5) {
            std::sort(Ks[b].begin(), Ks[b].end());
            gSim->SetPoint(gSim->GetN(), 10 * b + 5, Ks[b][Ks[b].size() / 2]);
         }
      }
   }

   auto *c = new TCanvas("c_dists", "16C(d,t)15C generated distributions", 1500, 900);
   c->Divide(3, 2, 0.004, 0.004);
   c->cd(1); hKEth->Draw("colz");
   auto *lcut = new TLine(56, 0, 56, 40);
   lcut->SetLineColor(kRed + 1); lcut->SetLineStyle(2); lcut->SetLineWidth(2); lcut->Draw();
   auto *tcut = new TLatex(40, 36, "#color[2]{analysis cut}");
   tcut->SetTextSize(0.04); tcut->Draw();
   c->cd(2); hKEcm->Draw("colz");
   c->cd(3);
   auto *fr = gPad->DrawFrame(0, 0, 70, 35);
   fr->SetTitle("KE trend: sim truth vs data;#theta_{cm} [deg];median KE_{t} [MeV]");
   gSim->SetMarkerStyle(20); gSim->SetMarkerColor(kAzure + 2); gSim->SetLineColor(kAzure + 2); gSim->SetLineWidth(2);
   gData->SetMarkerStyle(24); gData->SetMarkerColor(kBlack); gData->SetLineColor(kBlack); gData->SetLineWidth(2);
   gSim->Draw("PL same"); gData->Draw("PL same");
   auto *lg = new TLegend(0.15, 0.68, 0.60, 0.86);
   lg->SetFillStyle(0);
   lg->AddEntry(gSim, Form("simulation truth, Ex = %.3f", resEx), "lp");
   lg->AddEntry(gData, "data (all Ex, reconstructed)", "lp");
   lg->Draw();
   c->cd(4); hCM->Draw("hist");
   c->cd(5); hTh->Draw("hist");
   auto *l56 = new TLine(56, 0, 56, hTh->GetMaximum() * 1.05);
   l56->SetLineColor(kRed + 1); l56->SetLineStyle(2); l56->Draw();
   c->cd(6); hZ->Draw("hist");
   c->SaveAs("./data/dt_sim_dists.png");
   printf("\n  wrote data/dt_sim_dists.png -- canvas c_dists is interactive\n\n");
}
