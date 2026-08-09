/// @file plot_states_C16pd.C
/// @brief The four generated 16C(p,d)15C states: kinematics, and what the beam energy loss does.
///
/// Everything here is MC TRUTH. There is no digitisation, no reconstruction and no detector
/// resolution in these numbers, so any width seen is generated physics plus the analysis
/// convention -- not a measurement error. That is the point: it isolates one effect.
///
/// The effect is the beam losing energy before it reacts. The vertex is uniform along the metre
/// of gas, the 16C loses about 28 MeV crossing all of it, and the excitation energy is
/// reconstructed at the NOMINAL 192 MeV, which is what the data analysis does. So a reaction at
/// the far end of the chamber is reconstructed with a beam energy nearly 28 MeV too high and its
/// apparent Ex comes out too large. Since the vertex is uniform, that turns into a spread.
///
/// The bottom row shows what correcting it costs: nothing but the vertex position. Replacing the
/// constant 192 MeV with 192 - 28 z/1000 collapses the width from 0.84 MeV to 0.037 MeV and puts
/// every state on the energy it was generated at.
///
/// Note the correction here uses the TRUE vertex z. In data it would use the reconstructed one,
/// so the gain realised depends on the vertex resolution -- which the digitised simulation, not
/// this macro, has to establish.
///
///   root -b -q 'plot_states_C16pd.C()'

void plot_states_C16pd(TString dir = "/mnt/f/a1975_C16_pd_sim/", TString seed = "s1001",
                       Double_t Ebeam = 192.0, Double_t dEdrift = 28.0, TString tag = "")
{
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());

   const double u = 931.49401;
   const double mb = 16.0147 * u, mt = 1.0078250322 * u, md = 2.0141017781 * u, mr = 15.0105993 * u;

   const int NS = 4;
   const char *TAGS[NS] = {"gs", "ex1", "ex2", "ex3"};
   const double EXT[NS] = {0.0, 0.740, 3.100, 4.660};
   const char *NAME[NS] = {"g.s.", "0.740", "3.10", "4.66"};
   const int COL[NS] = {kBlack, kAzure + 2, kRed + 1, kGreen + 3};

   auto exOf = [&](double Eb, double ke, double th) {
      double E1 = Eb + mb, p1 = std::sqrt(Eb * Eb + 2 * Eb * mb);
      double E3 = ke + md, p3 = std::sqrt(E3 * E3 - md * md);
      double E4 = E1 + mt - E3, p4sq = p1 * p1 + p3 * p3 - 2 * p1 * p3 * std::cos(th);
      return std::sqrt(std::max(E4 * E4 - p4sq, 0.0)) - mr;
   };

   TH2D *hkt[NS];
   TH1D *hRaw[NS], *hCor[NS];
   TH2D *hz = new TH2D("hzgs", "g.s.: reconstructed E_{x} vs vertex z;vertex z [mm];E_{x} [MeV]", 40, 0, 1000, 60,
                       -1.5, 3.0);
   for (int i = 0; i < NS; ++i) {
      hkt[i] = new TH2D(TString::Format("hkt%d", i), "", 90, 0, 45, 120, 0, 60);
      hRaw[i] = new TH1D(TString::Format("hr%d", i), "", 220, -2, 9);
      hCor[i] = new TH1D(TString::Format("hc%d", i), "", 220, -2, 9);
   }

   printf("\n  state   generated   reconstructed at %.0f MeV     corrected with vertex z\n", Ebeam);
   for (int i = 0; i < NS; ++i) {
      TString fn = dir + TAGS[i] + "_" + seed + "_sim.root";
      TFile *f = TFile::Open(fn);
      if (!f || f->IsZombie()) {
         printf("  %-5s MISSING %s\n", TAGS[i], fn.Data());
         continue;
      }
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *mc = nullptr;
      t->SetBranchAddress("MCTrack", &mc);
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020)
               continue; // deuteron
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0)
               continue;
            double ke = std::sqrt(pp * pp + md * md) - md, th = std::acos(pz / pp), z = p->GetStartZ() * 10;
            hkt[i]->Fill(th * TMath::RadToDeg(), ke);
            double exRaw = exOf(Ebeam, ke, th);
            hRaw[i]->Fill(exRaw);
            hCor[i]->Fill(exOf(Ebeam - dEdrift * z / 1000.0, ke, th));
            if (i == 0)
               hz->Fill(z, exRaw);
         }
      }
      printf("  %-5s   %7.3f     %6.2f +- %.3f              %6.2f +- %.3f\n", NAME[i], EXT[i], hRaw[i]->GetMean(),
             hRaw[i]->GetRMS(), hCor[i]->GetMean(), hCor[i]->GetRMS());
      f->Close();
   }

   TCanvas *c = new TCanvas("cpd", "16C(p,d)15C states", 1500, 950);
   c->Divide(2, 2);

   // ---- 1. the kinematic loci
   c->cd(1);
   auto *fr = new TH2D("frk", "deuteron kinematics (MC truth);#theta_{lab} [deg];KE [MeV]", 1, 0, 45, 1, 0, 60);
   fr->Draw();
   auto *lg = new TLegend(0.62, 0.62, 0.89, 0.88);
   for (int i = 0; i < NS; ++i) {
      hkt[i]->SetMarkerStyle(1);
      hkt[i]->SetMarkerColor(COL[i]);
      hkt[i]->Draw("same");
      auto *m = new TMarker(0, 0, 20);
      m->SetMarkerColor(COL[i]);
      lg->AddEntry(m, TString::Format("%s MeV", NAME[i]), "p");
   }
   lg->Draw();

   // ---- 2. Ex as the analysis reconstructs it
   c->cd(2);
   double mx = 0;
   for (int i = 0; i < NS; ++i)
      mx = std::max(mx, hRaw[i]->GetMaximum());
   auto *f2 = new TH1D("fr2", TString::Format("E_{x} at the nominal %.0f MeV (what the analysis does);"
                                              "E_{x} [MeV];deuterons",
                                              Ebeam),
                       1, -2, 9);
   f2->SetMaximum(mx * 1.25);
   f2->Draw();
   for (int i = 0; i < NS; ++i) {
      hRaw[i]->SetLineColor(COL[i]);
      hRaw[i]->SetLineWidth(2);
      hRaw[i]->Draw("hist same");
      auto *l = new TLine(EXT[i], 0, EXT[i], mx * 1.25);
      l->SetLineColor(COL[i]);
      l->SetLineStyle(2);
      l->Draw();
   }
   lg->Draw();

   // ---- 3. the same, with the beam loss corrected by the vertex
   c->cd(3);
   mx = 0;
   for (int i = 0; i < NS; ++i)
      mx = std::max(mx, hCor[i]->GetMaximum());
   auto *f3 = new TH1D("fr3", "E_{x} corrected event-by-event with the vertex z;E_{x} [MeV];deuterons", 1, -2, 9);
   f3->SetMaximum(mx * 1.25);
   f3->Draw();
   for (int i = 0; i < NS; ++i) {
      hCor[i]->SetLineColor(COL[i]);
      hCor[i]->SetLineWidth(2);
      hCor[i]->Draw("hist same");
      auto *l = new TLine(EXT[i], 0, EXT[i], mx * 1.25);
      l->SetLineColor(COL[i]);
      l->SetLineStyle(2);
      l->Draw();
   }
   lg->Draw();

   // ---- 4. where the width comes from
   c->cd(4);
   gPad->SetLogz();
   hz->Draw("colz");
   auto *pf = hz->ProfileX("pfz");
   pf->SetLineColor(kRed + 1);
   pf->SetLineWidth(3);
   pf->SetMarkerColor(kRed + 1);
   pf->SetMarkerStyle(20);
   pf->Draw("same");
   auto *tx = new TLatex(60, 2.55, "dashed lines = generated energies");
   tx->SetTextSize(0.035);
   tx->Draw();

   TString png = here + "/plots/pd_states" + tag + ".png";
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(png);
   printf("\n  wrote %s\n\n", png.Data());
}
