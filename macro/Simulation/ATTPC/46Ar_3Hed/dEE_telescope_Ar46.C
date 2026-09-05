/// @file dEE_telescope_Ar46.C
/// @brief dE-E from the forward telescope for 46Ar(3He,d)47K.
///
///   root -b -q 'dEE_telescope_Ar46.C("gs.root,ex2020.root","g.s.,2.02 MeV")'
///
/// Sums the deposit per LAYER per TRACK. That is not bookkeeping fussiness: AtSiArray writes TWO
/// AtSiPoints per traversal -- an entry point carrying zero energy and an exit point carrying the
/// whole deposit -- so anything that averages over points rather than summing per track comes out
/// a factor two low. Reading them as one hit each is how a "90 % acceptance" got quoted here once
/// when the truth was 46 %.
///
/// Layers are identified by VOLUME NAME, which is the contract set in geometry/Ar46_telescope.C:
/// silicon_Ar46_dE, silicon_Ar46_E, CsI_Ar46_i_j.
///
/// Everything is keyed on the MC-truth track, so the species is known rather than inferred; this
/// is a detector-response plot, not a demonstration of particle identification. A real PID plot
/// would have to work without that.
/// @param dEum  dE thickness in um, FOR THE AXIS LABEL ONLY. It is not read from the geometry, so
///              passing the wrong value mislabels the plot without changing anything it shows --
///              which is exactly what happened when a 20 um run was drawn as "500 um".
void dEE_telescope_Ar46(TString files, TString labels, TString outPng = "plots/dEE_telescope_Ar46.png",
                        Double_t dEum = 500.)
{
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);

   std::unique_ptr<TObjArray> af(files.Tokenize(","));
   std::unique_ptr<TObjArray> al(labels.Tokenize(","));
   const int NF = af->GetEntries();

   // y range follows the dE thickness: a 20 um layer deposits ~40 MeV and would be an
   // invisible line at the bottom of a 700 MeV axis. Frame on the data.
   const double yHi = (dEum > 200.) ? 700. : 120.;
   auto *hAll = new TH2D("hAll", "", 200, 0, 700, 200, 0, yHi);
   std::vector<TGraph *> gr;
   int cols[4] = {kBlack, kRed + 1, kAzure + 2, kGreen + 3};

   printf("\n%-14s %8s %10s %10s %10s %10s\n", "file", "47K", "in dE", "in E", "in CsI", "mean dE");
   for (int i = 0; i < NF; ++i) {
      TString fn = ((TObjString *)af->At(i))->GetString();
      TFile *f = TFile::Open(fn);
      if (!f || f->IsZombie()) { printf("  cannot open %s\n", fn.Data()); continue; }
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *si = nullptr, *mc = nullptr;
      t->SetBranchAddress("AtSiArrayPoint", &si);
      t->SetBranchAddress("MCTrack", &mc);

      auto *g = new TGraph();
      long nK = 0, ndE = 0, nE = 0, nCs = 0;
      double sde = 0;
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         int kid = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *m = (AtMCTrack *)mc->At(k);
            if (m->GetMotherId() == -1 && m->GetPdgCode() == 1000190470) { kid = k; break; }
         }
         if (kid < 0) continue;
         ++nK;
         // sum per layer for THIS track (see the header: two points per traversal)
         double de = 0, ee = 0, cs = 0;
         for (int k = 0; k < si->GetEntriesFast(); ++k) {
            auto *p = (AtSiPoint *)si->At(k);
            if (p->GetTrackID() != kid) continue;
            double el = p->GetEnergyLoss() * 1000.0; // GeV -> MeV
            if (el <= 0) continue;
            TString v = p->GetVolName();
            if (v.Contains("_dE")) de += el;
            else if (v.Contains("CsI")) cs += el;
            else ee += el;
         }
         if (de > 0) { ++ndE; sde += de; }
         if (ee > 0) ++nE;
         if (cs > 0) ++nCs;
         // the abscissa is the energy left AFTER the dE, which is what a real dE-E plots
         if (de > 0) { g->SetPoint(g->GetN(), ee + cs, de); hAll->Fill(ee + cs, de); }
      }
      printf("%-14s %8ld %10ld %10ld %10ld %10.1f\n", gSystem->BaseName(fn), nK, ndE, nE, nCs,
             ndE ? sde / ndE : 0.0);
      g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetMarkerColor(cols[i % 4]);
      gr.push_back(g);
      f->Close();
   }

   TCanvas *c = new TCanvas("cdee", "dE-E", 1250, 560);
   c->Divide(2, 1, 0.001, 0.001);

   c->cd(1);
   gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
   TH1F *fr = gPad->DrawFrame(0, 0, 700, yHi);
   fr->SetTitle(Form("^{47}K in the telescope;E after the #DeltaE (E + CsI)  [MeV];#DeltaE  (%.0f #mum Si)  [MeV]", dEum));
   fr->GetYaxis()->SetTitleOffset(1.3);
   auto *lg = new TLegend(0.45, 0.72, 0.95, 0.90);
   lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.036);
   for (int i = 0; i < (int)gr.size(); ++i) {
      gr[i]->Draw("P same");
      lg->AddEntry(gr[i], (i < al->GetEntries()) ? ((TObjString *)al->At(i))->GetString() : "", "p");
   }
   lg->Draw();
   // the line every dE-E plot is read against: everything stopping in the dE sits on E = 0
   auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.033); tx->SetTextColor(kGray + 3);
   tx->DrawLatex(0.16, 0.20, "points on E = 0 stopped in the #DeltaE:");
   tx->DrawLatex(0.16, 0.155, "no #DeltaE-E information for those");

   c->cd(2);
   gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.13);
   gPad->SetLogz();
   hAll->SetTitle(Form("all states together;E after the #DeltaE (E + CsI)  [MeV];#DeltaE  (%.0f #mum Si)  [MeV]", dEum));
   hAll->GetYaxis()->SetTitleOffset(1.3);
   hAll->Draw("colz");

   gSystem->mkdir("plots", kTRUE);
   c->SaveAs(outPng);
   printf("\n  wrote %s\n\n", outPng.Data());
}
