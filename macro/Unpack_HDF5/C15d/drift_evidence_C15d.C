/// @file drift_evidence_C15d.C
/// @brief The four plots behind the drift working point, and why the far anchor is not measurable.
///
///   root -b -q 'drift_evidence_C15d.C()'
///
/// Panel 1  hit arrival time, all hits. The PAD-PLANE end is a sharp, stable lower edge.
/// Panel 2  the LATEST-arriving charge per track. This is where the window should show up, and
///          instead there is a SPIKE against the end of the 512-sample record: charge that should
///          arrive later is lost, so the apparent edge is the DAQ boundary, not the chamber.
/// Panel 3  beam energy from the 15C(d,d) elastic ridge as a function of the assumed drift scale.
///          The measured anchor (TBentrance 498) gives 163 MeV; the run log's ~190 MeV needs
///          TBentrance ~583, i.e. a drift of 563 time buckets in a 512-bucket record.
/// Panel 4  the ridge itself with the elastic locus at both working points.
///
/// The conclusion the plots carry: a truncated record means the far anchor is a LOWER BOUND, so
/// mm/TB is an upper bound and the beam energy a lower bound -- which is the direction and roughly
/// the size of the discrepancy. The drift velocity has to come from outside this data.

namespace {
const double kU = 931.49410242, kM1 = 15.0105993 * kU, kMD = 2.0141018 * kU;
double Trec(double thDeg, double Eb)
{
   const double E1 = Eb + kM1, p2 = E1 * E1 - kM1 * kM1, c = std::cos(thDeg * TMath::DegToRad());
   const double n = 2 * kMD * p2 * c * c, d = (E1 + kMD) * (E1 + kMD) - p2 * c * c;
   return d > 0 ? n / d : -1;
}
} // namespace

void drift_evidence_C15d(TString recoDir = "/home/yassid/C15d_reco/",
                         TString kinGlob = "/tmp/claude-1000/-home-yassid/d13cc5a6-0089-40bf-9a9f-1ca75ae973e7/scratchpad/qf2/*_kin_d.root",
                         TString outDir = "plots/", Double_t mmPerTB = 2.0912, Double_t zPP = 1000.0,
                         Double_t tbEnt = 498.0, Int_t nRecoRuns = 3)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->mkdir(outDir, kTRUE);
   gStyle->SetOptStat(0);

   // ---- panels 1 and 2: arrival-time distributions -----------------------------------------
   auto *hAll = new TH1D("hAll", "all hits;arrival time bucket;hits", 540, 0, 540);
   auto *hMax = new TH1D("hMax", "latest-arriving charge per track;time bucket;tracks", 540, 0, 540);
   TSystemDirectory dir(recoDir, recoDir);
   TList *fl = dir.GetListOfFiles();
   std::vector<TString> use;
   if (fl) {
      TIter next(fl);
      while (auto *o = dynamic_cast<TSystemFile *>(next())) {
         TString n = o->GetName();
         if (!o->IsDirectory() && n.EndsWith("_reco.root")) use.push_back(recoDir + n);
      }
   }
   std::sort(use.begin(), use.end());
   int done = 0;
   for (auto &f : use) {
      if (done >= nRecoRuns) break;
      TFile fi(f);
      auto *t = (TTree *)fi.Get("cbmsim");
      if (!t) continue;
      TClonesArray *pa = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pa);
      const Long64_t n = std::min((Long64_t)4000, t->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         t->GetEntry(i);
         if (!pa || pa->GetEntriesFast() == 0) continue;
         auto *p = (AtPatternEvent *)pa->At(0);
         if (!p) continue;
         for (auto &trk : p->GetTrackCand()) {
            double mx = -1e9;
            for (auto &h : trk.GetHitArray()) {
               const double tb = tbEnt - (zPP - h->GetPosition().Z()) / mmPerTB;
               hAll->Fill(tb);
               if (tb > mx) mx = tb;
            }
            if (mx > -1e8) hMax->Fill(mx);
         }
      }
      ++done;
   }

   // ---- panels 3 and 4: the ridge -----------------------------------------------------------
   TChain ch("kin");
   ch.Add(kinGlob);
   Double_t th, ke, c2;
   ch.SetBranchAddress("thetaXtr", &th);
   ch.SetBranchAddress("keXtr", &ke);
   ch.SetBranchAddress("chi2ndf", &c2);
   std::vector<double> PT, PZ;
   auto *hKT = new TH2D("hKT", "15C(d,d') deuterons;#theta_{lab} [deg];KE [MeV]", 150, 20, 90, 150, 0, 60);
   for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
      ch.GetEntry(i);
      if (!(ke > 0) || c2 > 5) continue;
      hKT->Fill(th, ke);
      if (ke < 3 || ke > 90 || th < 30 || th > 75) continue;
      const double p = std::sqrt((ke + kMD) * (ke + kMD) - kMD * kMD), t = th * TMath::DegToRad();
      PT.push_back(p * std::sin(t));
      PZ.push_back(p * std::cos(t));
   }
   auto *gEb = new TGraph();
   for (double k = 1.25; k >= 0.78; k -= 0.025) {
      double bE = 0;
      long best = -1;
      for (double E = 80; E <= 340; E += 1) {
         long n = 0;
         for (size_t i = 0; i < PT.size(); ++i) {
            const double pz = k * PZ[i], pt = PT[i], p = std::sqrt(pt * pt + pz * pz);
            const double m = Trec(std::atan2(pt, pz) * TMath::RadToDeg(), E);
            if (m > 0 && std::fabs(std::sqrt(p * p + kMD * kMD) - kMD - m) / m < 0.08) ++n;
         }
         if (n > best) { best = n; bE = E; }
      }
      gEb->SetPoint(gEb->GetN(), 20.0 + 1000.0 / (mmPerTB * k), bE); // TBentrance vs Ebeam
   }

   auto *c = new TCanvas("cdrift", "drift evidence", 1500, 1050);
   c->Divide(2, 2);

   c->cd(1);
   gPad->SetLogy();
   hAll->SetLineColor(kAzure + 2);
   hAll->Draw("hist");
   auto *l1 = new TLine(20, 0, 20, hAll->GetMaximum());
   l1->SetLineColor(kGreen + 2); l1->SetLineWidth(2); l1->Draw();
   auto *tx1 = new TLatex(); tx1->SetNDC(); tx1->SetTextSize(0.036);
   tx1->DrawLatex(0.18, 0.86, "pad plane: sharp, stable at TB ~20");

   c->cd(2);
   gPad->SetLogy();
   hMax->SetLineColor(kRed + 1);
   hMax->Draw("hist");
   auto *l2 = new TLine(512, 0, 512, hMax->GetMaximum());
   l2->SetLineColor(kBlack); l2->SetLineStyle(2); l2->SetLineWidth(2); l2->Draw();
   auto *tx2 = new TLatex(); tx2->SetNDC(); tx2->SetTextSize(0.034);
   tx2->DrawLatex(0.14, 0.86, "SPIKE against the 512-sample record end");
   tx2->DrawLatex(0.14, 0.81, "=> truncation, not the chamber window");

   c->cd(3);
   gEb->SetTitle("beam energy from the elastic ridge;assumed TB of the entrance window;E_{beam} [MeV]");
   gEb->SetMarkerStyle(20); gEb->SetLineWidth(2);
   gEb->Draw("ALP");
   auto *l3 = new TLine(gEb->GetXaxis()->GetXmin(), 190, gEb->GetXaxis()->GetXmax(), 190);
   l3->SetLineColor(kGreen + 2); l3->SetLineWidth(2); l3->Draw();
   auto *l4 = new TLine(512, gEb->GetYaxis()->GetXmin(), 512, gEb->GetYaxis()->GetXmax());
   l4->SetLineColor(kBlack); l4->SetLineStyle(2); l4->SetLineWidth(2); l4->Draw();
   auto *tx3 = new TLatex(); tx3->SetNDC(); tx3->SetTextSize(0.034);
   tx3->DrawLatex(0.16, 0.85, "run log ~190 MeV");
   tx3->DrawLatex(0.16, 0.80, "needs TB_{ent} ~583 > 512 record");

   c->cd(4);
   gPad->SetLogz();
   hKT->Draw("colz");
   int col[2] = {kRed + 1, kGreen + 2};
   double ebs[2] = {163, 190};
   for (int i = 0; i < 2; ++i) {
      auto *g = new TGraph();
      for (double a = 25; a <= 88; a += 1) {
         const double v = Trec(a, ebs[i]);
         if (v > 0 && v < 60) g->SetPoint(g->GetN(), a, v);
      }
      g->SetLineColor(col[i]); g->SetLineWidth(3);
      g->Draw("L same");
   }
   auto *tx4 = new TLatex(); tx4->SetNDC(); tx4->SetTextSize(0.034);
   tx4->SetTextColor(kRed + 1);  tx4->DrawLatex(0.55, 0.85, "163 MeV (measured anchor)");
   tx4->SetTextColor(kGreen + 2); tx4->DrawLatex(0.55, 0.80, "190 MeV (run log)");

   c->SaveAs(outDir + "drift_evidence_C15d.png");
   std::cout << "wrote " << outDir << "drift_evidence_C15d.png\n";
}
