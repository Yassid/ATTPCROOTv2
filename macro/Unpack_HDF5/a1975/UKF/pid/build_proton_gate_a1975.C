/// @file build_proton_gate_a1975.C
/// @brief Build a proton-band gate that ENCLOSES the band, data-driven, on the
///        DEFAULT (Spyral) sqrt(dEdx)-vs-brho plane.
///
/// For each sqrt(dEdx) slice it finds the proton-band brho peak and the brho range
/// where the density stays above a fraction of that peak (the band envelope). The
/// polygon = upper envelope (left->right) + lower envelope (right->left), padded a
/// little so the band is fully inside. Overlays it, reports the enclosed fraction,
/// and writes proton_band.json (AtParticleID, sqrtdEdx vs brho).
///
///   root -b -q 'build_proton_gate_a1975.C()'

#include <vector>

void build_proton_gate_a1975(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,"
                                               "run_0114,run_0115",
                            TString inDir = "/mnt/f/a1975/reco/", double frac = 0.25, double pad = 0.012,
                            Bool_t write = true)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double xlo = 3.0, xhi = 40.0;
   // upper-left cap: raise the upper envelope on the left so the steep high-brho
   // continuation of the proton band (low sqrt(dEdx)) is enclosed.
   const double xCapEnd = 9.0, capLeft = 1.05, capRight = 0.42;
   TH2F *h = new TH2F("h", "proton gate (data-driven);#sqrt{dEdx};B#rho [T m]", 360, 0, 40, 400, 0, 1.4);
   std::vector<std::array<double, 4>> rows; // event-light: store nothing; we use h

   // ---- fill the plane ----
   std::vector<float> sx, sb; // keep tracks to compute enclosed fraction later
   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString fn = inDir + run + "_pid.root";
      if (gSystem->AccessPathName(fn))
         continue;
      TFile *f = TFile::Open(fn);
      TTree *t = (TTree *)f->Get("cbmsim");
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPIDEvent", &pe);
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &r : ev->GetDefault())
            if (r.valid) {
               h->Fill(r.sqrtdEdx, r.brho);
               sx.push_back(r.sqrtdEdx);
               sb.push_back(r.brho);
            }
      }
      f->Close();
   }

   // ---- per sqrt(dEdx) slice: find proton-band peak + envelope ----
   std::vector<double> ex, eHi, eLo;             // slice center, upper brho, lower brho
   const int nSlice = 24;
   double dx = (xhi - xlo) / nSlice;
   for (int s = 0; s < nSlice; ++s) {
      double xc = xlo + (s + 0.5) * dx;
      int bx1 = h->GetXaxis()->FindBin(xlo + s * dx);
      int bx2 = h->GetXaxis()->FindBin(xlo + (s + 1) * dx) - 1;
      TH1D *py = h->ProjectionY("py", bx1, bx2);
      // proton band = brightest peak below brho 0.75 (the lower band)
      int bpk = 0;
      double pk = 0;
      int bSearchHi = py->GetXaxis()->FindBin(0.75);
      for (int b = py->GetXaxis()->FindBin(0.07); b <= bSearchHi; ++b)
         if (py->GetBinContent(b) > pk) {
            pk = py->GetBinContent(b);
            bpk = b;
         }
      if (pk < 4) {
         delete py;
         continue;
      } // too few -> skip slice
      double thr = frac * pk;
      int bUp = bpk, bDn = bpk;
      while (bUp < py->GetNbinsX() && py->GetBinContent(bUp + 1) > thr)
         ++bUp;
      while (bDn > 1 && py->GetBinContent(bDn - 1) > thr)
         --bDn;
      double bhi = py->GetXaxis()->GetBinUpEdge(bUp) + pad;
      double blo = py->GetXaxis()->GetBinLowEdge(bDn) - pad;
      if (blo < 0.05)
         blo = 0.05;
      // raise the upper edge on the left to enclose the high-brho continuation
      if (xc < xCapEnd) {
         double t = (xc - xlo) / (xCapEnd - xlo);
         double cap = capLeft + t * (capRight - capLeft);
         if (cap > bhi)
            bhi = cap;
      }
      ex.push_back(xc);
      eHi.push_back(bhi);
      eLo.push_back(blo);
      delete py;
   }

   // ---- assemble polygon: upper edge L->R, lower edge R->L ----
   std::vector<std::pair<double, double>> poly;
   for (size_t i = 0; i < ex.size(); ++i)
      poly.push_back({ex[i], eHi[i]});
   for (int i = (int)ex.size() - 1; i >= 0; --i)
      poly.push_back({ex[i], eLo[i]});
   AtTools::AtCut2D cut("proton_band", poly, "sqrtdEdx", "brho");

   long nIn = 0;
   for (size_t i = 0; i < sx.size(); ++i)
      if (cut.IsInside(sx[i], sb[i]))
         ++nIn;
   printf("\ntracks=%zu  enclosed by gate=%ld (%.1f%%)  (%zu polygon vertices)\n", sx.size(), nIn,
          sx.size() ? 100.0 * nIn / sx.size() : 0, poly.size());

   TGraph *g = new TGraph;
   for (size_t k = 0; k < poly.size(); ++k)
      g->SetPoint(k, poly[k].first, poly[k].second);
   g->SetPoint(poly.size(), poly[0].first, poly[0].second);
   g->SetLineColor(kRed);
   g->SetLineWidth(3);

   TCanvas *c = new TCanvas("c", "gate", 1000, 750);
   c->SetLogz();
   h->GetYaxis()->SetRangeUser(0, 1.0);
   h->Draw("colz");
   g->Draw("L same");
   c->SaveAs("pid/plots/proton_gate.png");
   printf("saved proton_gate.png\n");

   if (write) {
      AtTools::AtParticleID(cut, 1, 1).WriteJSON("proton_band.json");
      printf("wrote proton_band.json\n");
   }
}
