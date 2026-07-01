// Test AtPSAMultiFit on a synthetic pad with TWO OVERLAPPING GET pulses.
// Shows it recovers BOTH peaks (deconvolving the overlap) where AtPSAMax finds only one.
// Run: root -l 'test_multifit.C'   (ATTPCROOT classes autoload after gSystem->Load; no #includes)
#include <cmath>

// GET nominal response (same form as AtNominalResponse / AtPSAMultiFit)
static double Resp(double t, double t0, double tauTB)
{
   double u = (t - t0) / tauTB;
   if (u <= 0)
      return 0;
   return std::exp(-3.0 * u) * std::sin(u) * u * u * u;
}

void test_multifit()
{
   gSystem->Load("libAtReconstruction");

   // ---- 1. FairRun + AtDigiPar (so AtPSA::Init() works) ----
   TString dir = gSystem->Getenv("VMCWORKDIR");
   TString parFile = dir + "/parameters/ATTPC.a1975_deuterium.par";
   auto *run = new FairRunAna();
   auto *rtdb = run->GetRuntimeDb();
   auto *parIo = new FairParAsciiFileIo();
   parIo->open(parFile.Data(), "in");
   rtdb->setFirstInput(parIo);
   auto *par = dynamic_cast<AtDigiPar *>(rtdb->getContainer("AtDigiPar"));
   rtdb->initContainers(0);
   double TBTime = par->GetTBTime();                  // ns/TB
   double tauUs = 0.720;                              // peaking time
   double tauTB = tauUs * 1000.0 / TBTime;            // tau in TB
   // response peak (reduced time + value) for placing/labeling
   double du = 0.001, best = 0, bestU = 0;
   for (double u = du; u < M_PI; u += du) {
      double r = std::exp(-3 * u) * std::sin(u) * u * u * u;
      if (r > best) { best = r; bestU = u; }
   }
   double peakOff = bestU * tauTB;
   printf("\nTBTime=%.1f ns  tau=%.2f us = %.1f TB  peakOffset=%.1f TB\n", TBTime, tauUs, tauTB, peakOff);

   // ---- 2. synthetic pad: TWO overlapping pulses ----
   // pulse starts t0; observed PEAK is at t0+peakOff. Choose peaks ~1.2*tau apart (overlapping).
   double t0a = 150.0, Aa = 520.0;                    // pulse A: Aa = desired PEAK ADC
   double t0b = t0a + 2.8 * tauTB, Ab = 340.0;        // pulse B (overlaps A)
   double peakA = t0a + peakOff, peakB = t0b + peakOff;
   auto *rng = new TRandom3(42);
   AtPad pad(420);
   pad.SetPadCoord(ROOT::Math::XYPoint(25.0, -40.0)); // arbitrary valid (x,y)
   pad.SetSizeID(1);
   pad.SetPedestalSubtracted(kTRUE);
   for (int i = 0; i < 512; ++i) {
      // coeff = peakADC / responsePeakValue so the observed peak equals Aa,Ab
      double v = (Aa / best) * Resp(i, t0a, tauTB) + (Ab / best) * Resp(i, t0b, tauTB) + rng->Gaus(0, 4.0);
      pad.SetADC(i, v);
   }
   printf("TRUTH: pulse A peak@%.1f TB (A=%.0f)   pulse B peak@%.1f TB (A=%.0f)   |  peak separation %.1f TB\n",
          peakA, Aa, peakB, Ab, peakB - peakA);

   // ---- 3. run AtPSAMultiFit ----
   AtPSAMultiFit mf;
   mf.SetPeakingTime(tauUs);
   mf.SetThreshold(40);
   mf.SetMaxPeaks(4);
   mf.SetMinSeparation(4);
   mf.Init();
   auto mfHits = mf.AnalyzePad(&pad);

   printf("\n=== AtPSAMultiFit: %zu hits ===\n", mfHits.size());
   for (auto &h : mfHits)
      printf("   peakTB=%.1f   z=%.1f mm   amp=%.0f   Q=%.0f\n", (double)h->GetTimeStamp(), h->GetPosition().Z(),
             h->GetCharge(), h->GetTraceIntegral());

   // ---- 4. run AtPSAMax for comparison ----
   AtPSAMax mx;
   mx.SetThreshold(40);
   mx.Init();
   auto mxHits = mx.AnalyzePad(&pad);
   printf("=== AtPSAMax: %zu hit(s) ===\n", mxHits.size());
   for (auto &h : mxHits)
      printf("   peakTB=%.1f   z=%.1f mm   amp=%.0f\n", (double)h->GetTimeStamp(), h->GetPosition().Z(), h->GetCharge());

   // ---- 5. plot ----
   auto *gTrace = new TGraph();
   for (int i = 100; i < 260; ++i) gTrace->SetPoint(gTrace->GetN(), i, pad.GetADC(i));
   gTrace->SetTitle("AtPSAMultiFit on overlapping double-pulse;time bucket;ADC");
   gTrace->SetMarkerStyle(20); gTrace->SetMarkerSize(0.5);

   auto *c = new TCanvas("c", "c", 1000, 600);
   gTrace->Draw("AP");
   // fitted components (reconstruct from hits: A=amp/peakVal, t0=peakTB-peakOff)
   int col[] = {kRed, kGreen + 2, kMagenta, kOrange};
   int ci = 0;
   for (auto &h : mfHits) {
      double A = h->GetCharge() / best;
      double t0 = h->GetTimeStamp() - peakOff;
      auto *f = new TF1(Form("comp%d", ci), [A, t0, tauTB](double *x, double *) { return A * Resp(x[0], t0, tauTB); },
                        100, 260, 0);
      f->SetLineColor(col[ci % 4]); f->SetLineWidth(2); f->Draw("same"); ci++;
   }
   // total fit (sum)
   auto *fSum = new TF1("sum", [&](double *x, double *) {
      double s = 0;
      for (auto &h : mfHits) s += (h->GetCharge() / best) * Resp(x[0], h->GetTimeStamp() - peakOff, tauTB);
      return s; }, 100, 260, 0);
   fSum->SetLineColor(kBlue); fSum->SetLineStyle(2); fSum->Draw("same");
   // AtPSAMax single peak (vertical line)
   for (auto &h : mxHits) {
      auto *l = new TLine(h->GetTimeStamp(), 0, h->GetTimeStamp(), 350);
      l->SetLineColor(kGray + 2); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw();
   }
   auto *leg = new TLegend(0.62, 0.6, 0.89, 0.88);
   leg->AddEntry(gTrace, "trace (2 pulses+noise)", "p");
   leg->AddEntry((TObject *)nullptr, Form("MultiFit: %zu hits", mfHits.size()), "");
   leg->AddEntry(fSum, "fit sum", "l");
   leg->AddEntry((TObject *)nullptr, Form("Max: %zu hit (dashed)", mxHits.size()), "");
   leg->Draw();
   c->SaveAs("test_multifit.png");
   printf("\nwrote test_multifit.png\n");
}
