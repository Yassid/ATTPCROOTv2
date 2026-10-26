// Full MC truth for the tracking-efficiency and resolution measurements.
//
// Track 0 of every Geant4 entry is the BEAM; the two PRODUCTS follow it, and only in the
// odd entries.  All three carry fMotherId = -1, so position in the list is the only
// discriminator.  Angles are measured from the BEAM direction because the beam is tilted
// 6.47 deg -- measuring from z would bias every theta_cm by twice that.
//
// Merged event i comes from Geant4 entries 2i (beam) and 2i+1 (products).
//
//   root -l -b -q 'mc_truth2.C("sim.root","truth2.csv")'
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TSystem.h>
#include <fstream>
#include <vector>
#include <algorithm>

void mc_truth2(const char *in, const char *out)
{
   TFile *f = TFile::Open(in);
   if (!f || f->IsZombie()) { printf("ERROR cannot open %s\n", in); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("ERROR no cbmsim tree\n"); return; }
   t->SetEstimate(-1);
   auto col = [&](const char *e) {
      Long64_t n = t->Draw(e, "", "goff");
      if (n > t->GetEstimate()) { t->SetEstimate(n + 1); n = t->Draw(e, "", "goff"); }
      return std::vector<double>(t->GetV1(), t->GetV1() + n);
   };
   auto ev = col("Entry$+0*MCTrack.fPx");
   auto px = col("MCTrack.fPx"), py = col("MCTrack.fPy"), pz = col("MCTrack.fPz");
   auto vx = col("MCTrack.fStartX"), vy = col("MCTrack.fStartY"), vz = col("MCTrack.fStartZ");

   const double M = 3.727379;  // GeV, alpha
   std::ofstream o(out);
   o << "event,Ecm,theta_cm,E1,E2,theta1,theta2,"
        "bdx,bdy,bdz,d1x,d1y,d1z,d2x,d2y,d2z,vx,vy,vz\n";
   size_t i = 0; long written = 0;
   while (i < ev.size()) {
      long e = (long)ev[i];
      size_t j = i;
      while (j < ev.size() && (long)ev[j] == e) ++j;
      if (e % 2 == 1 && j - i >= 3) {
         double bn = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
         double bd[3] = {px[i]/bn, py[i]/bn, pz[i]/bn};
         std::vector<double> E, th; std::vector<std::vector<double>> d;
         for (size_t k = i + 1; k < j; ++k) {
            double p = sqrt(px[k]*px[k] + py[k]*py[k] + pz[k]*pz[k]);
            if (p <= 0) continue;
            E.push_back((sqrt(p*p + M*M) - M) * 1000.0);
            d.push_back({px[k]/p, py[k]/p, pz[k]/p});
            double c = d.back()[0]*bd[0] + d.back()[1]*bd[1] + d.back()[2]*bd[2];
            th.push_back(TMath::ACos(std::max(-1.0, std::min(1.0, c))) * TMath::RadToDeg());
         }
         if (E.size() >= 2) {
            int a = 0, b = 1;
            if (th[1] < th[0]) { a = 1; b = 0; }   // a = the forward product
            o << (e - 1) / 2 << ","
              << Form("%.4f,%.3f,%.4f,%.4f,%.3f,%.3f", (E[0]+E[1])/2.0, 2.0*th[a],
                      E[a], E[b], th[a], th[b]) << ","
              << Form("%.5f,%.5f,%.5f", bd[0], bd[1], bd[2]) << ","
              << Form("%.5f,%.5f,%.5f", d[a][0], d[a][1], d[a][2]) << ","
              << Form("%.5f,%.5f,%.5f", d[b][0], d[b][1], d[b][2]) << ","
              << Form("%.3f,%.3f,%.3f", vx[i+1], vy[i+1], vz[i+1]) << "\n";
            ++written;
         }
      }
      i = j;
   }
   o.close();
   printf("  %s : %ld reactions -> %s\n", gSystem->BaseName(in), written, gSystem->BaseName(out));
}
