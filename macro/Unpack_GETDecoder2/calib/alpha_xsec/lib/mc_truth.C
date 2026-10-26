// MC truth for the tracking-efficiency measurement.
//
// For every simulated reaction (the ODD Geant4 entries, which hold the two outgoing
// alphas) write the true kinematics: the beam energy at the vertex, E_cm, and theta_cm.
// Pairing with the merged events is by index -- merged event i comes from Geant4 entries
// 2i and 2i+1 -- so the tracker's answer can be compared with what was actually there.
//
//   root -l -b -q 'mc_truth.C("sim.root","truth.csv")'
//
// Columns: event,Ecm,theta_cm,E1,E2,ntracks_true
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <fstream>
#include <vector>

void mc_truth(const char *in, const char *out)
{
   TFile *f = TFile::Open(in);
   if (!f || f->IsZombie()) { printf("ERROR cannot open %s\n", in); return; }
   TTree *t = (TTree *)f->Get("cbmsim");
   if (!t) { printf("ERROR no cbmsim tree\n"); return; }
   t->SetEstimate(-1);

   // one Draw per column, with the entry number broadcast to track multiplicity
   auto col = [&](const char *e) {
      Long64_t n = t->Draw(e, "", "goff");
      if (n > t->GetEstimate()) { t->SetEstimate(n + 1); n = t->Draw(e, "", "goff"); }
      return std::vector<double>(t->GetV1(), t->GetV1() + n);
   };
   auto ev = col("Entry$+0*MCTrack.fPx");
   auto px = col("MCTrack.fPx");
   auto py = col("MCTrack.fPy");
   auto pz = col("MCTrack.fPz");
   auto mo = col("MCTrack.fMotherId");

   const double M = 3.727379;  // GeV, alpha
   std::ofstream o(out);
   o << "event,Ecm,theta_cm,E1,E2\n";
   // Track 0 of every entry is the BEAM -- it carries exactly the generated momentum and
   // the full 11.4 MeV.  The two PRODUCTS follow it, and only in the odd entries.  All
   // three have fMotherId = -1, so position in the list is the only discriminator.
   // theta is measured from the BEAM direction, not from z: the beam is tilted 6.47 deg.
   size_t i = 0;
   long written = 0;
   while (i < ev.size()) {
      long e = (long)ev[i];
      size_t j = i;
      while (j < ev.size() && (long)ev[j] == e) ++j;
      if (e % 2 == 1 && j - i >= 3) {                 // beam + two products
         double bx = px[i], by = py[i], bz = pz[i];   // the beam, first in the list
         double bn = sqrt(bx*bx + by*by + bz*bz);
         std::vector<double> E, th;
         for (size_t k = i + 1; k < j; ++k) {         // skip the beam
            double p2 = px[k]*px[k] + py[k]*py[k] + pz[k]*pz[k];
            if (p2 <= 0) continue;
            E.push_back((sqrt(p2 + M*M) - M) * 1000.0);
            double cosang = (px[k]*bx + py[k]*by + pz[k]*bz) / (sqrt(p2) * bn);
            th.push_back(TMath::ACos(TMath::Max(-1.0, TMath::Min(1.0, cosang))) * TMath::RadToDeg());
         }
         if (E.size() >= 2) {
            double Ebeam = E[0] + E[1];               // = the beam energy at the vertex
            double thlab = std::min(th[0], th[1]);    // the forward product
            o << (e - 1) / 2 << "," << Form("%.4f,%.2f,%.4f,%.4f",
                   Ebeam / 2.0, 2.0 * thlab, E[0], E[1]) << "\n";
            ++written;
         }
      }
      i = j;
   }
   o.close();
   printf("  %s : %ld reactions -> %s\n", gSystem->BaseName(in), written, gSystem->BaseName(out));
}

void mc_truth() { printf("usage: mc_truth.C(\"sim.root\",\"truth.csv\")\n"); }
