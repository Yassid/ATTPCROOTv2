/// @file ab_catima_diag.C
/// @brief Why is chi2/ndf degenerate in the CATIMA-off arm? Look at ndf and chi2 raw.
///
/// A -42% width change is exactly the size that a heavy fit-failure tail has faked before on
/// this data, so before quoting it the two arms have to be shown to be fitting comparably --
/// not one of them quietly producing unconverged tracks that the other does not.

#include <algorithm>
#include <cmath>
#include <vector>

namespace {
void dump(TString file, const char *tag)
{
   TFile *f = TFile::Open(file);
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);
   long n = 0, ndfZero = 0, ndfNeg = 0, chi2Nan = 0, chi2Zero = 0;
   std::vector<double> ndfs, chi2s, ratios;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!te || te->GetEntries() == 0)
         continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev)
         continue;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft)
            continue;
         const auto &m = ft->GetTrackMetadata();
         if (!m)
            continue;
         const double ndf = m->GetNdf(), chi2 = m->GetChi2();
         ++n;
         if (ndf == 0)
            ++ndfZero;
         if (ndf < 0)
            ++ndfNeg;
         if (std::isnan(chi2))
            ++chi2Nan;
         if (chi2 == 0)
            ++chi2Zero;
         ndfs.push_back(ndf);
         chi2s.push_back(chi2);
         if (ndf > 0 && !std::isnan(chi2))
            ratios.push_back(chi2 / ndf);
      }
   }
   auto med = [](std::vector<double> v) {
      if (v.empty())
         return 0.0;
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };
   printf("%-12s n=%-5ld  ndf==0: %-5ld  ndf<0: %-4ld  chi2 NaN: %-4ld  chi2==0: %-4ld\n", tag, n, ndfZero, ndfNeg,
          chi2Nan, chi2Zero);
   printf("%-12s median ndf=%8.2f   median chi2=%10.3f   median chi2/ndf=%8.4f  (n valid=%zu)\n", "", med(ndfs),
          med(chi2s), med(ratios), ratios.size());
   f->Close();
}
} // namespace

void ab_catima_diag(TString fOff, TString fOn)
{
   gSystem->Load("libAtReconstruction.so");
   dump(fOff, "CATIMA off");
   dump(fOn, "CATIMA on");
}
