// Step 1 -- export hits from the unpacked ROOT files to CSV.
//
// These files carry their own StreamerInfo, so PLAIN ROOT reads them: no ATTPCROOT
// libraries, no FairRoot, nothing to build.  Check you have ROOT with
//     root-config --version
//
//   root -l -b -q 'dump_hits.C("in.root","out.csv")'          all events
//   root -l -b -q 'dump_hits.C("in.root","out.csv",200)'      first 200, to try it
//
// Columns:  event,pad,tb,x,y,z,q
//   x, y   pad-plane coordinates in mm -- absolutely calibrated, trust these
//   tb     time bucket: the RAW drift time.  z is rebuilt from this downstream, because
//          it depends on the drift velocity, which is not finally settled.
//   z      the z ATTPCROOT itself computed, kept only for comparison
//   q      pulse amplitude, used to weight the line fits
//
// One TTree::Draw per column.  Drawing several at once and reading GetV1..GetV4 looks
// tempting but does not give a reliable 4th column -- it silently mis-pairs the rows.
// Draw twice: the first pass reports how many VALUES there are (which for a hit array is
// far more than the number of entries), the second pass runs with a buffer big enough to
// hold them.  SetEstimate(-1) is NOT enough -- it sizes the buffer to the entry count, so
// with ~160 hits per event the values silently wrap and the columns stop lining up.
std::vector<double> col(TTree *t, const char *expr, Long64_t nmax)
{
   Long64_t n = t->Draw(expr, "", "goff", nmax);
   if (n > t->GetEstimate()) {
      t->SetEstimate(n + 1);
      n = t->Draw(expr, "", "goff", nmax);
   }
   return std::vector<double>(t->GetV1(), t->GetV1() + n);
}

void dump_hits(const char *in, const char *out, Long64_t nEvents = -1)
{
   TFile f(in);
   if (f.IsZombie()) { printf("ERROR cannot open %s\n", in); return; }
   TTree *t = (TTree *)f.Get("cbmsim");
   if (!t) { printf("ERROR no cbmsim tree in %s\n", in); return; }
   Long64_t nmax = (nEvents < 0) ? t->GetEntries() : nEvents;

   const char *P = "AtEventH.fHitArray";
   // Entry$ on its own yields ONE value per entry, not one per hit.  Multiplying by
   // 0 * (a hit-array expression) forces it to be evaluated per hit, so it broadcasts to
   // the hit multiplicity and lines up with the other columns.
   auto ev  = col(t, Form("Entry$+0*%s.fTimeStamp", P), nmax);
   auto x   = col(t, Form("%s.fPosition.fX", P),    nmax);
   auto y   = col(t, Form("%s.fPosition.fY", P),    nmax);
   auto z   = col(t, Form("%s.fPosition.fZ", P),    nmax);
   auto tb  = col(t, Form("%s.fTimeStamp", P),      nmax);
   auto q   = col(t, Form("%s.fCharge", P),         nmax);
   auto pad = col(t, Form("%s.fPadNum", P),         nmax);

   size_t n = x.size();
   if (ev.size() != n || tb.size() != n || q.size() != n) {
      printf("ERROR column lengths disagree (ev %zu, x %zu, tb %zu, q %zu)\n",
             ev.size(), n, tb.size(), q.size());
      return;
   }
   std::ofstream o(out);
   o << "event,pad,tb,x,y,z,q\n";
   for (size_t i = 0; i < n; ++i)
      o << (long)ev[i] << "," << (long)pad[i] << "," << (long)tb[i] << ","
        << Form("%.3f,%.3f,%.3f,%.4g", x[i], y[i], z[i], q[i]) << "\n";
   o.close();
   printf("  %s : %lld events, %zu hits -> %s\n", gSystem->BaseName(in), nmax, n,
          gSystem->BaseName(out));
}

void dump_hits() { printf("usage: root -l -b -q 'dump_hits.C(\"in.root\",\"out.csv\")'\n"); }
