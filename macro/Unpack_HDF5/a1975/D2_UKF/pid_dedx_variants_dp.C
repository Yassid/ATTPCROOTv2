/// @file pid_dedx_variants_dp.C
/// @brief dE/dx definitions on the SAME 3D estimate per track: the adopted total charge / helix path
///        (s0) against helix-path TRUNCATED MEANS (sample length L, kept fraction f). Brho is common (b3),
///        so only the dE/dx definition changes. Reference branch names as pid_variants_dp.C, so
///        pid_proton_fixed.C / pid_dd_generic.C / pid_gate_check.C score it unchanged.
///   root -l -b -q 'pid_dedx_variants_dp.C("<reco>","out.root")'
#include "../UKF/pid/helix3d.h"
void pid_dedx_variants_dp(TString in, TString out, Long64_t maxEv = 1000000000)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   struct V { const char *name; double L, f; };
   std::vector<V> var = {{"total (adopted)", 0, 0}, {"trunc L8 f0.7", 8, 0.7}, {"trunc L8 f0.6", 8, 0.6},
                         {"trunc L8 f0.8", 8, 0.8}, {"trunc L12 f0.7", 12, 0.7}, {"trunc L5 f0.7", 5, 0.7}};
   const int K = var.size();
   TFile f(in);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pa = nullptr;
   t->SetBranchStatus("*", 0);
   t->SetBranchStatus("AtPatternEvent*", 1);
   t->SetBranchAddress("AtPatternEvent", &pa);
   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);
   TFile fo(out, "RECREATE");
   TTree o("w", "3D PID, dE/dx definitions");
   int nh, va, v3, ns[16];
   double br, sqa, pol, b3, sq3h, th3, sp3, s[16], b[16];
   o.Branch("nh", &nh); o.Branch("va", &va); o.Branch("br", &br); o.Branch("sqa", &sqa); o.Branch("pol", &pol);
   o.Branch("v3", &v3); o.Branch("b3", &b3); o.Branch("sq3h", &sq3h); o.Branch("th3", &th3); o.Branch("sp3", &sp3);
   for (int k = 0; k < K; ++k) {
      o.Branch(Form("s%d", k), &s[k]); o.Branch(Form("b%d", k), &b[k]); o.Branch(Form("ns%d", k), &ns[k]);
   }
   for (Long64_t ev = 0; ev < std::min(maxEv, t->GetEntries()); ++ev) {
      if (t->GetEntry(ev) <= 0 || !pa || !pa->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
         auto &H = tr.GetHitArray();
         nh = H.size();
         if (nh < 30) continue;
         std::vector<h3d::P> h;
         for (auto &x : H) { auto p = x->GetPosition(); h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()}); }
         auto sr = spy.Estimate(tr);
         va = sr.valid; br = sr.valid ? sr.brho : -1; sqa = sr.valid ? sr.sqrtdEdx : -1;
         pol = sr.valid ? sr.polar * TMath::RadToDeg() : -999;
         v3 = 0; b3 = sq3h = th3 = -1; sp3 = 1e9;
         for (int k = 0; k < K; ++k) { s[k] = b[k] = -1; ns[k] = 0; }
         auto r = h3d::estimate(h, 2.85, -1, 1000, 5, 4.0, -4.0, 8, 0.15, 0.5, h3d::EndOpt());
         if (r.ok) {
            int spyEnd = -1;
            if (sr.valid) spyEnd = (sr.direction == 0) == (h[r.order.front()].z < h[r.order.back()].z) ? 0 : 1;
            auto V = h3d::majorityEnd(r, spyEnd);
            double bb, ss, th, pp;
            if (h3d::firstArcPID(r, V.end, 152.0, bb, ss, th, pp)) {
               v3 = 1; b3 = bb; sq3h = ss; th3 = th * TMath::RadToDeg(); sp3 = r.length / r.nBackbone;
               s[0] = ss; b[0] = bb;
               for (int k = 1; k < K; ++k) {
                  int n;
                  if (h3d::firstArcPIDTrunc(r, V.end, 152.0, var[k].L, var[k].f, bb, ss, th, pp, n)) { s[k] = ss; b[k] = bb; ns[k] = n; }
               }
            }
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
