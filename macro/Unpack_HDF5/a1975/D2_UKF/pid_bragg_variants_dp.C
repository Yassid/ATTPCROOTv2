/// @file pid_bragg_variants_dp.C
/// @brief dE/dx segment kept OFF the Bragg rise, on the SAME 3D estimate per track: the adopted
///        first-arc segment (s0) against segments that stop L mm before the far end of the track, or at a
///        fraction f of its length from the vertex (AtHelix3DCore firstArcPID stopExcl / keepFrac).
///        Brho is common, only the dE/dx segment changes. Reference branch names as pid_variants_dp.C, so
///        pid_proton_fixed.C / pid_dd_generic.C / pid_gate_check.C score it unchanged.
///   root -l -b -q 'pid_bragg_variants_dp.C("<reco>","out.root")'
#include "../UKF/pid/helix3d.h"
void pid_bragg_variants_dp(TString in, TString out, Long64_t maxEv = 1000000000)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;
   struct V { const char *name; double L, f; };
   std::vector<V> var = {{"adopted", 0, 1}, {"stop-20mm", 20, 1}, {"stop-40mm", 40, 1},
                         {"stop-60mm", 60, 1}, {"keep 75%", 0, 0.75}, {"keep 60%", 0, 0.60}};
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
               for (int k = 1; k < K; ++k)
                  if (h3d::firstArcPID(r, V.end, 152.0, bb, ss, th, pp, var[k].L, var[k].f)) { s[k] = ss; b[k] = bb; ns[k] = 1; }
            }
         }
         o.Fill();
      }
   }
   fo.cd(); o.Write(); fo.Close();
   printf("%s: %lld tracks\n", out.Data(), o.GetEntries());
}
