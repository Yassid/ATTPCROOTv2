/// @file spyral_dir_audit.C
/// @brief Is AtSpyralPID's forward/backward call wrong on a1975 DATA? Scored with the opposite call.
///
/// AtSpyralPID::Direction() decides the direction from the winding sense of the azimuth along the
/// z-sorted point cloud (fraction of steps with increasing angle > fDirThreshold). Everything after
/// it -- the vertex end, the first arc, the polar sign, failCodes 8 and 9 -- hangs on that call, and
/// nothing checks it. failCodes 8+9 reject 32 % of H2 data tracks and ~0 % of the simulation.
///
/// THE OPPOSITE CALL WITHOUT TOUCHING THE CLASS: mirroring x -> -x reverses the winding sense and
/// leaves every other quantity Estimate() uses (z, rho, radius, charge, arc length) unchanged. So
/// Estimate(mirrored track) IS Estimate with the direction forced the other way.
///
/// THE REFEREE IS PHYSICS, NOT A PROXY: on H2, 16C(p,p) recoil protons have a known theta_lab-KE
/// locus. Whichever call puts a track on Ex(16C) ~ 0 is the right one. Data hit frame: theta_lab =
/// 180 - polar, vz = 1000 - z (census_arc_dp.C convention).
///
/// Production settings: class defaults + |B| = 2.85, exactly as pidPass_a1975.C builds it.
///
///   root -l -b -q 'pid/spyral_dir_audit.C(106,108,"pid/plots/spyral_dir_audit.root",4000)'

namespace {
const double kU = 931.49410242, kMe = 0.51099895;
const double kMp = 938.272088;
const double kMC16 = 16.014701 * kU - 6 * kMe;

// missing-mass Ex of 16C for a recoil proton (theta_lab rad, KE MeV) from 16C on p at Eb
double exPP(double Eb, double th, double ke)
{
   double Eb_t = Eb + kMC16, pb = std::sqrt(Eb_t * Eb_t - kMC16 * kMC16);
   double Ep = ke + kMp, pp = std::sqrt(Ep * Ep - kMp * kMp);
   double E4 = Eb_t + kMp - Ep;
   double px = -pp * std::sin(th), pz = pb - pp * std::cos(th);
   double m2 = E4 * E4 - px * px - pz * pz;
   return m2 > 0 ? std::sqrt(m2) - kMC16 : -999;
}

// winding fraction as Direction() computes it, without its spline pre-smoothing (diagnostic only)
double windFrac(std::vector<std::array<double, 3>> p)
{
   std::sort(p.begin(), p.end(), [](auto &a, auto &b) { return a[2] < b[2]; });
   int n = p.size();
   if (n < 3) return -1;
   double mx = 0, my = 0;
   for (auto &q : p) { mx += q[0]; my += q[1]; }
   mx /= n; my /= n;
   double Suu = 0, Svv = 0, Suv = 0, Suuu = 0, Svvv = 0, Suvv = 0, Suuv = 0;
   for (auto &q : p) {
      double u = q[0] - mx, v = q[1] - my;
      Suu += u * u; Svv += v * v; Suv += u * v;
      Suuu += u * u * u; Svvv += v * v * v; Suvv += u * v * v; Suuv += u * u * v;
   }
   double det = Suu * Svv - Suv * Suv;
   if (std::abs(det) < 1e-12) return -1;
   double bu = 0.5 * (Suuu + Suvv), bv = 0.5 * (Suuv + Svvv);
   double cx = (bu * Svv - bv * Suv) / det + mx, cy = (Suu * bv - Suv * bu) / det + my;
   int pos = 0;
   double prev = std::atan2(p[0][1] - cy, p[0][0] - cx);
   if (prev < 0) prev += 2 * M_PI;
   for (int i = 1; i < n; ++i) {
      double a = std::atan2(p[i][1] - cy, p[i][0] - cx);
      if (a < 0) a += 2 * M_PI;
      if (a - prev > 0) ++pos;
      prev = a;
   }
   return double(pos) / (n - 1);
}
} // namespace

void spyral_dir_audit(int runLo = 106, int runHi = 108, TString out = "pid/plots/spyral_dir_audit.root",
                      Long64_t maxEv = 4000, double Eb = 185.0, TString inDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   gErrorIgnoreLevel = kFatal;

   AtTools::AtSpyralPID spy;
   spy.SetBField(2.85);

   TFile fo(out, "RECREATE");
   TTree o("d", "Spyral direction audit: normal vs forced-opposite call");
   int run, tid, nh, geoBack;
   Long64_t ev;
   double geoTh, wf, qLo, qHi, rLo, rHi;
   // n = normal call, f = forced opposite (mirrored)
   int nv, nf, nd, fv, ff, fd;
   double nth, nbr, nde, nvz, nke, nex, fth, fbr, fde, fvz, fke, fex;
   o.Branch("run", &run); o.Branch("ev", &ev); o.Branch("tid", &tid); o.Branch("nh", &nh);
   o.Branch("geoTh", &geoTh); o.Branch("geoBack", &geoBack); o.Branch("wf", &wf);
   o.Branch("qLo", &qLo); o.Branch("qHi", &qHi); o.Branch("rLo", &rLo); o.Branch("rHi", &rHi);
   o.Branch("nv", &nv); o.Branch("nf", &nf); o.Branch("nd", &nd); o.Branch("nth", &nth); o.Branch("nbr", &nbr);
   o.Branch("nde", &nde); o.Branch("nvz", &nvz); o.Branch("nke", &nke); o.Branch("nex", &nex);
   o.Branch("fv", &fv); o.Branch("ff", &ff); o.Branch("fd", &fd); o.Branch("fth", &fth); o.Branch("fbr", &fbr);
   o.Branch("fde", &fde); o.Branch("fvz", &fvz); o.Branch("fke", &fke); o.Branch("fex", &fex);

   auto fill = [&](const AtTools::AtSpyralResult &r, int &v, int &fc, int &d, double &th, double &br, double &de, double &vz,
                   double &ke, double &ex) {
      v = r.valid; fc = r.failCode; d = r.direction;
      th = br = de = vz = ke = ex = -999;
      if (!r.valid) return;
      th = 180 - r.polar * TMath::RadToDeg();
      vz = 1000 - r.vertex.Z();
      br = r.brho; de = r.sqrtdEdx;
      double p = 299.792458 * r.brho;
      ke = std::sqrt(p * p + kMp * kMp) - kMp;
      ex = exPP(Eb, th * TMath::DegToRad(), ke);
   };

   for (run = runLo; run <= runHi; ++run) {
      TFile *f = TFile::Open(Form("%s/run_%04d_reco.root", inDir.Data(), run));
      if (!f || f->IsZombie()) continue;
      auto *t = (TTree *)f->Get("cbmsim");
      TClonesArray *pa = nullptr;
      t->SetBranchStatus("*", 0);
      t->SetBranchStatus("AtPatternEvent*", 1);
      t->SetBranchAddress("AtPatternEvent", &pa);
      Long64_t N = std::min(maxEv, t->GetEntries());
      for (ev = 0; ev < N; ++ev) {
         if (t->GetEntry(ev) <= 0 || !pa || pa->GetEntries() == 0) continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe) continue;
         for (auto &tr : pe->GetTrackCand()) {
            auto &hits = tr.GetHitArray();
            nh = hits.size();
            if (nh < 5) continue;
            tid = tr.GetTrackID();
            geoTh = tr.GetGeoTheta() * TMath::RadToDeg();
            geoBack = geoTh > 90;
            std::vector<std::array<double, 4>> pz;
            for (auto &h : hits) { auto p = h->GetPosition(); pz.push_back({p.X(), p.Y(), p.Z(), h->GetCharge()}); }
            std::sort(pz.begin(), pz.end(), [](auto &a, auto &b) { return a[2] < b[2]; });
            // charge and |xy| at the two z ends (15 % of points each): Bragg / beam-axis markers
            int k = std::max(1, nh * 15 / 100);
            qLo = qHi = rLo = rHi = 0;
            for (int i = 0; i < k; ++i) {
               qLo += pz[i][3]; rLo += std::hypot(pz[i][0], pz[i][1]);
               qHi += pz[nh - 1 - i][3]; rHi += std::hypot(pz[nh - 1 - i][0], pz[nh - 1 - i][1]);
            }
            qLo /= k; qHi /= k; rLo /= k; rHi /= k;
            std::vector<std::array<double, 3>> p3;
            for (auto &q : pz) p3.push_back({q[0], q[1], q[2]});
            wf = windFrac(p3);

            fill(spy.Estimate(tr), nv, nf, nd, nth, nbr, nde, nvz, nke, nex);
            AtTrack m;
            m.SetTrackID(tid);
            for (auto &h : hits) {
               AtHit c(*h);
               auto p = c.GetPosition();
               c.SetPosition({-p.X(), p.Y(), p.Z()});
               m.AddHit(c);
            }
            fill(spy.Estimate(m), fv, ff, fd, fth, fbr, fde, fvz, fke, fex);
            o.Fill();
         }
      }
      f->Close();
      printf("run %d done, %lld tracks so far\n", run, o.GetEntries());
   }
   fo.cd();
   o.Write();
   fo.Close();
}
