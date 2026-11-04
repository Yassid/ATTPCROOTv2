/// @file geom_vars_dp.C
/// Per event of the (d,p) g.s. sim: where the true proton goes in the gas, from the AtTpcPoint steps. Tree "g", lengths in mm, z in the Monte Carlo frame. Input of geom_corr_dp.C.
///   root -l -b -q 'geom_vars_dp.C("sim.root","out.root")'
namespace gv { const int NZONE = 10, NRING = 7; const double RE[NRING + 1] = {0, 30, 60, 100, 150, 200, 250, 400}; }
void geom_vars_dp(TString simFile = "/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root", TString out = "/mnt/f/a1975/sim_rebuild/dp_geom/truth_vars.root") {
   using namespace gv;
   gSystem->Load("libAtSimulationData.so"); gErrorIgnoreLevel = kError;
   TFile f(simFile); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *mc = nullptr, *pt = nullptr;
   t->SetBranchStatus("*", 0); t->SetBranchStatus("MCTrack*", 1); t->SetBranchStatus("AtTpcPoint*", 1); t->SetBranchAddress("MCTrack", &mc); t->SetBranchAddress("AtTpcPoint", &pt);
   TFile fo(out, "RECREATE"); TTree g("g", "true proton path per event");
   int ev; double tvz, tth, tke, R0, pitch, turns, path, zend, rend, rmax, zrmax, zapex, rapex, zret, rret, nturnz, dzturn, zone[NZONE], ring[NRING];
   g.Branch("ev", &ev); g.Branch("tvz", &tvz); g.Branch("tth", &tth); g.Branch("tke", &tke); g.Branch("R0", &R0); g.Branch("pitch", &pitch); g.Branch("turns", &turns); g.Branch("path", &path); g.Branch("zend", &zend); g.Branch("rend", &rend); g.Branch("rmax", &rmax); g.Branch("zrmax", &zrmax); g.Branch("zapex", &zapex); g.Branch("rapex", &rapex); g.Branch("zret", &zret); g.Branch("rret", &rret); g.Branch("nturnz", &nturnz); g.Branch("dzturn", &dzturn); g.Branch("zone", zone, Form("zone[%d]/D", NZONE)); g.Branch("ring", ring, Form("ring[%d]/D", NRING));
   for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); int ip = -1; AtMCTrack *m = nullptr;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) { auto *q = (AtMCTrack *)mc->At(j); if (q && q->GetPdgCode() == 2212) { ip = j; m = q; break; } }
      if (ip < 0) continue;
      const double p = m->GetP() * 1000; ev = i; tvz = m->GetStartZ() * 10; tth = acos(m->GetPz() / m->GetP()) * TMath::RadToDeg(); tke = sqrt(p * p + 938.272 * 938.272) - 938.272;
      const double st = sin(tth * TMath::DegToRad()), ct = cos(tth * TMath::DegToRad()); R0 = p * st / (0.299792 * 2.85); pitch = 2 * TMath::Pi() * R0 * fabs(ct / std::max(st, 1e-3)); nturnz = pitch > 0 ? tvz / pitch : 0;
      std::fill(zone, zone + NZONE, 0.); std::fill(ring, ring + NRING, 0.); path = 0; turns = 0; rmax = 0; zrmax = tvz; zapex = rapex = zret = rret = -1;
      double px = 0, py = 0, pz = 0, ph = 0, pr = 0, rmin0 = 0; bool first = true, hasph = false, up = false, gotapex = false, gotret = false; int n = 0;
      for (int k = 0; k < pt->GetEntriesFast(); ++k) { auto *a = (AtMCPoint *)pt->At(k); if (a->GetTrackID() != ip) continue;
         double x = a->GetX() * 10, y = a->GetY() * 10, z = a->GetZ() * 10, r = hypot(x, y);
         if (!first) { double dx = x - px, dy = y - py, dz = z - pz, ds = sqrt(dx * dx + dy * dy + dz * dz); path += ds; int iz = int(z / 100.); if (iz >= 0 && iz < NZONE) zone[iz] += ds; for (int b = 0; b < NRING; ++b) if (r >= RE[b] && r < RE[b + 1]) ring[b] += ds;
            if (hypot(dx, dy) > 0.2) { double h = atan2(dy, dx); if (hasph) { double d = h - ph; while (d > TMath::Pi()) d -= 2 * TMath::Pi(); while (d < -TMath::Pi()) d += 2 * TMath::Pi(); turns += fabs(d) / (2 * TMath::Pi()); } ph = h; hasph = true; }
            if (r > pr + 1e-3) up = true; if (!gotapex && up && r < pr - 0.05 && pr > 5) { gotapex = true; zapex = pz; rapex = pr; rmin0 = r; }
            if (gotapex && !gotret) { if (r < rmin0) rmin0 = r; else if (r > rmin0 + 0.05) { gotret = true; zret = pz; rret = pr; } } }
         if (r > rmax) { rmax = r; zrmax = z; } px = x; py = y; pz = z; pr = r; first = false; ++n; }
      zend = pz; rend = pr; dzturn = turns > 0 ? (tvz - zend) / turns : 0; if (n > 0) g.Fill(); }
   fo.cd(); g.Write(); printf("wrote %lld events to %s\n", g.GetEntries(), out.Data()); fo.Close(); }
