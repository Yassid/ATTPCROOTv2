/// @file acceptance_fullchain_dt.C
/// @brief 16C(d,t)15C acceptance through the COMPLETE chain the yields pass (2026-11-03): the numerator is a truth-matched fitted track with valid Spyral PID, inside a
///        PID gate, chi2/ndf below the cut and vz_fit inside the window, binned at the true theta_cm. Unlike acceptance_C16dt.C the fit is required.
///        The denominator window is on true z + zOffset (vz_fit = z_true + 19.1 mm, measured), so that both sides describe the same part of the target.
///        Prints the factor of each step per state, so a loss can be attributed. An empty gate skips the gate step. Gate files are arguments: none is adopted here.
///   root -l -b -q 'acceptance_fullchain_dt.C("/mnt/f/a1975/sim_rebuild/reco_g100000/","/mnt/f/a1975/sim_rebuild/fit_g100000/","triton.json","out.root")'
namespace afc {
double om2(double x, double y, double z) { return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z); }
double thetaCm(double m1, double m2, double m3, double K1, double thLab, double K3)
{
   // the true Ex is implied by the MC kinematics of each event, as in acceptance_C16dt.C
   double Et1 = K1 + m1, Et3 = K3 + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1, u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m42 = (std::cos(thLab) * om2(s, m1 * m1, m2 * m2) * om2(u, m2 * m2, m3 * m3) - (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) / (2 * m2 * m2) + s + u - m2 * m2;
   double t = m2 * m2 + m42 - 2 * m2 * Et4;
   return (TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m42) + (m1 * m1 - m2 * m2) * (m3 * m3 - m42)) / (om2(s, m1 * m1, m2 * m2) * om2(s, m3 * m3, m42)))) * 57.29578;
}
}
void acceptance_fullchain_dt(TString recoDir, TString fitDir, TString gateFile, TString outFile, TString statesCSV = "gs,ex1,ex2,ex3,ex4",
                             TString mcDir = "/mnt/f/a1975_C16_dt_sim/", TString mcSuffix = "_s3001_sim.root", double Ebeam = 184.17, double chi2Max = 4.5,
                             double vzLo = 25, double vzHi = 600, double zOffset = 19.1, double bFieldForFit = 0, int nBins = 28, double cmMax = 70, Long64_t maxEv = 0)
{
   gSystem->Load("libAtReconstruction.so"); gSystem->Load("libAtSimulationData.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR"); gErrorIgnoreLevel = kFatal;
   const double mC = 16.0147 * 931.49410, md = 2.01355321 * 931.49410, mt = 3.01550072 * 931.49410;
   const bool useGate = gateFile.Length() > 0;
   AtTools::AtParticleID gate;
   if (useGate) gate = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   TFile fo(outFile, "RECREATE");
   TObjArray *st = statesCSV.Tokenize(",");
   for (int is = 0; is < st->GetEntries(); ++is) {
      TString s = ((TObjString *)st->At(is))->GetString();
      TFile fm(mcDir + s + mcSuffix), fr(recoDir + "t_" + s + "_reco.root"), ff(fitDir + "t_" + s + "_genfitter_t.root");
      auto *tm = (TTree *)fm.Get("cbmsim"); auto *tr = (TTree *)fr.Get("cbmsim"); auto *tf = (TTree *)ff.Get("cbmsim");
      TClonesArray *mc = nullptr, *pa = nullptr, *te = nullptr, *pe = nullptr;
      tm->SetBranchStatus("*", 0); tm->SetBranchStatus("MCTrack*", 1); tm->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchStatus("*", 0); tr->SetBranchStatus("AtPatternEvent*", 1); tr->SetBranchAddress("AtPatternEvent", &pa);
      tf->SetBranchAddress("AtTrackingEvent", &te); tf->SetBranchAddress("AtPIDEvent", &pe);
      TH1D hG("hG_" + s, ";#theta_{cm} [deg]", nBins, 0, cmMax), hN("hN_" + s, ";#theta_{cm} [deg]", nBins, 0, cmMax);
      hG.Sumw2(); hN.Sumw2();
      long nGen = 0, nPat = 0, nVal = 0, nGate = 0, nFit = 0, nChi = 0, nVz = 0;
      Long64_t N = tf->GetEntries();
      // the fit tree repeats events; only entries carrying PID results count, and each event once
      std::set<int> seen;
      const Long64_t nEv = maxEv > 0 ? std::min(maxEv, tm->GetEntries()) : tm->GetEntries();
      for (Long64_t i = 0; i < nEv; ++i) {
         tm->GetEntry(i);
         double ke = -1, th = 0, z = 0;
         for (int j = 0; j < mc->GetEntriesFast(); ++j) {
            auto *m = (AtMCTrack *)mc->At(j);
            if (m && m->GetPdgCode() == 1000010030) { ke = (m->GetEnergy() - m->GetMass()) * 1000.; th = std::acos(m->GetPz() / m->GetP()); z = m->GetStartZ() * 10; break; }
         }
         if (ke < 0 || z + zOffset < vzLo || z + zOffset > vzHi) continue;
         const double cm = afc::thetaCm(mC, md, mt, Ebeam, th, ke);
         ++nGen; hG.Fill(cm);
         tr->GetEntry(i);
         if (!pa->GetEntries() || i >= N) continue;
         std::set<int> ids;
         for (auto &t : ((AtPatternEvent *)pa->At(0))->GetTrackCand()) {
            int n = 0, nt = 0;
            for (auto &h : t.GetHitArray()) { auto &ms = h->GetMCSimPointArray(); if (ms.empty()) continue; ++n; if (ms[0].A == 3 && ms[0].Z == 1) ++nt; }
            if (n && nt > n / 2) ids.insert(t.GetTrackID());
         }
         if (ids.empty()) continue;
         ++nPat;
         tf->GetEntry(i);
         if (!pe->GetEntries() || !te->GetEntries()) continue;
         auto *pid = (AtPIDEvent *)pe->At(0);
         if (!pid) continue;
         bool v = false, g = false, f = false, c = false, w = false;
         auto *ev = (AtTrackingEvent *)te->At(0);
         for (auto &sr : pid->GetSpyral()) {
            if (!sr.valid || !ids.count(sr.trackID)) continue;
            v = true;
            if (useGate && !gate.IsInside(sr.sqrtdEdx, sr.brho)) continue;
            g = true;
            for (auto &ft : ev->GetFittedTracks()) {
               if (!ft || ft->GetTrackID() != sr.trackID) continue;
               const double ndf = ft->GetTrackMetadata()->GetNdf();
               if (ft->GetKinematicsXtr().kineticEnergy <= 0 || ndf <= 0) continue;
               f = true;
               if (ft->GetTrackMetadata()->GetChi2() / ndf >= chi2Max) continue;
               c = true;
               const double vz = ft->GetVertex().Z();
               if (vz >= vzLo && vz <= vzHi) w = true;
            }
         }
         nVal += v; nGate += g; nFit += f; nChi += c; nVz += w;
         if (w) hN.Fill(cm);
      }
      auto rat = [](long a, long b) { return b ? double(a) / b : 0.; };
      printf("%-4s gen %6ld | pattern %.3f | Spyral valid %.3f | gate %.3f | fit %.3f | chi2 %.3f | vz %.3f || OVERALL %.3f\n", s.Data(), nGen, rat(nPat, nGen), rat(nVal, nPat), rat(nGate, nVal), rat(nFit, nGate), rat(nChi, nFit), rat(nVz, nChi), rat(nVz, nGen));
      fo.cd();
      TH1D hA(hN); hA.SetName("hAcc_" + s); hA.Divide(&hN, &hG, 1, 1, "B");
      hG.Write(); hN.Write(); hA.Write();
      for (int b = 1; b <= nBins; ++b) if (hG.GetBinContent(b) >= 20) printf("   %5.1f %7.0f %7.0f  %.3f +- %.3f\n", hG.GetBinCenter(b), hG.GetBinContent(b), hN.GetBinContent(b), hA.GetBinContent(b), hA.GetBinError(b));
   }
   fo.Close();
}
