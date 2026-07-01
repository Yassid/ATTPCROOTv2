// Event-aware Ex ntuple builder for the adaptive-TC comparison. Same selection as
// ex_dp_a1975.C (proton hypothesis, PID gate, chi2/ndf + theta + KE cuts, 2-body
// kinematics for 17C) but the ntuple ALSO stores the event number so default and
// loose results can be merged per event.
//
//   root -b -q 'build_ex_ntuple.C("run_0016","","cache_def.root")'        // default fit
//   root -b -q 'build_ex_ntuple.C("run_0016loose","","cache_loose.root")' // loose fit
#include <map>
#include <tuple>

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
// two-body kinematics (verbatim from ex_dp_a1975.C): returns {Ex, theta_cm[deg]}
static std::tuple<double, double> kine_2b(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et2 = m2, Et3 = K_eject + m3, Et4 = Et1 + Et2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * omega2(s, m1 * m1, m2 * m2) * omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (omega2(s, m1 * m1, m2 * m2) * omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

void build_ex_ntuple(TString run = "run_0016", TString fitSuffix = "", TString outCache = "cache_def.root",
                     TString inDir = "/mnt/f/a1975/reco_d2/", TString gateFile = "pid/proton_band_d2_v2.json",
                     double Ebeam = 192.0, double chi2Cut = 10.0, double keMax = 50.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.01410178 * u;
   const double m_p = 1.00782503 * u, m_C17 = 17.0225864 * u;

   bool useGate = (gateFile.Length() > 0);
   AtTools::AtParticleID pid;
   if (useGate) pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TFile *fcache = new TFile(outCache, "RECREATE");
   TNtuple *ntk = new TNtuple("pk", "candidate proton kinematics (d,p)",
                              "event:ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx");

   TString gf = inDir + run + "_genfitter_p" + fitSuffix + ".root";
   if (gSystem->AccessPathName(gf)) { printf("\033[1;31mmissing %s\033[0m\n", gf.Data()); return; }
   TFile *fu = TFile::Open(gf);
   TTree *tu = (TTree *)fu->Get("cbmsim");
   TClonesArray *te = nullptr, *pe = nullptr;
   tu->SetBranchAddress("AtTrackingEvent", &te);
   tu->SetBranchAddress("AtPIDEvent", &pe);

   long nCand = 0;
   Long64_t N = tu->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      tu->GetEntry(i);
      if (pe->GetEntries() == 0 || te->GetEntries() == 0) continue;
      auto *pidev = (AtPIDEvent *)pe->At(0);
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!pidev || !ev) continue;
      std::map<int, AtFittedTrack *> fmap;
      for (auto &ft : ev->GetFittedTracks())
         if (ft) fmap[ft->GetTrackID()] = ft.get();
      for (auto &sr : pidev->GetSpyral()) {
         if (!sr.valid) continue;
         if (useGate && !pid.IsInside(sr.sqrtdEdx, sr.brho)) continue;
         auto it = fmap.find(sr.trackID);
         if (it == fmap.end()) continue;
         auto *ft = it->second;
         auto &k = ft->GetKinematics();
         double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
         double c2n = ndf > 0 ? chi2 / ndf : 1e9;
         double ke = k.kineticEnergy, thRad = k.theta;
         double thDeg = thRad * TMath::RadToDeg();
         if (ke <= 0 || ke > keMax || c2n > chi2Cut) continue;
         if (thDeg < 10.0 || thDeg > 170.0) continue;
         auto v = ft->GetVertex();
         auto [ex, thcm] = kine_2b(m_C16, m_d, m_p, m_C17, Ebeam, thRad, ke);
         if (std::isnan(ex)) continue;
         ++nCand;
         float row[10] = {(float)i, (float)ke, (float)thDeg, (float)v.Z(), (float)thcm, (float)ex,
                          (float)c2n, (float)sr.brho, (float)sr.dEdx, (float)sr.sqrtdEdx};
         ntk->Fill(row);
      }
   }
   fcache->cd(); ntk->Write();
   printf("\033[1;32m%s\033[0m: %lld events, %ld proton candidates -> %s\n", run.Data(), N, nCand, outCache.Data());
   fcache->Close();
}
