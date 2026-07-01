/// @file ex_dp_a1975.C
/// @brief 16C excitation-energy spectrum from PROTON-hypothesis genfit fits — the
///        16C(p,p)17C neutron-stripping channel (a1975 D2-target runs 0016-0103).
///
/// Inverse kinematics: 16C beam + d target -> p + 17C; the detected light ejectile
/// is the PROTON (going largely BACKWARD in the lab, see project memory). The fit
/// must have been run with the proton hypothesis (fitGenfitter_a1975_deuterium.C
/// -> reco_d2/<run>_genfitter_p.root). Two-body kinematics with m1=16C, m2=d,
/// m3=p, m4=17C reconstructs the 16C excitation energy. A genuine (p,p) stripping
/// populates the 17C bound states (g.s. 3/2+, 0.217 1/2+, 0.335 5/2+); junk fit as
/// protons does NOT close kinematics there.
///
///   root -b -q 'ex_dp_a1975.C("run_0016","/mnt/f/a1975/reco_d2/")'
///   root -b -q 'ex_dp_a1975.C("run_0016,run_0017,run_0018,run_0019","/mnt/f/a1975/reco_d2/")'
///
/// PID selection: the in-fitter Spyral PID (AtPIDEvent, the SAME plane the gate was
/// built on) gated by pid/proton_band_d2.json, plus proton hypothesis + chi2/ndf cut
/// + kinematics. There is NO ion-chamber gate (the D2 reco has no _FRIB.root); the
/// proton gate + good chi2 + kinematic closure are the selection.
///
/// theta convention: the clean AtGenfitter stores theta in RADIANS (unlike the old
/// AtGenfit which stored degrees) -- kine_2b is fed k.theta directly, NO DegToRad.
///
/// Ebeam / exShift: Ebeam is the 16C energy at the reaction vertex (~192 MeV in this
/// chamber, same beam as the (p,d) work). exShift slides the g.s. to Ex=0; leave 0
/// for the first pass, then read the g.s. peak and set the calibration.

#include <map>
#include <tuple>
#include <vector>

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

// two-body kinematics (verbatim from C16_pp_ana.C): returns {Ex, theta_cm[deg]}
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

void ex_elastic_pp(TString runsCSV = "run_0016", TString inDir = "/mnt/f/a1975/reco_d2/",
                 TString fitSuffix = "", TString gateFile = "pid/proton_band_d2_v2.json", double Ebeam = 192.0,
                 double chi2Cut = 10.0, double exShift = 0.0, double keMax = 50.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.01410178 * u;
   const double m_p = 1.00782503 * u, m_C17 = 17.0225864 * u;

   bool useGate = (gateFile.Length() > 0);
   AtTools::AtParticleID pid;
   if (useGate)
      pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TFile *fcache = new TFile("proton_kin_pp.root", "RECREATE");
   TNtuple *ntk =
      new TNtuple("pk", "candidate proton kinematics (p,p)", "ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx");

   TH1F *hex = new TH1F(
      "hex", "^{16}C excitation energy (p,p), proton hyp;E_{x}(^{16}C) [MeV];proton candidates", 200, -10, 20);
   TH2F *hexcm = new TH2F("hexcm", "E_{x}(^{16}C) vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -10,
                          20);
   TH2F *hbt = new TH2F("hbt", "B#rho vs #theta_{lab} (proton cand.);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200,
                        0, 1.2);
   TH2F *hpid = new TH2F("hpid", "PID plane (proton cand.);#sqrt{dEdx};B#rho [T m]", 200, 0, 50, 200, 0, 1.2);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nCand = 0, nValPID = 0, nGated = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString gf = inDir + run + "_genfitter_p" + fitSuffix + ".root";
      if (gSystem->AccessPathName(gf)) {
         printf("skip %s (missing %s)\n", run.Data(), gf.Data());
         continue;
      }
      TFile *fu = TFile::Open(gf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      TClonesArray *te = nullptr, *pe = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      tu->SetBranchAddress("AtPIDEvent", &pe);

      Long64_t N = tu->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         tu->GetEntry(i);
         if (pe->GetEntries() == 0 || te->GetEntries() == 0)
            continue;
         auto *pidev = (AtPIDEvent *)pe->At(0);
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!pidev || !ev)
            continue;
         std::map<int, AtFittedTrack *> fmap; // fitted tracks by trackID
         for (auto &ft : ev->GetFittedTracks())
            if (ft)
               fmap[ft->GetTrackID()] = ft.get();
         for (auto &sr : pidev->GetSpyral()) {
            if (!sr.valid)
               continue;
            ++nValPID;
            if (useGate && !pid.IsInside(sr.sqrtdEdx, sr.brho))
               continue;
            ++nGated;
            auto it = fmap.find(sr.trackID);
            if (it == fmap.end())
               continue;
            auto *ft = it->second;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double c2n = ndf > 0 ? chi2 / ndf : 1e9;
            double ke = k.kineticEnergy, thRad = k.theta;
            double thDegCut = thRad * TMath::RadToDeg();
            if (ke <= 0 || ke > keMax || c2n > chi2Cut)
               continue;
            if (thDegCut < 10.0 || thDegCut > 170.0) // drop near-beam forward/backward
               continue;
            auto v = ft->GetVertex();
            auto [ex, thcm] = kine_2b(m_C16, m_p, m_p, m_C16, Ebeam, thRad, ke);
            if (std::isnan(ex))
               continue;
            ++nCand;
            double thDeg = thRad * TMath::RadToDeg();
            double exCal = ex + exShift;
            hex->Fill(exCal);
            hexcm->Fill(thcm, exCal);
            hbt->Fill(thDeg, sr.brho);
            hpid->Fill(sr.sqrtdEdx, sr.brho);
            ntk->Fill(ke, thDeg, v.Z(), thcm, ex, c2n, sr.brho, sr.dEdx, sr.sqrtdEdx);
         }
      }
      fu->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nvalid Spyral PID: %ld -> in proton gate: %ld -> proton candidates (p,p) -> 17C Ex: %ld\n", nValPID, nGated,
          nCand);
   fcache->cd();
   ntk->Write();
   printf("cached -> proton_kin_pp.root (ntuple pk: ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx)\n");

   TCanvas *c = new TCanvas("c", "ex_dp", 1550, 1040);
   c->Divide(2, 2);
   c->cd(1);
   hex->Draw();
   TLine *l0 = new TLine(0, 0, 0, hex->GetMaximum());
   l0->SetLineColor(kRed);
   l0->SetLineStyle(2);
   l0->Draw();
   c->cd(2);
   hexcm->Draw("colz");
   c->cd(3);
   hbt->Draw("colz");
   c->cd(4);
   hpid->Draw("colz");
   c->SaveAs("plots/ex_elastic_pp.png");
   printf("saved plots/ex_elastic_pp.png  (red dashed line in panel 1 = 16C g.s. Ex=0)\n");
}
