/// @file ex_pd_a1975.C
/// @brief 15C excitation-energy spectrum from DEUTERON-hypothesis UKF fits — the
///        16C(p,d)15C neutron-pickup channel.
///
/// Mirror of ex_a1975.C but for the (p,d) channel: the detected light ejectile is
/// a DEUTERON, so the UKF must have been run with the deuteron hypothesis
/// (fitUKF_a1975.C(...,"deuteron",...,"_d",...) -> <run>_ukf_d.root). Two-body
/// kinematics with m1=16C, m2=p, m3=d, m4=15C reconstructs the 15C excitation
/// energy. A GENUINE (p,d) pickup populates a sharp 15C ground-state peak at
/// Ex~0; protons mis-fit as deuterons do NOT close kinematics there. This is the
/// decisive test that the (p,d) channel is real (vs the dEdx "deuteron band"
/// which we showed is an arclength artifact).
///
///   root -b -q 'ex_pd_a1975.C("run_0106","/mnt/f/a1975/reco/")'
///
/// Default reads <run>_ukf_d.root (deuteron fits). No dEdx PID gate by default —
/// the deuteron hypothesis + good chi2 + IC beam gate + kinematics ARE the
/// selection. The surviving-track brho-vs-dedx plane is also drawn so a deuteron
/// gate can be built later if a real locus appears.

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

void ex_pd_a1975(TString runsCSV = "run_0106", TString inDir = "/mnt/f/a1975/reco/", TString fitSuffix = "_d",
                 TString gateFile = "pid/deuteron_band.json", double Ebeam = 192.0, double icMin = 950,
                 double icMax = 1350, int icTbLo = 1000, int icTbHi = 1350, double chi2Cut = 5.0, double bField = 2.85,
                 TString fitDir = "/mnt/f/a1975/reco_pd/", double exShift = -0.615)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double u = 931.49401;
   const double m_C16 = 16.0147 * u, m_p = 1.007825 * u;
   const double m_d = 2.01410178 * u, m_C15 = 15.0105993 * u;

   bool useGate = (gateFile.Length() > 0);
   AtTools::AtParticleID pid;
   if (useGate)
      pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   TFile *fcache = new TFile("deuteron_kin.root", "RECREATE");
   TNtuple *ntk =
      new TNtuple("dk", "candidate deuteron kinematics", "ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx");

   TH1F *hex =
      new TH1F("hex", "^{15}C excitation energy (p,d), calibrated;E_{x}(^{15}C) [MeV];deuteron candidates", 200, -10, 20);
   TH2F *hexcm = new TH2F("hexcm", "E_{x}(^{15}C) vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -10,
                          20);
   TH2F *hbt = new TH2F("hbt", "B#rho vs #theta_{lab} (deuteron cand.);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 200,
                        0, 2.0);
   TH2F *hpid = new TH2F("hpid", "PID plane (deuteron cand.);#sqrt{dEdx};B#rho [T m]", 200, 0, 40, 200, 0, 2.0);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nCand = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      // deuteron fits live in fitDir (reco_pd/), FRIB stays in inDir (reco/)
      TString fdir = (fitDir.Length() ? fitDir : inDir);
      TString uf = fdir + run + "_ukf" + fitSuffix + ".root", ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(uf) || gSystem->AccessPathName(ff)) {
         printf("skip %s (missing %s or %s)\n", run.Data(), uf.Data(), ff.Data());
         continue;
      }
      TFile *fu = TFile::Open(uf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      TClonesArray *te = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      TFile *fc = TFile::Open(ff);
      TTree *tc = (TTree *)fc->Get("cbmsim");
      TClonesArray *re = nullptr;
      tc->SetBranchAddress("AtRawEvent", &re);

      Long64_t N = std::min(tu->GetEntries(), tc->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         tc->GetEntry(i);
         double ic = -1;
         if (re->GetEntries() > 0) {
            auto *raw = (AtRawEvent *)re->At(0);
            if (raw && !raw->GetGenTraces().empty()) {
               auto &adc = raw->GetGenTraces()[0]->GetADC();
               double mx = -1e9;
               for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b)
                  mx = std::max(mx, adc[b]);
               ic = mx;
            }
         }
         if (ic < icMin || ic > icMax)
            continue; // IC gate (16C beam)
         tu->GetEntry(i);
         if (te->GetEntries() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         std::vector<AtTrack> orig = ev->GetTrackArray();
         std::map<int, AtTrack *> byID;
         for (auto &tr : orig)
            byID[tr.GetTrackID()] = &tr;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double c2n = ndf > 0 ? chi2 / ndf : 1e9;
            double ke = k.kineticEnergy, thRad = k.theta;
            if (ke <= 0 || ke > 1000 || c2n > chi2Cut)
               continue;
            auto it = byID.find(ft->GetTrackID());
            if (it == byID.end())
               continue;
            auto r = spy.Estimate(*it->second);
            if (!r.valid)
               continue;
            if (useGate && !pid.IsInside(r.sqrtdEdx, r.brho))
               continue;
            auto [ex, thcm] = kine_2b(m_C16, m_p, m_d, m_C15, Ebeam, thRad, ke);
            if (std::isnan(ex))
               continue;
            ++nCand;
            double thDeg = thRad * TMath::RadToDeg();
            double exCal = ex + exShift; // g.s.->0 calibration (see pd/cal_pd_a1975.C)
            hex->Fill(exCal);
            hexcm->Fill(thcm, exCal);
            hbt->Fill(thDeg, r.brho);
            hpid->Fill(r.sqrtdEdx, r.brho);
            // cache stores RAW ex (downstream view_pd/cal_pd apply their own shift)
            ntk->Fill(ke, thDeg, r.vertex.Z(), thcm, ex, c2n, r.brho, r.dEdx, r.sqrtdEdx);
         }
      }
      fu->Close();
      fc->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\ndeuteron candidates -> 15C Ex: %ld\n", nCand);
   fcache->cd();
   ntk->Write();
   printf("cached -> deuteron_kin.root (ntuple dk: ke:theta:vertexz:thcm:ex:chi2ndf:brho:dedx:sqrtdedx)\n");

   TCanvas *c = new TCanvas("c", "ex_pd", 1550, 1040);
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
   c->SaveAs("pd/plots/ex_pd_spectrum.png");
   printf("saved ex_pd_spectrum.png  (red dashed line in panel 1 = 15C g.s. Ex=0)\n");
}
