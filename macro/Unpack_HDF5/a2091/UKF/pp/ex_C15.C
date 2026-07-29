/// @file ex_C15.C
/// @brief 15C excitation-energy spectrum from UKF-fitted protons: 15C(p,p')15C*.
///
/// Ported/simplified from a1975 UKF/pp/ex_a1975.C. Uses the UKF proton (KE, theta) +
/// two-body kinematics to reconstruct the 15C excitation energy. Elastic -> Ex~0;
/// excited states appear as peaks.
///
/// FIRST-PASS SIMPLIFICATIONS (vs a1975 ex_a1975.C) -- all flagged for later:
///   * NO FRIB IC beam gate  -- FRIB aux files not yet produced for a2091 (needs the
///     _FRIBDAQ unpacker). Add once <run>_FRIB.root exist.
///   * NO proton PID gate    -- no 15C proton_band.json yet. Build from the sqrt(dEdx)
///     vs Brho plane (AtSpyralPID), then gate as in a1975.
///   * fit-quality cut only  -- chi2/ndf < chi2Cut and physical-KE window.
///
/// >>> Ebeam DEFAULT (195 MeV) is a PLACEHOLDER taken from the old a2091 C15_pp_anaFit.C
///     (Ebeam_buff = 195; ~13 MeV/u for 15C). CONFIRM the real a2091 15C beam energy and
///     pass it explicitly; the ABSOLUTE Ex scale depends on it. <<<
///
///   root -b -q 'ex_C15.C("run_0138", "/home/yassid/a2091_C15_reco/", 195.0)'
///   root -b -q 'ex_C15.C("run_0138,run_0058,run_0060", "/home/yassid/a2091_C15_reco/", 195.0)'

#include <map>
#include <tuple>
#include <vector>

static double omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

// two-body kinematics (verbatim from C16_pp_ana.C / ex_a1975.C): returns {Ex, theta_cm[deg]}
static std::tuple<double, double> kine_2b(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
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

/// Ejectile/residual masses default to (p,p'); for 15C(p,d)14C pass
///   15C(p,d)14C -> mEjectAmu = 2.013553 (deuteron), mResidAmu = 14.003242 (14C).
void ex_C15(TString runsCSV = "run_0138", TString inDir = "/home/yassid/a2091_C15_reco/", double Ebeam = 195.0,
             double chi2Cut = 5.0, TString outTag = "", double mEjectAmu = 1.007825, double mResidAmu = 15.0105993,
             TString chanLabel = "", TString fitter = "ukf")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   std::cout << "\033[1;31m>>> Ebeam = " << Ebeam
             << " MeV -- CONFIRM the a2091 15C beam energy (placeholder default is 195).\033[0m" << std::endl;

   const double u = 931.49401;
   const double m_C15 = 15.0105993 * u, m_p = 1.007825 * u;  // beam 15C, target p
   const double m_eject = mEjectAmu * u, m_resid = mResidAmu * u;
   if (chanLabel.IsNull())
      chanLabel = (mEjectAmu > 1.5) ? "15C(p,d)14C" : "15C(p,p')";
   printf("channel: %s   fitter: %s   (m3 = %.6f u, m4 = %.6f u)\n", chanLabel.Data(), fitter.Data(),
          mEjectAmu, mResidAmu);

   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a2091/UKF/pp/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   TFile *fcache = new TFile((plotDir + "proton_kin" + outTag + ".root").Data(), "RECREATE");
   TNtuple *ntk = new TNtuple("pk", "proton kinematics", "ke:theta:vertexz:thcm:ex:chi2ndf");

   TH1F *hex = new TH1F("hex", Form("%s excitation energy;E_{x} [MeV];ejectiles", chanLabel.Data()), 200, -5, 25);
   TH2F *hexcm = new TH2F("hexcm", "E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -5, 25);
   TH2F *hket = new TH2F("hket", "proton KE vs #theta_{lab} (handedness check);#theta_{lab} [deg];KE [MeV]", 180, 0,
                         180, 200, 0, 40);
   TH2F *hexvz =
      new TH2F("hexvz", "E_{x} vs vertex z;vertex z [mm];E_{x} [MeV]", 200, -100, 1100, 150, -5, 25);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nProton = 0, nEvtWithTrk = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString uf = inDir + run + "_" + fitter + ".root";
      if (gSystem->AccessPathName(uf)) {
         printf("skip %s (no %s)\n", run.Data(), uf.Data());
         continue;
      }
      TFile *fu = TFile::Open(uf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      TClonesArray *te = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      Long64_t N = tu->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         tu->GetEntry(i);
         if (te->GetEntries() == 0)
            continue;
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev)
            continue;
         bool counted = false;
         for (auto &ft : ev->GetFittedTracks()) {
            if (!ft)
               continue;
            auto &k = ft->GetKinematics();
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double c2n = ndf > 0 ? chi2 / ndf : 1e9;
            double ke = k.kineticEnergy, thRad = k.theta;
            if (ke <= 0 || ke > 1000 || c2n > chi2Cut)
               continue;
            if (!counted) {
               ++nEvtWithTrk;
               counted = true;
            }
            double thDeg = thRad * TMath::RadToDeg();
            hket->Fill(thDeg, ke);
            auto [ex, thcm] = kine_2b(m_C15, m_p, m_eject, m_resid, Ebeam, thRad, ke);
            if (std::isnan(ex))
               continue;
            ++nProton;
            double vz = ft->GetVertex(0).Z(); // reaction vertex z [mm] -- was hardcoded 0
            hex->Fill(ex);
            hexcm->Fill(thcm, ex);
            hexvz->Fill(vz, ex);
            ntk->Fill(ke, thDeg, vz, thcm, ex, c2n);
         }
      }
      fu->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nevents with a good track: %ld | protons -> Ex: %ld\n", nEvtWithTrk, nProton);
   fcache->cd();
   ntk->Write();

   TCanvas *c = new TCanvas("c", "ex", 1550, 520);
   c->Divide(3, 1);
   c->cd(1);
   hex->Draw();
   c->cd(2);
   hket->Draw("colz");
   c->cd(3);
   hexcm->Draw("colz");
   TString png = plotDir + "ex_C15" + outTag + ".png";
   c->SaveAs(png);
   printf("saved %s\n", png.Data());
}
