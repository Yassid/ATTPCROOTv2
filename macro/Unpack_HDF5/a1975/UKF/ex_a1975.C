/// @file ex_a1975.C
/// @brief 16C excitation-energy spectrum from the clean UKF protons (16C(p,p')16C*).
///
/// Uses the UKF-fitted proton (KE, theta) + two-body kinematics to reconstruct the
/// 16C excitation energy, on the IC-gated 16C beam + proton-PID-gated + good-fit
/// sample (same selection as ukf_clean_a1975.C). Elastic -> Ex~0; 16C excited states
/// appear as peaks. Beam energy is the nominal 192 MeV (no vertex energy-loss
/// correction yet -- that refinement broadens/shifts Ex; flagged for later).
///
///   root -b -q 'ex_a1975.C("run_0106,...","/mnt/f/a1975/reco/")'

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

void ex_a1975(TString runsCSV = "run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113,run_0114,"
                                "run_0115",
             TString inDir = "/mnt/f/a1975/reco/", TString gateFile = "proton_band.json", double Ebeam = 192.0,
             double icMin = 950, double icMax = 1350, int icTbLo = 1000, int icTbHi = 1350, double chi2Cut = 5.0,
             double bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double u = 931.49401;
   const double m_C16 = 16.0147 * u, m_p = 1.007825 * u;

   auto pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   // flat cache of clean protons so the vertex correction / resolution work can be
   // iterated WITHOUT re-reading the slow FRIB files.
   TFile *fcache = new TFile("proton_kin.root", "RECREATE");
   TNtuple *ntk = new TNtuple("pk", "clean proton kinematics", "ke:theta:vertexz:thcm:ex:chi2ndf");

   TH1F *hex = new TH1F("hex", "16C excitation energy (p,p');E_{x} [MeV];protons", 200, -5, 25);
   TH2F *hexcm = new TH2F("hexcm", "E_{x} vs #theta_{cm};#theta_{cm} [deg];E_{x} [MeV]", 180, 0, 180, 150, -5, 25);
   TH2F *hexvz = new TH2F("hexvz", "E_{x} vs vertex z (uncorrected);vertex z [mm];E_{x} [MeV]", 200, 0, 1000, 150, -5,
                          25);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nProton = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString uf = inDir + run + "_ukf.root", ff = inDir + run + "_FRIB.root";
      if (gSystem->AccessPathName(uf) || gSystem->AccessPathName(ff)) {
         printf("skip %s\n", run.Data());
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
            continue; // IC gate (16C)
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
            if (!r.valid || !pid.IsInside(r.sqrtdEdx, r.brho))
               continue;
            auto [ex, thcm] = kine_2b(m_C16, m_p, m_p, m_C16, Ebeam, thRad, ke);
            if (std::isnan(ex))
               continue;
            ++nProton;
            hex->Fill(ex);
            hexcm->Fill(thcm, ex);
            hexvz->Fill(r.vertex.Z(), ex);
            ntk->Fill(ke, thRad * TMath::RadToDeg(), r.vertex.Z(), thcm, ex, c2n);
         }
      }
      fu->Close();
      fc->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nclean protons -> Ex: %ld\n", nProton);
   fcache->cd();
   ntk->Write();
   printf("cached -> proton_kin.root (ntuple pk: ke:theta:vertexz:thcm:ex:chi2ndf)\n");

   TCanvas *c = new TCanvas("c", "ex", 1550, 520);
   c->Divide(3, 1);
   c->cd(1);
   hex->Draw();
   c->cd(2);
   hexcm->Draw("colz");
   c->cd(3);
   hexvz->Draw("colz");
   // profile: mean Ex per vertex-z slice (the systematic to flatten)
   TProfile *prof = hexvz->ProfileX("prof");
   prof->SetLineColor(kRed);
   prof->SetLineWidth(2);
   prof->Draw("same");
   c->SaveAs("ex_spectrum.png");
   printf("saved ex_spectrum.png\n");
}
