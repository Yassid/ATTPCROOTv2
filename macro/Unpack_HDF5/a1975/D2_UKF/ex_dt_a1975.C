/// @file ex_dt_a1975.C
/// @brief 15C excitation-energy spectrum from TRITON-hypothesis genfit fits — the
///        16C(d,t)15C neutron-pickup channel (a1975 D2-target runs 0016-0103).
///
/// Inverse kinematics: 16C beam + d target -> t + 15C. The triton is FORWARD in the
/// lab: the g.s. locus runs from (theta~50 deg, KE~3 MeV) up to (theta~25 deg,
/// KE~70 MeV), turning over at theta_max ~ 56 deg. Sibling of ex_dp_a1975.C
/// (the (d,p) proton channel); reads reco_d2/<run>_genfitter_t<suffix>.root.
///
///   root -b -q 'ex_dt_a1975.C("run_0016","/mnt/f/a1975/reco_d2/")'
///
/// This macro exists to REPLACE the throw-away cache that produced ~/a1975_panels/
/// dt_kin.root, which stored only (ke, theta, chi2ndf, ex) — with no vertex, no
/// PID and no run number there is no way to tell a mis-identified deuteron from a
/// badly fitted triton. Everything needed for that call now goes into the ntuple.
///
/// The old (and reportedly much better) analysis of this channel was
/// a1975/C16_dt_anaFit.C. Its selection chain, for reference:
///   - ion-chamber beam gate 900 < E_IC < 1300 (FRIB) -- NOT available in reco_d2
///   - one track per event: the largest theta among tracks with theta < 90 deg
///   - hand-drawn TCutG on (eLossADC, brho), not on the Spyral sqrt(dEdx) plane
///   - 25 < theta < 56 deg, 2 < KE_t < 20 MeV, 0 < vertex z < 50 cm
///   - Ebeam 180 MeV (192 is commented out in that macro)
///   - deuteron-hypothesis fit with KE rescaled by m_d/m_t = 2/3, not a triton fit
/// The defaults here follow the NEW pipeline; pass old-style cuts to compare.

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
   double arg = (std::cos(thetalab) * omega2(s, m1 * m1, m2 * m2) * omega2(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) +
                s + u - m2 * m2;
   if (arg < 0)
      return {std::nan(""), std::nan("")};
   double m4_ex = std::sqrt(arg);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (omega2(s, m1 * m1, m2 * m2) * omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

/// icDir: directory of the per-run ion-chamber files written by unpackIC_d2.C. The IC
/// amplitude is STORED in the ntuple (branch `ic`), never cut on here, so the beam gate
/// stays a knob at analysis time. ic = -1 means "no IC file for this run".
void ex_dt_a1975(TString runsCSV = "run_0016", TString inDir = "/mnt/f/a1975/reco_d2/", TString fitSuffix = "",
                 TString gateFile = "pid/triton_d2.json", double Ebeam = 180.0, double chi2Cut = 10.0,
                 double exShift = 0.0, double keMax = 90.0, double thMinDeg = 10.0, double thMaxDeg = 90.0,
                 TString cacheOut = "triton_kin_dt.root", TString plotOut = "plots/ex_dt_spectrum.png",
                 TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.0135532 * u;
   const double m_t = 3.01550072 * u, m_C15 = 15.0105993 * u;

   bool useGate = (gateFile.Length() > 0);
   AtTools::AtParticleID pid;
   if (useGate)
      pid = AtTools::AtParticleID::LoadJSON(gateFile.Data());

   TFile *fcache = new TFile(cacheOut, "RECREATE");
   // ke/theta/ex are the BACK-EXTRAPOLATED values (AtGenfitter's Xtr slot: vertex state on the
   // beam axis, plus any manual dE/dx). kefit/thetafit/exfit are what the fit itself returned at
   // the first measurement point, kept so the two can be compared on the same tracks. Before
   // AtGenfitter wrote them separately both slots held the corrected value, so an older cache
   // has kefit == ke by construction.
   TNtuple *ntk = new TNtuple("pk", "candidate triton kinematics (d,t)",
                              "ke:theta:vertexz:vertexr:thcm:ex:chi2ndf:brho:dedx:sqrtdedx:ncl:ntrk:run:ic:"
                              "kefit:thetafit:exfit");

   auto *hex = new TH1F("hex", "^{15}C excitation energy (d,t), triton hyp;E_{x}(^{15}C) [MeV];tritons", 200, -10, 25);
   auto *hexAll = new TH1F("hexAll", "Ex, no PID gate", 200, -10, 25);
   hexAll->SetLineColor(kGray + 2);
   auto *hkt = new TH2F("hkt", "KE vs #theta_{lab} (gated);#theta_{lab} [deg];KE_{t} [MeV]", 120, 0, 90, 200, 0, 90);
   auto *hpidPk = new TH2F("hpidPk", "PID plane, |E_{x}|<3;#sqrt{dEdx};B#rho [T m]", 200, 0, 40, 200, 0, 2.5);
   auto *hpidBl = new TH2F("hpidBl", "PID plane, E_{x}>10;#sqrt{dEdx};B#rho [T m]", 200, 0, 40, 200, 0, 2.5);
   auto *hvzPk = new TH1F("hvzPk", "vertex z;z [mm];", 120, -100, 1200);
   auto *hvzBl = new TH1F("hvzBl", "vertex z", 120, -100, 1200);
   hvzBl->SetLineColor(kRed);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nFit = 0, nValPID = 0, nGated = 0, nMatch = 0, nCand = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString gf = inDir + run + "_genfitter_t" + fitSuffix + ".root";
      if (gSystem->AccessPathName(gf)) {
         printf("skip %s (missing %s)\n", run.Data(), gf.Data());
         continue;
      }
      // The fit file may carry a stage suffix (run_0017_multifit), but the run number and the IC
      // file are keyed on the bare run_NNNN. Taking the LAST four characters silently returned
      // runNo = 0 and looked for run_0017_multifit_IC.root, so the whole cache came out with
      // ic = -1 -- i.e. no beam gate at all, which is worth ~1.7 MeV on the reconstructed Ex.
      TString runTag = run;
      Ssiz_t rp = runTag.Index("run_");
      if (rp != kNPOS && runTag.Length() >= rp + 8)
         runTag = runTag(rp, 8); // "run_NNNN"
      int runNo = TString(runTag(runTag.Length() - 4, 4)).Atoi();
      TFile *fu = TFile::Open(gf);
      TTree *tu = (TTree *)fu->Get("cbmsim");
      TClonesArray *te = nullptr, *pe = nullptr;
      tu->SetBranchAddress("AtTrackingEvent", &te);
      tu->SetBranchAddress("AtPIDEvent", &pe);

      Long64_t N = tu->GetEntries();

      // ion chamber, entry-aligned with the fit tree (see unpackIC_d2.C on the alignment)
      std::vector<float> icv;
      if (icDir.Length()) {
         TString icf = icDir + runTag + "_IC.root";
         // RETRY. On WSL the data drive is drvfs, and with several of these running in parallel
         // AccessPathName/TFile::Open fail INTERMITTENTLY on a file that is plainly there: the
         // same run reports "IC for 41610 events" alone and "no IC file" inside a 6-way batch.
         // The old code took that at face value and fell through to ic = -1, i.e. it silently
         // produced an UNGATED cache -- 27 of 47 runs in one build, which then looks like a
         // physics difference in the track counts. Never treat a missing IC as routine.
         TFile *fi = nullptr;
         for (int att = 0; att < 5 && !fi; ++att) {
            if (att)
               gSystem->Sleep(200);
            if (!gSystem->AccessPathName(icf))
               fi = TFile::Open(icf);
            if (fi && fi->IsZombie()) {
               fi->Close();
               fi = nullptr;
            }
         }
         if (fi) {
            TTree *ti = (TTree *)fi->Get("ic");
            float icm = 0;
            ti->SetBranchAddress("icmax", &icm);
            icv.reserve(ti->GetEntries());
            for (Long64_t j = 0; j < ti->GetEntries(); ++j) {
               ti->GetEntry(j);
               icv.push_back(icm);
            }
            fi->Close();
            printf("  %s: IC for %zu events\n", run.Data(), icv.size());
         } else {
            printf("\n  *** WARNING %s: NO IC after 5 attempts (%s) -- this run is UNGATED.\n"
                   "      Every track gets ic = -1 and the beam gate then drops all of them.\n"
                   "      Do not merge this into a cache that other runs are gated in.\n\n",
                   run.Data(), icf.Data());
         }
      }
      for (Long64_t i = 0; i < N; ++i) {
         float icNow = (i < (Long64_t)icv.size()) ? icv[i] : -1.f;
         tu->GetEntry(i);
         if (pe->GetEntries() == 0 || te->GetEntries() == 0)
            continue;
         auto *pidev = (AtPIDEvent *)pe->At(0);
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!pidev || !ev)
            continue;
         std::map<int, AtFittedTrack *> fmap;
         for (auto &ft : ev->GetFittedTracks())
            if (ft)
               fmap[ft->GetTrackID()] = ft.get();
         nFit += fmap.size();
         int ntrk = fmap.size();
         for (auto &sr : pidev->GetSpyral()) {
            if (!sr.valid)
               continue;
            ++nValPID;
            bool inGate = !useGate || pid.IsInside(sr.sqrtdEdx, sr.brho);
            auto it = fmap.find(sr.trackID);
            if (it == fmap.end())
               continue;
            ++nMatch;
            if (inGate)
               ++nGated;
            auto *ft = it->second;
            auto &k = ft->GetKinematicsXtr(); // back-extrapolated: the value the analysis uses
            auto &kf = ft->GetKinematics();   // raw fit at the first measurement point
            double ndf = ft->GetTrackMetadata()->GetNdf(), chi2 = ft->GetTrackMetadata()->GetChi2();
            double c2n = ndf > 0 ? chi2 / ndf : 1e9;
            double ke = k.kineticEnergy, thRad = k.theta, thDeg = thRad * TMath::RadToDeg();
            double keF = kf.kineticEnergy, thRadF = kf.theta;
            if (ke <= 0 || ke > keMax || c2n > chi2Cut)
               continue;
            if (thDeg < thMinDeg || thDeg > thMaxDeg)
               continue;
            auto [ex, thcm] = kine_2b(m_C16, m_d, m_t, m_C15, Ebeam, thRad, ke);
            if (std::isnan(ex))
               continue;
            double exCal = ex + exShift;
            hexAll->Fill(exCal);
            if (!inGate)
               continue;
            ++nCand;
            auto v = ft->GetVertex();
            double vr = std::sqrt(v.X() * v.X() + v.Y() * v.Y());
            hex->Fill(exCal);
            hkt->Fill(thDeg, ke);
            if (std::fabs(exCal) < 3) {
               hpidPk->Fill(sr.sqrtdEdx, sr.brho);
               hvzPk->Fill(v.Z());
            } else if (exCal > 10) {
               hpidBl->Fill(sr.sqrtdEdx, sr.brho);
               hvzBl->Fill(v.Z());
            }
            // same kinematics, evaluated on the uncorrected fit value; NaN if it does not close
            auto [exF, thcmF] = kine_2b(m_C16, m_d, m_t, m_C15, Ebeam, thRadF, keF);
            float row[17] = {(float)ke,
                             (float)thDeg,
                             (float)v.Z(),
                             (float)vr,
                             (float)thcm,
                             (float)exCal,
                             (float)c2n,
                             (float)sr.brho,
                             (float)sr.dEdx,
                             (float)sr.sqrtdEdx,
                             (float)sr.nClusters,
                             (float)ntrk,
                             (float)runNo,
                             icNow,
                             (float)keF,
                             (float)(thRadF * TMath::RadToDeg()),
                             (float)(std::isnan(exF) ? -999.0 : exF + exShift)};
            ntk->Fill(row);
         }
      }
      fu->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nfitted tracks %ld | valid Spyral %ld | matched to a fit %ld | in triton gate %ld | candidates %ld\n", nFit,
          nValPID, nMatch, nGated, nCand);
   fcache->cd();
   ntk->Write();
   printf("cached -> %s\n", cacheOut.Data());

   // g.s. kinematic locus (theta_lab, KE) of the triton, from the old kinematics table
   TString kinFile = TString(gSystem->Getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1975/C16_dt_gs.txt";
   std::vector<double> gx, gy;
   std::ifstream in(kinFile.Data());
   double a, b, c, d, e;
   while (in >> a >> b >> c >> d >> e) {
      gx.push_back(b);
      gy.push_back(c);
   }
   TGraph *gr = gx.empty() ? nullptr : new TGraph(gx.size(), gx.data(), gy.data());
   if (gr) {
      gr->SetLineColor(kRed);
      gr->SetLineWidth(2);
   }

   auto *cv = new TCanvas("cv", "ex_dt", 1600, 1100);
   cv->Divide(3, 2);
   cv->cd(1);
   hexAll->Draw();
   hex->Draw("same");
   cv->cd(2);
   gPad->SetLogz();
   hkt->Draw("colz");
   if (gr)
      gr->Draw("l same");
   cv->cd(3);
   gPad->SetLogz();
   hpidPk->Draw("colz");
   cv->cd(4);
   gPad->SetLogz();
   hpidBl->Draw("colz");
   cv->cd(5);
   hvzPk->Draw();
   hvzBl->Draw("same");
   cv->cd(6);
   hex->Draw();
   cv->SaveAs(plotOut);
   printf("saved %s\n", plotOut.Data());
}
