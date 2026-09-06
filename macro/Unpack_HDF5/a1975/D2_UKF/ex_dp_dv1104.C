/// @file ex_dp_dv1104.C
/// @brief 17C excitation-energy cache from PROTON-hypothesis genfit fits -- the
///        16C(d,p)17C neutron-STRIPPING channel (a1975 D2-target runs 0016-0103).
///
/// Inverse kinematics: 16C beam + d target -> p + 17C. Unlike the (d,t) triton, which
/// is confined forward of theta_max ~ 56 deg, the proton is emitted over the FULL
/// angular range, and the transfer strength sits BACKWARD: g.s. Brho falls from
/// 1.02 T m at 20 deg to 0.19 T m at 150 deg (KE 48.6 -> 1.8 MeV). Two consequences
/// that do not arise in (d,t): backwardSeedFix must be on in the production, and the
/// low-KE end lands below genfit's beta*gamma = 0.05 floor (KE 1.17 MeV for a proton),
/// which is what proton_D2_300torr.txt is for.
///
///   root -b -q 'ex_dp_dv1104.C("run_0016_multifit","/mnt/f/a1975/gf_dp_cateloss/")'
///
/// THIS IS A PORT OF ex_dt_a1975.C AND THE `pk` TREE IS COLUMN-IDENTICAL TO IT, on
/// purpose: every (d,t) tool that reads a pk cache then works unchanged on this one.
/// Only the ejectile/residual masses and the defaults differ.
///
/// IT SUPERSEDES ex_dp_a1975.C, which cannot be used for a production cache: it
/// hardcodes the output name "proton_kin_dp.root" (so per-run parallel caching writes
/// 47 runs onto one file), stores no run number, no IC amplitude and no material-effects
/// provenance, cuts at chi2/ndf < 10 / KE < 50 / theta 10-170 by DEFAULT, and defaults
/// to Ebeam 192.0 and /mnt/f/a1975/reco_d2/ -- a purge-list energy and a deleted
/// directory. Its physics (the kinematics and the atomic mass set) was correct.
///
/// MASS CONVENTION. m1/m4 ATOMIC, m2/m3 NUCLEAR -- identical to ex_dt_a1975.C, and it
/// is not a free choice: the electrons must cancel side to side (6 per carbon each
/// side). Mixing the two is 0.511 MeV per electron and looks exactly like a calibration
/// offset. Checked: Q = -1.498 MeV against Sn(17C) - B(d) = 0.729 - 2.225 = -1.496.
/// (The same check found the bug in pid/gate_draw_dp.C's mass table on 2026-09-05.)

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
///
/// EVERY CUT DEFAULTS TO "OFF" HERE, unlike ex_dt_a1975.C, whose 10 / 90 / 10-90 defaults
/// are the ones its own cache script has to override. This is a NEW channel: no selection
/// is inherited, each is chosen on this data and applied downstream off the stored columns.
/// chi2Cut = 1e9 KEEPS collapsed fits (c2n is set to 1e9 when ndf <= 0, and 1e9 > 1e9 is
/// false); they are recovered downstream with chi2ndf < 1e8, which is why 1e9 must not be
/// used as a "no cut" value anywhere a cut is actually wanted.
void ex_dp_dv1104(TString runsCSV = "run_0016_multifit", TString inDir = "/mnt/f/a1975/gf_dp_cateloss/",
                 TString fitSuffix = "", TString gateFile = "pid/proton_d2_dv1104.json", double Ebeam = 184.25,
                 double chi2Cut = 1e9, double exShift = 0.0, double keMax = 1e9, double thMinDeg = 0.0,
                 double thMaxDeg = 180.0, TString cacheOut = "proton_kin_dp.root",
                 TString plotOut = "plots/ex_dp_spectrum.png", TString icDir = "/mnt/f/a1975/ic_d2/")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   const double u = 931.49401;
   const double m_C16 = 16.0147013 * u, m_d = 2.0135532 * u;      // ATOMIC 16C, NUCLEAR d
   const double m_p = 1.00727646688 * u, m_C17 = 17.0225864 * u;  // NUCLEAR p, ATOMIC 17C

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
   TNtuple *ntk = new TNtuple("pk", "candidate proton kinematics (d,p)",
                              "ke:theta:vertexz:vertexr:thcm:ex:chi2ndf:brho:dedx:sqrtdedx:ncl:ntrk:run:ic:"
                              "kefit:thetafit:exfit:matfx:matfb");
   // matfx/matfb are the MATERIAL-EFFECTS PROVENANCE of each fit and they are not cosmetic.
   // AtGenfitter defaults fMatEffectsFallback=kTRUE: when a material-effects fit throws (a
   // stopping triton's RK extrapolation does this often) the track is silently refitted with
   // setNoEffects(true) and kept. Without these columns a "material effects" production is a
   // BLEND of two different physics models in one spectrum, which inflates the width and
   // fakes a shift. Require matfx==1 to quote a pure matFX sample; matfb==1 marks the retries.

   auto *hex = new TH1F("hex", "^{17}C excitation energy (d,p), proton hyp;E_{x}(^{17}C) [MeV];protons", 200, -10, 25);
   auto *hexAll = new TH1F("hexAll", "Ex, no PID gate", 200, -10, 25);
   hexAll->SetLineColor(kGray + 2);
   // 0-180 in theta and 0-60 in KE. The (d,t) axes (0-90, 0-90) would CROP the (d,p) transfer
   // region outright: the strength is backward, where the g.s. proton is 1.8 MeV at 150 deg.
   auto *hkt = new TH2F("hkt", "KE vs #theta_{lab} (gated);#theta_{lab} [deg];KE_{p} [MeV]", 180, 0, 180, 200, 0, 60);
   // Brho 0-1.6, not the tritons' 0-2.5: the proton band lives at 0.05-1.03 and a 2.5 axis hides
   // it in the bottom sixth of the pad. Same reason yMax is 1.6 in the gate drawer.
   auto *hpidPk = new TH2F("hpidPk", "PID plane, |E_{x}|<3;#sqrt{dEdx};B#rho [T m]", 200, 0, 40, 200, 0, 1.6);
   auto *hpidBl = new TH2F("hpidBl", "PID plane, E_{x}>10;#sqrt{dEdx};B#rho [T m]", 200, 0, 40, 200, 0, 1.6);
   auto *hvzPk = new TH1F("hvzPk", "vertex z;z [mm];", 120, -100, 1200);
   auto *hvzBl = new TH1F("hvzBl", "vertex z", 120, -100, 1200);
   hvzBl->SetLineColor(kRed);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nFit = 0, nValPID = 0, nGated = 0, nMatch = 0, nCand = 0;
   for (int ri = 0; ri < runs->GetEntries(); ++ri) {
      TString run = ((TObjString *)runs->At(ri))->GetString();
      TString gf = inDir + run + "_genfitter_p" + fitSuffix + ".root";
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
            auto [ex, thcm] = kine_2b(m_C16, m_d, m_p, m_C17, Ebeam, thRad, ke);
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
            auto [exF, thcmF] = kine_2b(m_C16, m_d, m_p, m_C17, Ebeam, thRadF, keF);
            float row[19] = {(float)ke,
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
                             (float)(std::isnan(exF) ? -999.0 : exF + exShift),
                             (float)ft->GetTrackMetadata()->GetMatEffects(),
                             (float)ft->GetTrackMetadata()->GetMatEffectsFallback()};
            ntk->Fill(row);
         }
      }
      fu->Close();
      printf("processed %s\n", run.Data());
   }
   printf("\nfitted tracks %ld | valid Spyral %ld | matched to a fit %ld | in proton gate %ld | candidates %ld\n", nFit,
          nValPID, nMatch, nGated, nCand);
   fcache->cd();
   ntk->Write();
   printf("cached -> %s\n", cacheOut.Data());

   // g.s. kinematic locus (theta_lab, KE) of the PROTON. The (d,t) macro reads a precomputed
   // C16_dt_gs.txt; no such table exists for this channel, and computing it from the SAME masses
   // and Ebeam the Ex above was built with is better anyway -- a stale table is exactly how a
   // locus and a spectrum end up disagreeing without anyone noticing.
   std::vector<double> gx, gy;
   {
      const double m4x = m_C17; // Ex = 0
      double E1 = Ebeam + m_C16, sM = m_C16 * m_C16 + m_d * m_d + 2 * m_d * E1, rs = std::sqrt(sM);
      if (rs > m_p + m4x) {
         double pf = omega2(sM, m_p * m_p, m4x * m4x) / (2 * rs);
         double E3cm = std::sqrt(pf * pf + m_p * m_p);
         double beta = std::sqrt(E1 * E1 - m_C16 * m_C16) / (E1 + m_d);
         double gam = 1.0 / std::sqrt(1 - beta * beta);
         for (int i = 1; i < 1800; ++i) {
            double tcm = TMath::Pi() * i / 1800.0;
            double pz = gam * (pf * std::cos(tcm) + beta * E3cm), pt = pf * std::sin(tcm);
            double pp = std::hypot(pz, pt);
            gx.push_back(std::atan2(pt, pz) * TMath::RadToDeg());
            gy.push_back(std::sqrt(pp * pp + m_p * m_p) - m_p);
         }
      }
   }
   TGraph *gr = gx.empty() ? nullptr : new TGraph(gx.size(), gx.data(), gy.data());
   if (gr) {
      gr->SetLineColor(kRed);
      gr->SetLineWidth(2);
   }

   auto *cv = new TCanvas("cv", "ex_dp", 1600, 1100);
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
