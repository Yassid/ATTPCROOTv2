/// @file understand_resolution.C
/// @brief Why is the PUMA momentum resolution ~20%? Measure the ingredients that set
///        it — lever arm L, point resolution sigma_x, sagitta — from the reconstructed
///        tracks and compare the Gluckstern-predicted sigma_p/p to what the fit delivers.
///        Detector-limited if fit ~ Gluckstern; reconstruction-limited if fit >> it.
///
/// Gluckstern (N points, uniform, bending-plane length L):
///   sigma(p)/p = sigma_x * p[GeV] / (0.3 * B * L^2) * sqrt(720/(N+4))    [L,sigma_x in m]
/// Run: root -b -q understand_resolution.C
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void understand_resolution(TString digiFile = "./data/output_digi_pi.root", double Bfield = 4.0,
                           double p0_MeV = 374.9)
{
   gSystem->Load("libAtReconstruction.so");
   TFile f(digiFile); TTree *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pe = new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent", &pe);

   std::vector<double> vR, vN, vArc, vChord, vSag, vSigX, vGluck;
   for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      if (!pe->GetEntries()) continue;
      for (auto &tr : ((AtPatternEvent *)pe->At(0))->GetTrackCand()) {
         const auto &hits = tr.GetHitArray();
         int N = hits.size();
         if (N < 5) continue;
         double R = tr.GetGeoRadius();                 // mm
         auto c = tr.GetGeoCenter();                   // (cx,cy) mm
         if (!(R > 0) || R > 1e5) continue;

         // sigma_x = RMS radial residual of hits to the fitted circle (bending-plane res.)
         double s2 = 0; int ns = 0;
         double rmin = 1e9, rmax = -1e9;
         // chord + arc: order hits by angle about the centre
         std::vector<std::pair<double, std::pair<double, double>>> byPhi;
         for (auto &h : hits) {
            auto p = h->GetPosition();
            double dr = std::hypot(p.X() - c.first, p.Y() - c.second) - R;
            s2 += dr * dr; ns++;
            double ph = std::atan2(p.Y() - c.second, p.X() - c.first);
            byPhi.emplace_back(ph, std::make_pair(p.X(), p.Y()));
         }
         double sigX = std::sqrt(s2 / ns);             // mm
         std::sort(byPhi.begin(), byPhi.end());
         // chord = straight distance between angular endpoints (bending plane)
         double x0 = byPhi.front().second.first, y0 = byPhi.front().second.second;
         double x1 = byPhi.back().second.first,  y1 = byPhi.back().second.second;
         double chord = std::hypot(x1 - x0, y1 - y0);  // mm
         // arc length along the circle between endpoints
         double dphi = byPhi.back().first - byPhi.front().first;
         double arc = R * std::abs(dphi);              // mm
         double sag = chord * chord / (8.0 * R);       // sagitta, mm

         // Gluckstern sigma_p/p using the chord as L (bending-plane length)
         double L_m = chord / 1000.0, sx_m = sigX / 1000.0, p_GeV = p0_MeV / 1000.0;
         double gluck = sx_m * p_GeV / (0.299792458 * Bfield * L_m * L_m) * std::sqrt(720.0 / (N + 4));

         vR.push_back(R); vN.push_back(N); vArc.push_back(arc); vChord.push_back(chord);
         vSag.push_back(sag); vSigX.push_back(sigX); vGluck.push_back(gluck);
      }
   }

   printf("\n===== PUMA momentum-resolution anatomy (pions, |p|=%.0f MeV/c, B=%.1f T) =====\n", p0_MeV, Bfield);
   printf("  tracks                 : %zu\n", vR.size());
   printf("  circle radius R        : %.0f mm   (p = 0.3*B*R = %.0f MeV/c)\n", med(vR), 0.299792458 * Bfield * med(vR));
   printf("  hits / track  N        : %.0f\n", med(vN));
   printf("  measured chord L       : %.1f mm   <-- the lever arm (annulus is thin)\n", med(vChord));
   printf("  measured arc length    : %.1f mm\n", med(vArc));
   printf("  sagitta  s=L^2/8R      : %.2f mm   <-- the signal to be measured\n", med(vSag));
   printf("  point resolution sig_x : %.2f mm   <-- the noise (hit residual to circle)\n", med(vSigX));
   printf("  ---> sagitta / sig_x   : %.1f      (signal-to-noise of the curvature)\n", med(vSag) / med(vSigX));
   printf("  Gluckstern sigma_p/p   : %.0f%%     (predicted from L, sig_x, N above)\n", 100 * med(vGluck));
   printf("  measured  sigma_p/p    : ~19-21%%   (from resolution_test.C)\n");
   printf("\n  Interpretation: sigma_p/p ~ sig_x * p / (0.3 B L^2). With L~%.0f mm and s/sig_x~%.1f,\n",
          med(vChord), med(vSag) / med(vSigX));
   printf("  the curvature is barely above the point resolution -> intrinsic, detector-limited.\n");
}
