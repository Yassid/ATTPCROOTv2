/// @file decompose_sigmax.C
/// @brief Decompose the digi hit position resolution sigma_x (~0.83 mm, the thing
///        that limits momentum). Match each reconstructed hit to its nearest MC
///        point (xy) and split the residual into RADIAL (the bending/momentum
///        direction) and AZIMUTHAL, then compare each to the local PUMA pad size
///        (equal-area annular map: 16 rings x 256 pads). If sigma ~ pad/sqrt(12)
///        it is pad-quantization-limited; if larger, PSA/diffusion adds on top.
/// Run: root -b -q decompose_sigmax.C
double iqr(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }
double rms(const std::vector<double>&v){ if(v.empty())return 0; double m=0,s=0; for(double x:v)m+=x; m/=v.size(); for(double x:v)s+=(x-m)*(x-m); return std::sqrt(s/v.size()); }
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

void decompose_sigmax(TString digiFile = "./data/output_digi_pi.root", TString simFile = "./data/attpcsim.root")
{
   gSystem->Load("libAtReconstruction.so");
   const double kRin = 62.9, kRout = 121.1; const int kNring = 16, kNpad = 256;
   // equal-area ring boundaries: r_{i+1}^2 = r_i^2 + (Rout^2-Rin^2)/Nring
   std::vector<double> redge(kNring + 1); redge[0] = kRin;
   double dA = (kRout * kRout - kRin * kRin) / kNring;
   for (int i = 0; i < kNring; ++i) redge[i + 1] = std::sqrt(redge[i] * redge[i] + dA);
   auto ringOf = [&](double r) { for (int i = 0; i < kNring; ++i) if (r >= redge[i] && r < redge[i + 1]) return i; return kNring - 1; };

   TFile fD(digiFile); TTree *tD = (TTree *)fD.Get("cbmsim");
   TFile fS(simFile);  TTree *tS = (TTree *)fS.Get("cbmsim");
   TClonesArray *pe = new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent", &pe);
   TClonesArray *mcP = new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint", &mcP);

   std::vector<double> dRad, dAzi, dTot, padRad, padAzi;
   Long64_t nE = std::min(tD->GetEntries(), tS->GetEntries());
   for (Long64_t e = 0; e < nE; ++e) {
      tD->GetEntry(e); tS->GetEntry(e);
      if (!pe->GetEntries()) continue;
      int nMC = mcP->GetEntries();
      std::vector<double> mx(nMC), my(nMC);
      for (int k = 0; k < nMC; ++k) { auto *mp = (AtMCPoint *)mcP->At(k); mx[k] = mp->GetX() * 10; my[k] = mp->GetY() * 10; }
      for (auto &tr : ((AtPatternEvent *)pe->At(0))->GetTrackCand()) {
         for (auto &h : tr.GetHitArray()) {
            auto p = h->GetPosition();
            double hx = p.X(), hy = p.Y();
            // nearest MC point in xy
            double best = 25.0; int bk = -1; // 5mm gate
            for (int k = 0; k < nMC; ++k) { double d2 = (hx - mx[k]) * (hx - mx[k]) + (hy - my[k]) * (hy - my[k]); if (d2 < best) { best = d2; bk = k; } }
            if (bk < 0) continue;
            double dx = hx - mx[bk], dy = hy - my[bk];
            double r = std::hypot(hx, hy); if (r < 1) continue;
            double rhx = hx / r, rhy = hy / r;        // radial unit
            double dr = dx * rhx + dy * rhy;            // radial residual
            double da = -dx * rhy + dy * rhx;           // azimuthal residual
            dRad.push_back(dr); dAzi.push_back(da); dTot.push_back(std::hypot(dx, dy));
            int ri = ringOf(r);
            padRad.push_back(redge[ri + 1] - redge[ri]);         // ring radial width
            padAzi.push_back(2 * M_PI * r / kNpad);              // azimuthal pitch at r
         }
      }
   }
   double sr = iqr(dRad), sa = iqr(dAzi), st = iqr(dTot);
   double prad = med(padRad), pazi = med(padAzi);
   printf("\n===== sigma_x decomposition (%zu matched hits) =====\n", dRad.size());
   printf("  RADIAL   residual sigma (IQR) : %.2f mm   [ring width %.2f mm -> pad/sqrt12 = %.2f mm]\n", sr, prad, prad / std::sqrt(12.));
   printf("  AZIMUTHAL residual sigma      : %.2f mm   [az pitch  %.2f mm -> pad/sqrt12 = %.2f mm]\n", sa, pazi, pazi / std::sqrt(12.));
   printf("  TOTAL |xy| residual sigma     : %.2f mm   (this is the sigma_x the fit sees)\n", st);
   printf("  medians: radial %+.2f, azimuthal %+.2f mm\n", med(dRad), med(dAzi));
   printf("\n  Radial is the momentum-sensitive direction. If radial sigma ~ ring/sqrt12, it is\n");
   printf("  pad-quantization-limited (needs finer radial segmentation); if >> that, PSA/diffusion.\n");
}
