/// @file ex_slices_3Hed.C
/// @brief Pre-fit Ex(47K) sliced in theta_lab, with the three states overlaid in every slice.
///
/// The question this answers: is there an angular window where the 0.36 MeV state separates from
/// the ground state? The angle-integrated spectrum says no (median +0.049 vs +0.337, IQR ~2.4 MeV
/// each), but the resolution is strongly angle dependent -- on the FITTED ground-state sample the
/// Ex IQR runs 3.39 MeV over theta_lab 55-90 and 0.84 MeV over 100-140. If that trend holds per
/// state, the useful measurement lives in a window, not over the whole range.
///
/// PRE-FIT, from the AtSpyralPID circle at fMinPoints = 15, over all six samples. It is the
/// pre-fit version rather than the fitted one because only gs_s3001 has been fitted so far, and
/// comparing states requires all three. The machinery is identical to kinematics_3Hed.C -- same
/// masses, same mirrored vertex, same two-body inversion -- so the slices are directly comparable
/// to the integrated numbers that macro prints.
///
/// Each panel is AREA-NORMALISED per state. The three states carry similar statistics overall but
/// not slice by slice, and the question here is about shape and position, not yield.
///
///   root -b -q 'ex_slices_3Hed.C()'

void ex_slices_3Hed(TString dir = "/mnt/f/ar46_3hed",
                    TString tags = "gs_s3001,gs_s3002,360_s3011,360_s3012,2020_s3021,2020_s3022",
                    TString png = "plots/ex_slices_3Hed.png", Int_t minPoints = 15, Double_t bField = 2.85,
                    Double_t dThetaMax = 10.0, Double_t driftLength = 100.0)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double M_b = 42809.757, M_t = 2808.392, M_R = 43734.759, M_e = 1875.613;
   const double Tb0 = 598.0, dEdz = 0.957, c = 299.792458;

   // Slice edges follow the structure, not round numbers: 83-96 is the band AtSpyralPID cannot
   // reconstruct at all (it is its own slice so the gap is visible rather than diluted), and the
   // 100-140 slice is the one the fitted sample showed to be 4x better than the forward half.
   const double edges[] = {58, 70, 83, 96, 110, 140};
   const int NS = 5;
   const char *stateName[3] = {"g.s. 1/2^{+}", "0.36 MeV 3/2^{+}", "2.02 MeV 7/2^{-}"};
   const double stateEx[3] = {0.0, 0.360, 2.020};
   const int stateCol[3] = {kBlack, kAzure + 2, kRed + 1};

   TH1D *h[3][NS];
   for (int s = 0; s < 3; ++s)
      for (int k = 0; k < NS; ++k) {
         h[s][k] = new TH1D(TString::Format("h_%d_%d", s, k), "", 90, -4, 6);
         h[s][k]->SetLineColor(stateCol[s]);
         h[s][k]->SetLineWidth(2);
      }

   AtTools::AtSpyralPID spy;
   spy.SetBField(std::abs(bField));
   if (minPoints > 0)
      spy.SetMinPoints(minPoints);

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      if (tg.IsNull())
         continue;
      int si = tg.BeginsWith("gs") ? 0 : (tg.BeginsWith("360") ? 1 : 2);
      TString fs = dir + "/" + tg + "_sim.root", fr = dir + "/" + tg + "_reco.root";
      if (gSystem->AccessPathName(fs) || gSystem->AccessPathName(fr)) {
         printf("  skip %-12s (missing)\n", tg.Data());
         continue;
      }
      TFile *Fs = TFile::Open(fs), *Fr = TFile::Open(fr);
      TTree *ts = (TTree *)Fs->Get("cbmsim"), *tr = (TTree *)Fr->Get("cbmsim");
      TClonesArray *mc = nullptr, *pa = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pa);

      Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
      for (Long64_t i = 0; i < N; ++i) {
         ts->GetEntry(i);
         tr->GetEntry(i);
         double thTrue = -1;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020)
               continue;
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp > 0)
               thTrue = std::acos(pz / pp) * TMath::RadToDeg();
            break;
         }
         if (thTrue < 0 || !pa || !pa->GetEntriesFast())
            continue;
         auto *pe = (AtPatternEvent *)pa->At(0);
         if (!pe)
            continue;
         double bd = 1e9, bTh = 0, bT = 0, bZ = 0;
         bool got = false;
         for (auto &track : pe->GetTrackCand()) {
            auto res = spy.Estimate(const_cast<AtTrack &>(track));
            if (!res.valid)
               continue;
            double th = 180.0 - res.polar * TMath::RadToDeg();
            double p = c * res.brho;
            double d = std::fabs(th - thTrue);
            if (d < bd) {
               bd = d;
               bTh = th;
               bT = std::sqrt(p * p + M_e * M_e) - M_e;
               bZ = res.vertex.Z() / 10.0;
               got = true;
            }
         }
         if (!got || bd > dThetaMax)
            continue;
         int sl = -1;
         for (int k = 0; k < NS; ++k)
            if (bTh >= edges[k] && bTh < edges[k + 1])
               sl = k;
         if (sl < 0)
            continue;

         double Tb = Tb0 - dEdz * (driftLength - bZ); // vertex is mirrored, as everywhere here
         if (Tb < 50 || Tb > Tb0 + 20)
            continue;
         double Eb = Tb + M_b, pb = std::sqrt(Tb * (Tb + 2 * M_b));
         double Ed = bT + M_e, pd = std::sqrt(bT * (bT + 2 * M_e));
         double th = bTh * TMath::DegToRad();
         double ER = Eb + M_t - Ed;
         double pRz = pb - pd * std::cos(th), pRt = pd * std::sin(th);
         double m2 = ER * ER - pRz * pRz - pRt * pRt;
         if (m2 <= 0)
            continue;
         h[si][sl]->Fill(std::sqrt(m2) - M_R);
      }
      Fs->Close();
      Fr->Close();
      printf("  %-12s done\n", tg.Data());
   }
   delete ta;

   printf("\n  theta_lab slice     state        n    median     IQR   (true Ex)\n");
   for (int k = 0; k < NS; ++k) {
      for (int s = 0; s < 3; ++s) {
         TH1D *hh = h[s][k];
         if (hh->GetEntries() < 50) {
            printf("   %3.0f - %3.0f deg      %-12s %5.0f      --       --\n", edges[k], edges[k + 1],
                   s == 0 ? "gs" : (s == 1 ? "0.36" : "2.02"), hh->GetEntries());
            continue;
         }
         double qp[3] = {0.25, 0.50, 0.75}, qv[3];
         hh->GetQuantiles(3, qv, qp);
         printf("   %3.0f - %3.0f deg      %-12s %5.0f  %+7.3f  %6.3f   (%.2f)\n", edges[k], edges[k + 1],
                s == 0 ? "gs" : (s == 1 ? "0.36" : "2.02"), hh->GetEntries(), qv[1], qv[2] - qv[0], stateEx[s]);
      }
      printf("\n");
   }

   TCanvas *cv = new TCanvas("cS", "Ex per angle slice", 1500, 900);
   cv->Divide(3, 2);
   for (int k = 0; k < NS; ++k) {
      cv->cd(k + 1);
      double mx = 0;
      for (int s = 0; s < 3; ++s)
         if (h[s][k]->Integral() > 0) {
            h[s][k]->Scale(1.0 / h[s][k]->Integral());
            mx = std::max(mx, h[s][k]->GetMaximum());
         }
      for (int s = 0; s < 3; ++s) {
         h[s][k]->SetTitle(TString::Format("#theta_{lab} %.0f - %.0f deg;E_{x} [MeV];fraction / bin", edges[k],
                                           edges[k + 1]));
         h[s][k]->SetMaximum(1.25 * mx);
         h[s][k]->Draw(s == 0 ? "hist" : "hist same");
      }
      for (int s = 0; s < 3; ++s) {
         auto *ln = new TLine(stateEx[s], 0, stateEx[s], 1.25 * mx);
         ln->SetLineColor(stateCol[s]);
         ln->SetLineStyle(3);
         ln->Draw();
      }
   }
   cv->cd(6);
   auto *lg = new TLegend(0.1, 0.35, 0.9, 0.75);
   lg->SetBorderSize(0);
   lg->SetHeader("dotted vertical = true excitation");
   for (int s = 0; s < 3; ++s)
      lg->AddEntry(h[s][0], stateName[s], "l");
   lg->Draw();
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   cv->SaveAs(png);
   printf("  wrote %s\n\n", png.Data());
}
