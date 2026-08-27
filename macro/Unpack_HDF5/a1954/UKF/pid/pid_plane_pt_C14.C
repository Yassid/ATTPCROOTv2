/// @file pid_plane_pt_C14.C
/// @brief The a1954 14C PID plane over its FULL range, to look for a triton locus.
///
/// The (p,p) analysis draws its proton polygon on Brho 0.04-0.96 Tm and never looks above it.
/// 14C(p,t)12C puts the ejectile somewhere else entirely: Q = -4.64 MeV, and at
/// E(14C) = 159.75 MeV the triton emerges at theta_lab 8-24 deg with 14-56 MeV, giving
/// Brho ~ 1.5 Tm -- stiff enough that its helix radius (~52 cm) exceeds the chamber radius, so it
/// leaves rather than curls. Its dE/dx is about 1.7x a 20 MeV proton's, since dE/dx ~ Z^2/beta^2
/// and beta^2 is 0.024 against 0.041. So a triton locus should sit ABOVE AND TO THE RIGHT of the
/// proton one, and the question this macro answers is simply whether anything is there.
///
/// AtSpyralPID takes +2.85 here, the same as pipeline/gate_events_C14.C. The -2.85 in the
/// fitGenfit calls is genfit's convention and is a different parameter; passing it gives a
/// SIGNED Brho that is uniformly negative and an empty plot.
///
///   root -b -q 'pid_plane_pt_C14.C("run_0055")'
/// Peak counter, identical to the one in pipeline/gate_events_C14.C: peaks above thr with
/// hysteresis at thr/2. Exactly one pulse is required, which is what rejects pile-up.
static int ppCountPulses(const std::vector<Double_t> &adc, double thr, int tbLo, int tbHi)
{
   int n = 0; bool in = false;
   for (int i = tbLo; i <= tbHi && i < (int)adc.size(); ++i) {
      if (!in && adc[i] > thr) { in = true; ++n; }
      else if (in && adc[i] < 0.5 * thr) in = false;
   }
   return n;
}

void pid_plane_pt_C14(TString runsCSV = "run_0055",
                      TString inDir = "/mnt/f/a1954_C14_reco_hdb_slim/",
                      Double_t bField = +2.85, Long64_t maxEvt = -1,
                      Int_t icTbLo = 1050, Int_t icTbHi = 1250,
                      Double_t peakThr = 200, Int_t pkTbLo = 800, Int_t pkTbHi = 1500)
{
   gSystem->Load("libAtTools.so");
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   AtTools::AtSpyralPID spy;
   spy.SetBField(bField);

   // The per-track POINTS, not just the histograms: gate_draw_pt_C14.C needs one row per track to
   // draw a polygon and to test what falls inside it against the kinematic locus.
   TString sav = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/pid_plane_pt.root";
   TFile *fo = new TFile(sav, "RECREATE");
   Float_t v_sqrtdedx, v_brho, v_polar, v_ic;
   auto *pts = new TTree("pts", "PID points, one row per pattern track");
   pts->Branch("sqrtdedx", &v_sqrtdedx);
   pts->Branch("brho", &v_brho);
   pts->Branch("polar", &v_polar);
   pts->Branch("ic", &v_ic);

   auto *h = new TH2D("hpid", "a1954 ^{14}C PID plane;#sqrt{dE/dx};B#rho [Tm]", 300, 0, 30, 300, 0, 3.0);
   auto *hb = new TH1D("hb", "B#rho;B#rho [Tm];tracks", 120, 0, 3.0);
   auto *hp = new TH2D("hpol", "polar vs B#rho;#theta_{lab}^{reco} [deg];B#rho [Tm]", 180, 0, 180, 150, 0, 3.0);
   long long nt = 0, nv = 0;
   TObjArray *a = runsCSV.Tokenize(",");
   for (auto o : *a) {
      TString r = ((TObjString *)o)->String();
      // ---- the ion chamber, per ENTRY INDEX -------------------------------------------------
      // The IC trace lives in <run>_FRIB.root and the pattern tracks in <run>_slim.root; the two
      // are matched by entry index, exactly as pipeline/gate_events_C14.C matches them. Without
      // this the PID plane superimposes every beam species in the cocktail and a gate drawn on it
      // is a gate drawn on more than one beam at once.
      std::vector<float> icOf;
      {
         TFile *fF = TFile::Open(inDir + r + "_FRIB.root");
         TTree *tF = fF && !fF->IsZombie() ? (TTree *)fF->Get("cbmsim") : nullptr;
         if (!tF) {
            printf("\033[1;33m  %s: no FRIB file -- ic left at -1, the plane will be UNGATED\033[0m\n",
                   r.Data());
         } else {
            TClonesArray *ra = nullptr;
            tF->SetBranchAddress("AtRawEvent", &ra);
            icOf.assign(tF->GetEntries(), -1.f);
            long nPile = 0;
            for (Long64_t i = 0; i < tF->GetEntries(); ++i) {
               tF->GetEntry(i);
               if (!ra || ra->GetEntries() == 0) continue;
               auto *raw = (AtRawEvent *)ra->At(0);
               if (!raw || raw->GetGenTraces().empty()) continue;
               const auto &adc = raw->GetGenTraces()[0]->GetADC();
               double mx = -1e9;
               for (int b = icTbLo; b < icTbHi && b < (int)adc.size(); ++b) mx = std::max(mx, (double)adc[b]);
               // exactly one pulse, or the event is pile-up and is marked unusable rather than
               // being allowed through with a meaningless amplitude
               if (ppCountPulses(adc, peakThr, pkTbLo, pkTbHi) != 1) { ++nPile; continue; }
               icOf[i] = (float)mx;
            }
            printf("  %s: IC read for %lld entries, %ld rejected as pile-up\n",
                   r.Data(), (Long64_t)icOf.size(), nPile);
         }
         if (fF) fF->Close();
      }

      TFile *f = TFile::Open(inDir + r + "_slim.root");
      TTree *t = f && !f->IsZombie() ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t) { printf("  skip %s\n", r.Data()); if (f) f->Close(); continue; }
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPatternEvent", &pe);
      Long64_t N = maxEvt > 0 ? std::min(maxEvt, t->GetEntries()) : t->GetEntries();
      for (Long64_t i = 0; i < N; ++i) {
         t->GetEntry(i);
         if (!pe || pe->GetEntries() == 0) continue;
         auto *p = (AtPatternEvent *)pe->At(0);
         if (!p) continue;
         for (auto &trk : p->GetTrackCand()) {
            AtTrack &tr = const_cast<AtTrack &>(trk);
            ++nt;
            auto s = spy.Estimate(tr);
            if (!s.valid) continue;
            ++nv;
            h->Fill(s.sqrtdEdx, s.brho);
            hb->Fill(s.brho);
            hp->Fill(s.polar * TMath::RadToDeg(), s.brho);
            v_sqrtdedx = s.sqrtdEdx; v_brho = s.brho;
            v_polar = s.polar * TMath::RadToDeg();
            v_ic = (i < (Long64_t)icOf.size()) ? icOf[i] : -1.f;
            pts->Fill();
         }
      }
      printf("  %s: %lld entries\n", r.Data(), N);
      f->Close();
   }
   printf("\n  %lld tracks, %lld with a valid PID estimate\n\n", nt, nv);
   printf("  Brho band populations:\n");
   for (int k = 0; k < 12; ++k) {
      double lo = k * 0.25, hi = lo + 0.25;
      double n = hb->Integral(hb->FindBin(lo + 1e-6), hb->FindBin(hi - 1e-6));
      printf("    %.2f-%.2f  %8.0f  %s\n", lo, hi, n,
             std::string((int)std::min(60.0, n / std::max(1.0, hb->GetMaximum() / 60)), '#').c_str());
   }

   auto *c = new TCanvas("cpid", "", 1500, 500); c->Divide(3, 1);
   c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.13); h->Draw("colz");
   auto *box = new TBox(1.673, 0.040, 16.0, 0.964);        // the proton polygon's bounding box
   box->SetFillStyle(0); box->SetLineColor(kRed + 1); box->SetLineWidth(3); box->Draw("l");
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040); tx.SetTextColor(kRed + 1);
   tx.DrawLatex(0.15, 0.32, "proton gate");
   tx.SetTextColor(kWhite);
   tx.DrawLatex(0.15, 0.86, "triton expected here:");
   tx.DrawLatex(0.15, 0.80, "B#rho #approx 1.5 Tm");
   c->cd(2); gPad->SetLogy(); gPad->SetGridx(); hb->SetLineWidth(2); hb->Draw("hist");
   auto *l1 = new TLine(0.964, 0.5, 0.964, hb->GetMaximum());
   l1->SetLineColor(kRed + 1); l1->SetLineWidth(2); l1->SetLineStyle(2); l1->Draw();
   auto *l2 = new TLine(1.5, 0.5, 1.5, hb->GetMaximum());
   l2->SetLineColor(kGreen + 2); l2->SetLineWidth(3); l2->Draw();
   TLatex t2; t2.SetNDC(); t2.SetTextSize(0.042);
   t2.SetTextColor(kRed + 1); t2.DrawLatex(0.38, 0.60, "proton gate ends");
   t2.SetTextColor(kGreen + 2); t2.DrawLatex(0.55, 0.50, "triton expected");
   c->cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.13); hp->Draw("colz");
   // SAVE THE HISTOGRAMS. Reading the point clouds and running AtSpyralPID over them takes ~10
   // minutes a run; every question asked afterwards (dE/dx in a Brho slice, the polar distribution
   // of a candidate group) is a projection of what is already here, and should not cost another
   // pass over the data.
   { fo->cd(); pts->Write(); h->Write(); hb->Write(); hp->Write();
     TNamed prov("provenance", Form("runs=%s bField=%+.2f inDir=%s", runsCSV.Data(), bField, inDir.Data()));
     prov.Write(); printf("  saved %lld points + histograms to %s\n", pts->GetEntries(), sav.Data());
     fo->Close(); }
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c->SaveAs(out + "16_pid_plane_pt.png");
   printf("\n  wrote %s16_pid_plane_pt.png\n\n", out.Data());
}
