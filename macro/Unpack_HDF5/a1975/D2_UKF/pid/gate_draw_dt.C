/// @file gate_draw_dt.C
/// @brief Buttoned gate drawer for the a1975 (d,t) Spyral PID plane at dv 1.10424.
///        Sibling of 16C_pp_a1975/gate_draw_sim.C, same GUI, different input.
///
/// Reads the ungated points cache from pid/pid_plane_dt.C -- pattern tracks, NO FIT NEEDED, so
/// the gate can be drawn before anything is fitted. That matters here: the fit needs the gate.
///
/// The polygon is drawn with ROOT's own CUTG tool via WaitPrimitive, not a hand-written click
/// handler. A hand-rolled handler picks up the spurious button event that arrives when the
/// canvas is first mapped and silently plants a phantom vertex.
///
/// [Draw new gate] -> click vertices, DOUBLE-CLICK to close -> saves immediately, then reports
/// what it keeps.  [Evaluate] re-reports.  [Save JSON] rewrites with the current name.
/// [Locus check] is the one that matters: it shows Brho against polar angle for the tracks
/// INSIDE the polygon. A gate in (sqrt(dE/dx), Brho) is only trustworthy if what it selects
/// lands on a kinematic locus -- pick the wrong band and the selected cloud leaves its curve.
///
///   cd .../a1975/D2_UKF
///   root -l 'pid/gate_draw_dt.C'
///   root -l 'pid/gate_draw_dt.C("pid/triton_d2_dv1104.json")'

#include <vector>

/// Brho(theta_lab) for a two-body channel at Ebeam, so the tracks that MUST be selected can be
/// marked on the PID plane itself. Without this the gate is drawn on a band whose identity is a
/// guess -- which is exactly how a low-rigidity polygon got fitted for a whole night.
static double dtBrhoAt(double thlab, double Eb = 184.17)
{
   const double u = 931.49401, m1 = 16.0147013 * u, m2 = 2.0135532 * u;
   const double m3 = 3.01550072 * u, m4 = 15.0105993 * u;
   double E1 = Eb + m1, pb = std::sqrt(E1 * E1 - m1 * m1), thr = thlab * TMath::DegToRad(), best = -1;
   for (double ke = 0.2; ke < 170; ke += 0.05) {
      double E3 = ke + m3, p3 = std::sqrt(E3 * E3 - m3 * m3), E4 = E1 + m2 - E3;
      double px = p3 * std::sin(thr), pz = p3 * std::cos(thr);
      double m4x2 = E4 * E4 - (px * px + (pb - pz) * (pb - pz));
      if (m4x2 < 0) continue;
      if (std::fabs(std::sqrt(m4x2) - m4) < 0.06) best = ke;
   }
   if (best < 0) return -1;
   return std::sqrt(2 * m3 * best + best * best) / 299.792458;
}

class DtGateDraw : public TObject {
public:
   DtGateDraw(TString outJson, TString cache, TString refGate, double xMax, double yMax, double icLo,
              double icHi, bool showLocus)
      : fOut(outJson), fXmax(xMax), fYmax(yMax), fIcLo(icLo), fIcHi(icHi), fShowLocus(showLocus)
   {
      fH = new TH2F("hdt",
                    Form("^{16}C(d,t) Spyral PID, dv 1.10424, IC gated [%.0f,%.0f]"
                         ";#sqrt{dE/dx};B#rho [T#upointm]", icLo, icHi),
                    300, 0, xMax, 300, 0, yMax);
      Load(cache);
      if (fOnLocus && fShowLocus)
         printf("=== %d tracks lie on the (d,t) g.s. locus (+-%.0f%%), drawn as red dots ===\n",
                fOnLocus->GetN(), 100 * fLocusTol);
      LoadRef(refGate);
      MakeGui();
      Redraw();
   }

   void DrawGate()
   {
      fCanvas->cd();
      printf("\n>>> click polygon vertices, DOUBLE-CLICK to close\n");
      TCutG *c = (TCutG *)fCanvas->WaitPrimitive("CUTG", "CutG");
      if (!c) {
         printf("no cut drawn\n");
         return;
      }
      fNew = (TCutG *)c->Clone("newgate");
      fNew->SetLineColor(kGreen + 2);
      fNew->SetLineWidth(3);
      // The drawn cut belongs to the pad: unlink before deleting, or the next Redraw walks a
      // dangling pointer -- the "realloc(): invalid next size" that has eaten a finished gate.
      if (fCanvas->GetListOfPrimitives())
         fCanvas->GetListOfPrimitives()->Remove(c);
      delete c;
      Save(); // save FIRST -- a polygon that exists only in memory is one segfault from gone
      Evaluate();
      Redraw();
   }

   void Evaluate()
   {
      if (!fNew) {
         printf("draw a gate first\n");
         return;
      }
      long in = 0;
      double pmin = 1e9, pmax = -1e9;
      for (size_t i = 0; i < fX.size(); ++i)
         if (fNew->IsInside(fX[i], fY[i])) {
            ++in;
            pmin = std::min(pmin, (double)fPol[i]);
            pmax = std::max(pmax, (double)fPol[i]);
         }
      TString m = Form("  keeps %ld / %zu tracks = %.1f %%   polar %.0f-%.0f deg", in, fX.size(),
                       fX.empty() ? 0.0 : 100.0 * in / fX.size(), in ? pmin : 0., in ? pmax : 0.);
      fLabel->SetText(m);
      printf("%s\n", m.Data());
   }

   /// Brho vs polar angle for the SELECTED tracks: the check that says whether the band is the
   /// right one. Shown against the whole sample so the selection is visible in context.
   void LocusCheck()
   {
      if (!fNew) {
         printf("draw a gate first\n");
         return;
      }
      auto *c2 = new TCanvas("clocus", "Brho vs polar -- inside the gate", 1000, 780);
      // Binning: the selected sample is a few thousand tracks over a narrow Brho band, so the
      // 180x200 grid used before put ~0.2 counts in a typical bin and the locus read as noise.
      // 1 deg in polar and a Brho range clipped to what the GATE actually spans (padded 20 %)
      // puts the structure where it can be seen. colz, not markers -- markers hide density.
      double ylo = 1e9, yhi = -1e9;
      for (size_t i = 0; i < fX.size(); ++i)
         if (fNew->IsInside(fX[i], fY[i])) {
            ylo = std::min(ylo, (double)fY[i]);
            yhi = std::max(yhi, (double)fY[i]);
         }
      if (ylo > yhi) {
         printf("gate selects nothing\n");
         return;
      }
      double pad = 0.2 * (yhi - ylo) + 1e-3;
      ylo = std::max(0.0, ylo - pad);
      yhi = yhi + pad;
      auto *hin = new TH2F("hin",
                           Form("B#rho vs polar, INSIDE the gate  (%zu tracks)"
                                ";#theta_{polar} [deg];B#rho [T#upointm]",
                                fX.size()),
                           180, 0, 180, 120, ylo, yhi);
      auto *hall = new TH2F("hall", "", 180, 0, 180, 120, ylo, yhi);
      long nin = 0;
      for (size_t i = 0; i < fX.size(); ++i) {
         hall->Fill(fPol[i], fY[i]);
         if (fNew->IsInside(fX[i], fY[i])) {
            hin->Fill(fPol[i], fY[i]);
            ++nin;
         }
      }
      hin->SetTitle(Form("B#rho vs polar, INSIDE the gate  (%ld tracks)"
                         ";#theta_{polar} [deg];B#rho [T#upointm]",
                         nin));
      c2->Divide(1, 2);
      c2->cd(1);
      gPad->SetLogz();
      gPad->SetRightMargin(0.13);
      hin->Draw("colz");
      c2->cd(2);
      gPad->SetRightMargin(0.13);
      // the polar projection of the selection: where in angle the gate actually lives
      auto *px = hin->ProjectionX("px");
      px->SetTitle("polar angle of the selected tracks;#theta_{polar} [deg];tracks");
      px->SetLineColor(kRed + 1);
      px->SetLineWidth(2);
      px->Draw("hist");
      c2->Modified();
      c2->Update();
      printf("locus check drawn: grey = all tracks, red = inside the gate\n");
   }

   void Save()
   {
      if (!fNew) {
         printf("no gate drawn -- nothing to save\n");
         return;
      }
      FILE *f = fopen(fOut.Data(), "w");
      if (!f) {
         printf("cannot write %s\n", fOut.Data());
         return;
      }
      // Z and A are NOT decoration: without them the loader cannot tell one polygon from another
      // when several are open, and the wrong one gets applied.
      fprintf(f, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrtdedx\",\n    \"yaxis\": \"brho\",\n",
              fName->GetText());
      fprintf(f, "    \"Z\": 1,\n    \"A\": 3,\n    \"vertices\": [\n");
      int np = fNew->GetN();
      for (int i = 0; i < np; ++i) {
         double x, y;
         fNew->GetPoint(i, x, y);
         fprintf(f, "        [%.4f, %.4f]%s\n", x, y, (i == np - 1) ? "" : ",");
      }
      fprintf(f, "    ]\n}\n");
      fclose(f);
      printf("saved %d vertices -> %s\n", np, fOut.Data());
      fLabel->SetText(Form("  saved %d vertices -> %s", np, fOut.Data()));
   }

   void SavePNG()
   {
      gSystem->mkdir("pid/plots", kTRUE);
      fCanvas->SaveAs("pid/plots/gate_draw_dt.png");
      printf("saved pid/plots/gate_draw_dt.png\n");
   }

   void Redraw()
   {
      fCanvas->cd();
      fCanvas->SetLogz();
      fCanvas->SetRightMargin(0.13);
      fH->Draw("colz");
      if (fOnLocus && fShowLocus) { // opt-in: tracks already known to lie on the (d,t) curve
         fOnLocus->SetMarkerColor(kRed);
         fOnLocus->SetMarkerStyle(20);
         fOnLocus->SetMarkerSize(0.45);
         fOnLocus->Draw("P same");
      }
      if (fRef)
         fRef->Draw("L");
      if (fNew)
         fNew->Draw("L");
      fCanvas->Modified();
      fCanvas->Update();
      gSystem->ProcessEvents();
   }

   void Quit() { gApplication->Terminate(0); }

private:
   void Load(TString cache)
   {
      TFile *f = TFile::Open(cache);
      if (!f || f->IsZombie()) {
         printf("cannot open %s -- run pid/pid_plane_dt.C first\n", cache.Data());
         return;
      }
      auto *t = (TTree *)f->Get("pts");
      if (!t) {
         printf("no tree 'pts' in %s\n", cache.Data());
         return;
      }
      float x, y, pol, ic = -1;
      t->SetBranchAddress("sqrtdedx", &x);
      t->SetBranchAddress("brho", &y);
      t->SetBranchAddress("polar", &pol);
      if (t->GetBranch("ic"))
         t->SetBranchAddress("ic", &ic);
      long nAll = 0, nCut = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         ++nAll;
         // The D2 beam is a cocktail: 16C near IC 1150 and a second component near 2060 in
         // comparable numbers. A gate drawn on the mixture is drawn on two beams at once.
         if (ic < fIcLo || ic > fIcHi)
            continue;
         ++nCut;
         fX.push_back(x);
         fY.push_back(y);
         fPol.push_back(pol);
         fH->Fill(x, y);
         // mark the tracks sitting on the (d,t) g.s. locus: theta_lab = 180 - theta_polar
         double b = dtBrhoAt(180.0 - pol);
         if (b > 0 && std::fabs(y - b) / b < fLocusTol) {
            if (!fOnLocus) fOnLocus = new TGraph();
            fOnLocus->SetPoint(fOnLocus->GetN(), x, y);
         }
      }
      printf("=== %s: %ld tracks, %ld inside IC [%.0f,%.0f] (%.1f%%) ===\n", cache.Data(), nAll, nCut,
             fIcLo, fIcHi, nAll ? 100.0 * nCut / nAll : 0.0);
      f->Close();
   }

   /// Overlay a reference polygon. It is SHOWN, never enforced: the old triton_d2.json belongs
   /// to the dv 1.136 plane, and Brho moves with the drift velocity.
   void LoadRef(TString ref)
   {
      if (ref.IsNull())
         return;
      std::ifstream in(ref.Data());
      if (!in)
         return;
      std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      auto p = s.find('[', s.find("vertices"));
      if (p == std::string::npos)
         return;
      std::vector<double> v;
      const char *q = s.c_str() + p, *e = s.c_str() + s.size();
      int depth = 0;
      while (q < e) {
         if (*q == '[')
            depth++;
         if (*q == ']') {
            depth--;
            if (depth <= 0)
               break;
         }
         char *np = nullptr;
         double d = strtod(q, &np);
         if (np != q) {
            v.push_back(d);
            q = np;
         } else
            ++q;
      }
      if (v.size() < 6)
         return;
      fRef = new TCutG("refgate", v.size() / 2);
      for (size_t i = 0, k = 0; i + 1 < v.size(); i += 2, ++k)
         fRef->SetPoint(k, v[i], v[i + 1]);
      fRef->SetLineColor(kGray + 2);
      fRef->SetLineStyle(2);
      fRef->SetLineWidth(2);
      printf("reference overlaid (dashed grey, OLD dv -- for orientation only): %s\n", ref.Data());
   }

   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1150, 900);
      main->SetWindowName("a1975 (d,t) PID  --  draw the triton gate");
      auto *bar = new TGHorizontalFrame(main);
      const char *lbl[] = {"Draw new gate", "Evaluate", "Locus check", "Save JSON", "Save PNG", "Redraw", "Quit"};
      const char *slot[] = {"DrawGate()", "Evaluate()", "LocusCheck()", "Save()", "SavePNG()", "Redraw()", "Quit()"};
      for (int i = 0; i < 7; ++i) {
         auto *b = new TGTextButton(bar, lbl[i]);
         b->Connect("Clicked()", "DtGateDraw", this, slot[i]);
         bar->AddFrame(b, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      }
      bar->AddFrame(new TGLabel(bar, "name:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      fName = new TGTextEntry(bar, "triton_d2_dv1104");
      fName->Resize(180, 22);
      bar->AddFrame(fName, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      fLabel = new TGLabel(main, "  draw a gate around the triton band, then use Locus check");
      main->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 6, 6, 2, 2));
      auto *ec = new TRootEmbeddedCanvas("ec", main, 1120, 800);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();
      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      main->MapWindow();
   }

   TString fOut;
   double fXmax, fYmax;
   TH2F *fH = nullptr;
   TCutG *fNew = nullptr, *fRef = nullptr;
   TCanvas *fCanvas = nullptr;
   TGLabel *fLabel = nullptr;
   TGTextEntry *fName = nullptr;
   std::vector<float> fX, fY, fPol;
   TGraph *fOnLocus = nullptr;
   double fLocusTol = 0.15;
   bool fShowLocus = false; // off by default -- the dots are a hint, not part of the data
   double fIcLo = -1e9, fIcHi = 1e9;
};

void gate_draw_dt(TString outJson = "pid/triton_d2_dv1104.json",
                  TString cache = "pid/pid_plane_dt_dv1104.root", TString refGate = "pid/triton_d2.json",
                  double xMax = 55.0, double yMax = 4.0, double icLo = 900, double icHi = 1300,
                  bool showLocus = false)
{
   new DtGateDraw(outJson, cache, refGate, xMax, yMax, icLo, icHi, showLocus);
}
