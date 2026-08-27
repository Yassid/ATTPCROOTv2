/// @file gate_draw_pt_C14.C
/// @brief Buttoned gate drawer for the 14C(p,t)12C triton gate, on the a1954 Spyral PID plane.
///        Adapted from a1975 D2_UKF/pid/gate_draw_dt.C -- same GUI, 14C kinematics.
///
/// Reads the ungated points cache from pid_plane_pt_C14.C (pattern tracks, NO FIT NEEDED), so the
/// gate can be drawn before anything is fitted. That matters here: the fit needs the gate, and the
/// (p,p') gated inputs are proton-only -- every triton was discarded at that step.
///
/// WHAT MAKES THIS DIFFERENT FROM A BARE TCutG, and why it is the tool to use:
///   * the 14C(p,t)12C g.s. Brho(theta) locus is drawn ON the PID plane, so the band being
///     selected has an identity instead of being a guess;
///   * [Locus check] shows Brho against polar angle for the tracks INSIDE the polygon, with the
///     same theoretical curve overlaid. A gate in (sqrt(dE/dx), Brho) is only trustworthy if what
///     it selects lands on a kinematic locus -- pick the wrong band and the cloud leaves its curve;
///   * the existing proton and deuteron gates are overlaid for orientation, so it is obvious which
///     band is which before a single vertex is placed.
///
/// [Draw new gate] -> click vertices, DOUBLE-CLICK to close -> saves IMMEDIATELY, then reports what
/// it keeps.  [Evaluate] re-reports.  [Locus check] is the one that matters.
///
/// RUN INTERACTIVELY -- it needs a display, so `root -l`, never `-b`.
///   cd .../a1954_Be12/UKF
///   root -l 'pid/gate_draw_pt_C14.C'
///
/// BEAM ENERGY: the locus is drawn at Ebeam = 159.75 MeV. Unlike the 12Be case this value is
/// anchored, not guessed -- it comes from scanning the beam energy until the resolved 6.094 MeV
/// level sits at its known excitation energy, and it is checked out-of-sample on the isolated
/// 8.317 MeV 2+. Expect the locus to be good to well under a percent in Brho, so a band that
/// misses it badly is NOT the triton.
///
/// WHAT TO EXPECT, AND THE REASON THIS MAY FIND NOTHING. 14C(p,t)12C has Q = -4.64 MeV, so at
/// Ebeam = 159.75 the triton emerges at theta_lab 8-24 deg with 14-56 MeV and Brho ~ 1.5 Tm --
/// far above the proton gate, which ends at 0.964. At that rigidity a triton carries 7.5x a
/// proton's dE/dx, so the two bands should be widely separated: look for a locus HIGH and to the
/// RIGHT, at reconstructed polar 156-172 deg (the reconstruction reports 180 - theta).
/// But a Brho of 1.5 Tm is a helix radius of ~52 cm in a chamber of ~29 cm, so these tracks LEAVE
/// rather than curl, and their rigidity comes from a partial arc. If the high-Brho population
/// turns out to be a smooth tail with no dE/dx structure, it is mis-fitted tracks and not tritons
/// -- which is the outcome the sibling 16C(p,t) analysis reached.
///
#include <vector>

/// Brho(theta_lab) for a two-body channel at Ebeam, so the tracks that MUST be selected can be
/// marked on the PID plane itself. Without this the gate is drawn on a band whose identity is a
/// guess -- which is exactly how a low-rigidity polygon got fitted for a whole night.
/// Brho(theta_lab) for 14C(p,t)12C at Ebeam, so the tracks that MUST be selected can be marked
/// on the PID plane itself. exStar lets a locus be drawn for an excited 12C state as well.
static double ptBrhoAtC14(double thlab, double Eb = 159.75, double exStar = 0.0)
{
   const double u = 931.49401;
   const double m1 = 14.003242 * u;   // 14C beam
   const double m2 = 1.007825 * u;    // proton target
   const double m3 = 3.016049 * u;    // triton ejectile
   const double m4 = 12.000000 * u + exStar;  // 12C residual
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

/// The locus as a curve, in whichever polar convention is asked for. flip = true means the stored
/// polar angle is 180 - theta_lab, which is what the a1975 drawer assumed; it is an ARGUMENT here
/// because the convention is exactly the kind of thing that is wrong silently.
static TGraph *ptLocusGraph(bool flip, double Eb, double exStar, int colour)
{
   auto *g = new TGraph();
   for (double th = 1; th <= 179; th += 1.0) {
      double b = ptBrhoAtC14(th, Eb, exStar);
      if (b <= 0) continue;
      g->SetPoint(g->GetN(), flip ? 180.0 - th : th, b);
   }
   g->SetLineColor(colour);
   g->SetLineWidth(3);
   return g;
}

class PtGateDraw : public TObject {
public:
   PtGateDraw(TString outJson, TString cache, TString refP, TString refD, double xMax, double yMax,
              double icLo, double icHi, bool showLocus, double eBeam, bool flipPolar)
      : fOut(outJson), fXmax(xMax), fYmax(yMax), fIcLo(icLo), fIcHi(icHi), fShowLocus(showLocus),
        fEbeam(eBeam), fFlip(flipPolar)
   {
      fH = new TH2F("hpt",
                    Form("a1954 ^{12}Be Spyral PID, IC gated [%.0f,%.0f]  --  draw the TRITON gate"
                         ";#sqrt{dE/dx};B#rho [T#upointm]", icLo, icHi),
                    300, 0, xMax, 300, 0, yMax);
      Load(cache);
      if (fOnLocus && fShowLocus)
         printf("=== %d tracks lie on the (p,t) g.s. locus (+-%.0f%%), drawn as red dots ===\n",
                fOnLocus->GetN(), 100 * fLocusTol);
      fRefP = LoadRef(refP, kRed + 1, "proton gate");
      fRefD = LoadRef(refD, kGreen + 2, "deuteron gate");
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
      // THE POINT OF THIS PANEL: the selected cloud must lie ON a kinematic curve. Both polar
      // conventions are drawn because which one the data uses is exactly what silently goes wrong
      // -- solid is the convention in use (flipPolar), dashed is the other one. If the selection
      // follows the DASHED curve, the flipPolar argument is set the wrong way round.
      auto *gUse = ptLocusGraph(fFlip, fEbeam, 0.0, kRed + 1);
      auto *gAlt = ptLocusGraph(!fFlip, fEbeam, 0.0, kGray + 2);
      gAlt->SetLineStyle(2);
      gAlt->Draw("L same");
      gUse->Draw("L same");
      // 12C 3.368 MeV 2+ : the first excited state, so the g.s. band can be told from its neighbour
      auto *gEx = ptLocusGraph(fFlip, fEbeam, 3.368, kMagenta + 1);
      gEx->SetLineStyle(7);
      gEx->Draw("L same");
      auto *leg = new TLegend(0.60, 0.72, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(gUse, Form("(p,t) g.s., E_{b}=%.0f", fEbeam), "l");
      leg->AddEntry(gEx, "(p,t) 12C 3.368 (2^{+})", "l");
      leg->AddEntry(gAlt, "other polar convention", "l");
      leg->Draw();
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
      fCanvas->SaveAs("pid/plots/gate_draw_pt_C14.png");
      printf("saved pid/plots/gate_draw_pt_C14.png\n");
   }

   /// Toggle the red locus markers. They are OFF by default: they are an orientation aid, and
   /// while placing vertices they sit exactly where the cursor needs to be and obscure the band
   /// being enclosed. The kinematic CURVE is a separate object and stays drawn either way, so
   /// turning these off does not remove the band's identity.
   /// Refill the displayed plane from the points already in memory, keeping only tracks inside
   /// the laboratory-angle window. This is a DISPLAY filter and a diagnostic: the clutter on the
   /// plane is dominated by tracks the physics analysis never uses -- the recoil proton of
   /// 14C(p,p') emerges at theta_lab < 90 deg, and 14C(p,t) tritons at 8-24 deg, so anything
   /// outside that is either a beam-like track, a fragment, or a piece of a track that
   /// pattern recognition split. Narrowing the window shows how much of the plane is which.
   void ApplyAngle()
   {
      if (fEThLo) fThLo = fEThLo->GetNumber();
      if (fEThHi) fThHi = fEThHi->GetNumber();
      fH->Reset();
      long kept = 0;
      for (size_t i = 0; i < fX.size(); ++i) {
         double thLab = fFlip ? 180.0 - fPol[i] : fPol[i];
         if (thLab < fThLo || thLab > fThHi) continue;
         fH->Fill(fX[i], fY[i]);
         ++kept;
      }
      printf("theta_lab in [%.1f, %.1f] deg : %ld of %zu tracks kept (%.1f%%)\n",
             fThLo, fThHi, kept, fX.size(), 100.0 * kept / std::max<size_t>(1, fX.size()));
      if (fLabel)
         fLabel->SetText(Form("  theta_lab %.0f-%.0f deg : %ld of %zu tracks (%.0f%%)",
                              fThLo, fThHi, kept, fX.size(), 100.0 * kept / std::max<size_t>(1, fX.size())));
      Redraw();
   }

   void ToggleLocus()
   {
      fShowLocus = !fShowLocus;
      printf("locus markers %s\n", fShowLocus ? "ON" : "OFF");
      Redraw();
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
      if (fRefP) fRefP->Draw("L");
      if (fRefD) fRefD->Draw("L");
      if (fNew) fNew->Draw("L");
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
      // NO ic BRANCH MEANS NO IC CUT, not an IC of -1. Left as it was, a points file without the
      // branch put every track at ic = -1, outside any window, and the drawer opened on an empty
      // plane reporting "0 inside IC (0.0%)" -- which reads as "this beam has no such particles"
      // rather than as a missing column.
      const bool hasIC = t->GetBranch("ic") != nullptr;
      if (hasIC)
         t->SetBranchAddress("ic", &ic);
      else
         printf("\033[1;33mno ic branch in the points file -- IC gate NOT applied, the plane is "
                "UNGATED and carries every beam species\033[0m\n");
      long nAll = 0, nCut = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         ++nAll;
         // The D2 beam is a cocktail: 16C near IC 1150 and a second component near 2060 in
         // comparable numbers. A gate drawn on the mixture is drawn on two beams at once.
         if (hasIC && (ic < fIcLo || ic > fIcHi))
            continue;
         ++nCut;
         fX.push_back(x);
         fY.push_back(y);
         fPol.push_back(pol);
         fH->Fill(x, y);
         // mark the tracks sitting on the (p,t) g.s. locus. Whether the stored polar is theta_lab
         // or 180 - theta_lab is a CONVENTION, so it is an argument rather than an assumption.
         double b = ptBrhoAtC14(fFlip ? 180.0 - pol : pol, fEbeam);
         if (b > 0 && std::fabs(y - b) / b < fLocusTol) {
            if (!fOnLocus) fOnLocus = new TGraph();
            fOnLocus->SetPoint(fOnLocus->GetN(), x, y);
         }
      }
      printf("=== %s: %ld tracks, %ld inside IC [%.0f,%.0f] (%.1f%%) ===\n", cache.Data(), nAll, nCut,
             fIcLo, fIcHi, nAll ? 100.0 * nCut / nAll : 0.0);
      f->Close();
   }

   /// Overlay a reference polygon. Shown, never enforced -- these are the EXISTING proton and
   /// deuteron gates, drawn so it is obvious which band is which before the first vertex.
   TCutG *LoadRef(TString ref, int colour, const char *what)
   {
      if (ref.IsNull()) return nullptr;
      std::ifstream in(ref.Data());
      if (!in) { printf("no %s at %s (skipped)\n", what, ref.Data()); return nullptr; }
      std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      auto p = s.find('[', s.find("vertices"));
      if (p == std::string::npos) return nullptr;
      std::vector<double> v;
      const char *q = s.c_str() + p, *e = s.c_str() + s.size();
      int depth = 0;
      while (q < e) {
         if (*q == '[') depth++;
         if (*q == ']') { depth--; if (depth <= 0) break; }
         char *np = nullptr;
         double d = strtod(q, &np);
         if (np != q) { v.push_back(d); q = np; } else ++q;
      }
      if (v.size() < 6) return nullptr;
      auto *g = new TCutG(Form("ref_%s", what), v.size() / 2);
      for (size_t i = 0, k = 0; i + 1 < v.size(); i += 2, ++k) g->SetPoint(k, v[i], v[i + 1]);
      g->SetLineColor(colour); g->SetLineStyle(2); g->SetLineWidth(2);
      printf("%s overlaid (dashed): %s\n", what, ref.Data());
      return g;
   }

   void MakeGui()
   {
      auto *main = new TGMainFrame(gClient->GetRoot(), 1150, 900);
      main->SetWindowName("a1954 14C(p,t)12C PID  --  draw the triton gate");
      auto *bar = new TGHorizontalFrame(main);
      const char *lbl[] = {"Draw new gate", "Evaluate", "Locus check", "Locus dots on/off",
                           "Save JSON", "Save PNG", "Redraw", "Quit"};
      const char *slot[] = {"DrawGate()", "Evaluate()", "LocusCheck()", "ToggleLocus()",
                            "Save()", "SavePNG()", "Redraw()", "Quit()"};
      for (int i = 0; i < 8; ++i) {
         auto *b = new TGTextButton(bar, lbl[i]);
         b->Connect("Clicked()", "PtGateDraw", this, slot[i]);
         bar->AddFrame(b, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      }
      bar->AddFrame(new TGLabel(bar, "name:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      fName = new TGTextEntry(bar, "triton_14C");
      fName->Resize(180, 22);
      bar->AddFrame(fName, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      // laboratory-angle window. Shown in the TRUE lab convention (180 - stored polar when
      // flipPolar is set), because that is the frame the kinematics are quoted in: the (p,p')
      // recoil proton lives below 90 deg and the (p,t) triton between 8 and 24.
      bar->AddFrame(new TGLabel(bar, "  theta_lab:"),
                    new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      fEThLo = new TGNumberEntry(bar, 0.0, 5, -1, TGNumberFormat::kNESRealOne,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0.0, 180.0);
      fEThLo->Resize(70, 22);
      fEThLo->Connect("ValueSet(Long_t)", "PtGateDraw", this, "ApplyAngle()");
      bar->AddFrame(fEThLo, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
      bar->AddFrame(new TGLabel(bar, "to"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 4, 4));
      fEThHi = new TGNumberEntry(bar, 180.0, 5, -1, TGNumberFormat::kNESRealOne,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0.0, 180.0);
      fEThHi->Resize(70, 22);
      fEThHi->Connect("ValueSet(Long_t)", "PtGateDraw", this, "ApplyAngle()");
      bar->AddFrame(fEThHi, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      {
         auto *ba = new TGTextButton(bar, "Apply angle");
         ba->Connect("Clicked()", "PtGateDraw", this, "ApplyAngle()");
         bar->AddFrame(ba, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      }
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      fLabel = new TGLabel(main, "  Draw around the band high and right of the proton gate, then Locus check. "
                                 "[Locus dots on/off] marks tracks on the (p,t) g.s. curve.");
      main->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 6, 6, 2, 2));
      auto *ec = new TRootEmbeddedCanvas("ec", main, 1120, 800);
      main->AddFrame(ec, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
      fCanvas = ec->GetCanvas();
      main->MapSubwindows();
      main->Resize(main->GetDefaultSize());
      // PLACE THE WINDOW EXPLICITLY. Launched detached under WSLg the frame was mapped at
      // +3790+285 -- off the visible desktop, so the GUI was running and invisible, which looks
      // exactly like a crash. Pinning it near the origin costs nothing and removes that failure.
      main->Move(40, 40);
      main->MapWindow();
      main->RaiseWindow();
   }

   TString fOut;
   double fXmax, fYmax;
   TH2F *fH = nullptr;
   TCutG *fNew = nullptr, *fRefP = nullptr, *fRefD = nullptr;
   TCanvas *fCanvas = nullptr;
   TGNumberEntry *fEThLo = nullptr, *fEThHi = nullptr;
   double fThLo = 0.0, fThHi = 180.0;
   TGLabel *fLabel = nullptr;
   TGTextEntry *fName = nullptr;
   std::vector<float> fX, fY, fPol;
   TGraph *fOnLocus = nullptr;
   double fLocusTol = 0.15;
   bool fShowLocus = false; // off by default -- the dots are a hint, not part of the data
   double fIcLo = -1e9, fIcHi = 1e9;
   double fEbeam = 159.75;
   bool fFlip = true;
   TGraph *fLocusPlane = nullptr;
};

/// @param outJson   where the gate is written (autosaved the moment the polygon closes)
/// @param cache     pts tree from pid_plane_pt_C14.C
/// @param refP/refD existing proton / deuteron gates, overlaid dashed for orientation
/// @param icLo/icHi 14C beam window on the FRIB IC amplitude. 625-750 is the 14C peak; the
///                  dominant ~1900 component is a CONTAMINANT, not the beam.
/// @param eBeam     only the locus depends on it, not the gate. 155 is the VOID July value.
/// @param flipPolar stored polar = 180 - theta_lab. Check it in the Locus check panel.
void gate_draw_pt_C14(TString outJson = "pid/triton_14C.json",
                       TString cache = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/pid_plane_pt.root",
                       TString refP = "pid/proton_14C.json", TString refD = "",
                       double xMax = 30.0, double yMax = 3.0, double icLo = 950, double icHi = 1350,
                       bool showLocus = false, double eBeam = 159.75, bool flipPolar = true)
{
   new PtGateDraw(outJson, cache, refP, refD, xMax, yMax, icLo, icHi, showLocus, eBeam, flipPolar);
}
