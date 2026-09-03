/// @file gate_draw_C15p.C
/// @brief Buttoned gate drawer for the C15p PID plane. Ported from the a1954 14C(p,t) drawer
///        (a1954/UKF/pid/gate_draw_pt_C14.C) -- same GUI, this analysis's data and kinematics.
///
/// Reads pid/points_C15p.root (built by pid/make_points_C15p.C from the pattern-recognition PID
/// caches), so the gate can be drawn before anything is fitted -- which matters, because the fit
/// needs the gate.
///
/// WHAT THIS GIVES OVER A BARE TCutG:
///   * the two-body Brho(theta) locus can be drawn ON the plane, so a band has an identity rather
///     than being a guess;
///   * [Locus check] plots Brho against polar angle for the tracks INSIDE the polygon with the
///     same theoretical curve overlaid. A gate in (sqrt(dE/dx), Brho) is only trustworthy if what
///     it selects lands on a kinematic locus;
///   * previously drawn gates can be overlaid for orientation before a vertex is placed;
///   * axis ranges, binning and the IC window are adjustable in-GUI.
///
/// [Draw new gate] -> click vertices, DOUBLE-CLICK to close -> saves IMMEDIATELY, then reports what
/// it keeps.  [Evaluate] re-reports.  [Locus check] is the one that matters.
///
/// RUN INTERACTIVELY -- it needs a display, so `root -l`, never `-b`:
///   cd .../macro/Unpack_HDF5/C15p
///   root -l 'pid/gate_draw_C15p.C'
///
/// ★ THE POINTS FILE IS ALREADY GAIN MATCHED (make_points_C15p.C applies the per-run table). Do
/// not apply it again, and do not mix a gate drawn here with a raw plane: the per-run drift over
/// this run set is ~29 %.
///
/// ★ IC IS OFF BY DEFAULT (icLo = icHi = -1) because the beam window for this experiment has not
/// been chosen yet -- the IC spectrum has to be looked at first. Until then the plane carries
/// EVERY beam species, and a band drawn on it may be a mixture. Pass a window once one is known.
///
/// ★ EBEAM DEFAULTS TO 0, WHICH DISABLES THE LOCUS, and that is deliberate rather than lazy: the
/// 15C+p beam energy is not established here, and a locus drawn at a guessed energy is worse than
/// no locus at all -- it invites a gate to be drawn onto a curve that is in the wrong place. Pass
/// a real eBeam (MeV, lab) once it is calibrated and the locus and [Locus check] become available.
///
/// Default channel is 15C(d,p)16C -- the proton ejectile. m3Amu/m4Amu select another:
///   (d,d) elastic : m3Amu 2.014102,  m4Amu 15.0105993
///   (d,t) 14C     : m3Amu 3.0160493, m4Amu 14.003242
#include <vector>

/// Brho(theta_lab) for a two-body channel at Ebeam, so the tracks that MUST be selected can be
/// marked on the PID plane itself. Without this the gate is drawn on a band whose identity is a
/// guess -- which is exactly how a low-rigidity polygon got fitted for a whole night.
/// Brho(theta_lab) for 14C(p,t)12C at Ebeam, so the tracks that MUST be selected can be marked
/// on the PID plane itself. exStar lets a locus be drawn for an excited 12C state as well.
/// Brho(theta_lab) for any two-body channel on the 14C beam. m3Amu/m4Amu default to the (p,t)
/// triton and 12C; passing the deuteron and 13C draws the 14C(p,d)13C band instead, which is the
/// contaminant sitting next to the tritons on this plane.
static double dpBrhoC15p(double thlab, double Eb = 159.75, double exStar = 0.0,
                          double m3Amu = 1.007825, double m4Amu = 15.0105993)
{
   const double u = 931.49401;
   const double m1 = 15.0105993 * u;  // 15C beam
   const double m2 = 1.00782503 * u;  // PROTON target (ATOMIC masses throughout: the electrons
                                      // cancel side to side, and mixing atomic with nuclear here
                                      // is a ~0.5 MeV error that looks like a calibration offset)
   const double m3 = m3Amu * u;       // ejectile
   const double m4 = m4Amu * u + exStar;      // residual
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
static TGraph *dpLocusGraph(bool flip, double Eb, double exStar, int colour,
                            double m3Amu = 1.007825, double m4Amu = 15.0105993)
{
   // Eb <= 0 means the beam energy is not established, and a locus drawn at a guessed energy is
   // worse than none: it invites a gate onto a curve that is in the wrong place. Return nothing.
   if (!(Eb > 0))
      return nullptr;
   auto *g = new TGraph();
   for (double th = 1; th <= 179; th += 1.0) {
      double b = dpBrhoC15p(th, Eb, exStar, m3Amu, m4Amu);
      if (b <= 0) continue;
      g->SetPoint(g->GetN(), flip ? 180.0 - th : th, b);
   }
   g->SetLineColor(colour);
   g->SetLineWidth(3);
   return g;
}


/// Ejectile and residual masses for the species the drawer is gating, so [Locus check] compares a
/// band against ITS OWN channel. Everything below hangs off one fact: the compound is 15C + p, so
/// Z = 6+1 = 7 and A = 15+2 = 17, and fixing the ejectile fixes the residual.
///     ejectile p (1,1) -> 16C      d (1,2) -> 15C      t (1,3) -> 14C
///              3He (2,3) -> 14B    alpha (2,4) -> 13B
/// Before this, every locus used the (d,p) proton/16C default whatever species was being drawn --
/// so checking a deuteron gate against "the locus" silently checked it against the proton channel.
static void c15pChannelMasses(int Z, int A, double &m3Amu, double &m4Amu)
{
   struct E { int Z, A; double m3, m4; };
   static const E tab[] = {
      {1, 1, 1.00782503, 15.0105993},   // (p,p')  -> 15C   elastic / inelastic
      {1, 2, 2.01410178, 14.0032420},   // (p,d)   -> 14C
      {1, 3, 3.01604928, 13.0033548},   // (p,t)   -> 13C
      {2, 3, 3.01602932, 13.0177800},   // (p,3He) -> 13B
      {2, 4, 4.00260325, 12.0143521},   // (p,a)   -> 12B
   };
   for (const auto &e : tab)
      if (e.Z == Z && e.A == A) { m3Amu = e.m3; m4Amu = e.m4; return; }
   m3Amu = A;
   m4Amu = 16 - A;
}

class C15pGateDraw : public TObject {
public:
   C15pGateDraw(TString outJson, TString cache, TString refP, TString refD, double xMax, double yMax,
              double icLo, double icHi, bool showLocus, double eBeam, bool flipPolar, int Z, int A)
      : fOut(outJson), fXmax(xMax), fYmax(yMax), fIcLo(icLo), fIcHi(icHi), fShowLocus(showLocus),
        fEbeam(eBeam), fFlip(flipPolar), fZ(Z), fA(A)
   {
      c15pChannelMasses(fZ, fA, fM3, fM4);
      printf("channel: ejectile Z=%d A=%d (m3 %.6f u), residual m4 %.6f u\n", fZ, fA, fM3, fM4);

      fH = new TH2F("hpt",
                    Form("C15p PID, gain matched, %s  --  draw the gate"
                         ";#sqrt{dE/dx};B#rho [T#upointm]",
                         (icLo >= 0 && icHi > icLo) ? Form("IC [%.0f,%.0f]", icLo, icHi) : "IC gate OFF"),
                    300, 0, xMax, 300, 0, yMax);
      fXhi = xMax; fYhi = yMax;
      Load(cache);
      if (fOnLocus && fShowLocus)
         printf("=== %d tracks lie on the (d,p) g.s. locus (+-%.0f%%), drawn as red dots ===\n",
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
      auto *gUse = dpLocusGraph(fFlip, fEbeam, 0.0, kRed + 1, fM3, fM4);
      auto *gAlt = dpLocusGraph(!fFlip, fEbeam, 0.0, kGray + 2, fM3, fM4);
      gAlt->SetLineStyle(2);
      gAlt->Draw("L same");
      gUse->Draw("L same");
      // 12C 3.368 MeV 2+ : the first excited state, so the g.s. band can be told from its neighbour
      auto *gEx = dpLocusGraph(fFlip, fEbeam, 3.368, kMagenta + 1, fM3, fM4);
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
      // ★ THE FILENAME FOLLOWS THE NAME FIELD, so there is ONE source of truth. Previously the
      // path came from the launcher argument while the name field was independent, so typing
      // "proton_C15p" while the drawer had been opened for the deuteron wrote a gate called
      // proton_C15p into deuteron_C15p.json -- with the deuteron's Z and A. Both files then
      // disagreed with their own contents, and only a band-position measurement could tell which
      // was which.
      TString outPath = fOut;
      TString typed = fName->GetText();
      // ...and the SPECIES follows it too. Deriving the filename from the name while leaving Z/A
      // on the launcher's value still lets proton_C15p.json be written with the deuteron's A=2,
      // which is exactly the half-fix that hides the problem instead of removing it.
      {
         TString low = typed;
         low.ToLower();
         struct { const char *key; int Z, A; } kSpecies[] = {
            {"proton", 1, 1}, {"deuteron", 1, 2}, {"triton", 1, 3}, {"3he", 2, 3}, {"alpha", 2, 4},
            {"4he", 2, 4}};
         for (auto &sp : kSpecies)
            if (low.Contains(sp.key)) {
               if (fZ != sp.Z || fA != sp.A) {
                  printf("\033[1;33mnote: name '%s' implies Z=%d A=%d; overriding the Z=%d A=%d this "
                         "drawer was opened with.\033[0m\n", typed.Data(), sp.Z, sp.A, fZ, fA);
                  fZ = sp.Z;
                  fA = sp.A;
               }
               break;
            }
      }
      typed = typed.Strip(TString::kBoth);
      if (typed.Length()) {
         TString dirPart = gSystem->DirName(fOut.Data());
         outPath = (dirPart.Length() && dirPart != ".") ? dirPart + "/" + typed + ".json" : typed + ".json";
         if (outPath != fOut)
            printf("\033[1;33mnote: writing to %s (from the name field), not %s\033[0m\n",
                   outPath.Data(), fOut.Data());
      }
      // ★ NEVER SILENTLY REPLACE A DIFFERENT GATE. Following the name field means a Save can now
      // land on a file the drawer was not opened for -- which is how a triton polygon, saved once
      // under the name "proton_C15p" while the name was being corrected, overwrote the real proton
      // gate. The old polygon is kept as .bak and the replacement is announced.
      if (!gSystem->AccessPathName(outPath.Data())) {
         TString bak = outPath + ".bak";
         gSystem->CopyFile(outPath.Data(), bak.Data(), kTRUE);
         printf("\033[1;31mNOTE: %s already existed and has been REPLACED. Previous version kept as "
                "%s -- check you did not mean a different file.\033[0m\n", outPath.Data(), bak.Data());
      }
      FILE *f = fopen(outPath.Data(), "w");
      if (!f) {
         printf("cannot write %s\n", outPath.Data());
         return;
      }
      // Z and A are NOT decoration: without them the loader cannot tell one polygon from another
      // when several are open, and the wrong one gets applied.
      // ★ Z, A AND THE AXIS NAME WERE HARDCODED TO THE a1954 (p,t) TRITON in the macro this was
      // ported from, so the first gate drawn here saved as a triton with the wrong axis label
      // while the plot on screen said proton. They are parameters now.
      //
      // xaxis is "sqrt_dEdx": that is the spyral_utils Cut2D spelling AtCut2D and
      // apply_gate_C15p.C expect. "sqrtdedx" is a different string and makes the consumer warn
      // that the gate may be for another observable.
      fprintf(f, "{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrt_dEdx\",\n    \"yaxis\": \"brho\",\n",
              fName->GetText());
      fprintf(f, "    \"Z\": %d,\n    \"A\": %d,\n    \"vertices\": [\n", fZ, fA);
      int np = fNew->GetN();
      for (int i = 0; i < np; ++i) {
         double x, y;
         fNew->GetPoint(i, x, y);
         fprintf(f, "        [%.4f, %.4f]%s\n", x, y, (i == np - 1) ? "" : ",");
      }
      fprintf(f, "    ]\n}\n");
      fclose(f);
      printf("saved %d vertices (Z=%d A=%d) -> %s\n", np, fZ, fA, outPath.Data());
      fLabel->SetText(Form("  saved %d vertices, Z=%d A=%d -> %s", np, fZ, fA, outPath.Data()));
   }

   void SavePNG()
   {
      gSystem->mkdir("pid/plots", kTRUE);
      fCanvas->SaveAs("pid/plots/gate_draw_C15p.png");
      printf("saved pid/plots/gate_draw_C15p.png\n");
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
   /// Refill the displayed plane from the points already in memory. Both the angle window and the
   /// axis controls go through here, so the two cannot get out of step -- an earlier version had
   /// the angle cut refill and the axis change rebuild, and a rebuild silently dropped the angle
   /// cut that was showing on screen.
   long FillPlane()
   {
      fH->Reset();
      long kept = 0;
      for (size_t i = 0; i < fX.size(); ++i) {
         double thLab = fFlip ? 180.0 - fPol[i] : fPol[i];
         if (thLab < fThLo || thLab > fThHi) continue;
         fH->Fill(fX[i], fY[i]);
         ++kept;
      }
      return kept;
   }

   /// Rebuild the histogram on new limits or binning, then refill. The POINTS are untouched: this
   /// only changes how they are displayed, so it is free to use as often as wanted and it can
   /// never alter what a gate selects.
   void ApplyAxes()
   {
      if (fEXlo) fXlo = fEXlo->GetNumber();
      if (fEXhi) fXhi = fEXhi->GetNumber();
      if (fEYlo) fYlo = fEYlo->GetNumber();
      if (fEYhi) fYhi = fEYhi->GetNumber();
      if (fENbx) fNbx = (int)fENbx->GetNumber();
      if (fENby) fNby = (int)fENby->GetNumber();
      if (fXhi <= fXlo || fYhi <= fYlo || fNbx < 1 || fNby < 1) {
         printf("\033[1;33maxis limits must increase and bins must be >= 1 -- ignored\033[0m\n");
         return;
      }
      TString ttl = fH->GetTitle();
      delete fH;
      fH = new TH2F("hpt", ttl, fNbx, fXlo, fXhi, fNby, fYlo, fYhi);
      long kept = FillPlane();
      printf("axes: x [%.2f, %.2f] / %d bins, y [%.3f, %.3f] / %d bins -- %ld tracks shown\n",
             fXlo, fXhi, fNbx, fYlo, fYhi, fNby, kept);
      Redraw();
   }

   void ApplyAngle()
   {
      if (fEThLo) fThLo = fEThLo->GetNumber();
      if (fEThHi) fThHi = fEThHi->GetNumber();
      long kept = FillPlane();
      printf("theta_lab in [%.1f, %.1f] deg : %ld of %zu tracks kept (%.1f%%)\n",
             fThLo, fThHi, kept, fX.size(), 100.0 * kept / std::max<size_t>(1, fX.size()));
      if (fLabel)
         fLabel->SetText(Form("  theta_lab %.0f-%.0f deg : %ld of %zu tracks (%.0f%%)",
                              fThLo, fThHi, kept, fX.size(), 100.0 * kept / std::max<size_t>(1, fX.size())));
      Redraw();
   }

   /// Back to the full plane, so a zoom can always be undone without restarting.
   void ResetAxes()
   {
      fXlo = 0; fXhi = fXmax; fYlo = 0; fYhi = fYmax; fNbx = 300; fNby = 300;
      if (fEXlo) fEXlo->SetNumber(fXlo);
      if (fEXhi) fEXhi->SetNumber(fXhi);
      if (fEYlo) fEYlo->SetNumber(fYlo);
      if (fEYhi) fEYhi->SetNumber(fYhi);
      if (fENbx) fENbx->SetNumber(fNbx);
      if (fENby) fENby->SetNumber(fNby);
      ApplyAxes();
   }

   /// Show or hide the (p,d) deuteron band.
   void ToggleDeuteron()
   {
      fShowDeut = !fShowDeut;
      printf("deuteron (p,d) band %s\n", fShowDeut ? "ON (violet, dashed)" : "OFF");
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
         if (fOnLocus->GetN() > 0) fOnLocus->Draw("P same");
      }
      // the 14C(p,d)13C deuteron band, drawn so the triton gate can be placed AWAY from it rather
      // than by eye. Deuterons are the natural contaminant here: same Z, so the same charge state,
      // and a rigidity that runs through the same region of the plane.
      if (fShowDeut) {
         if (!fDeutLoc) fDeutLoc = dpLocusGraph(fFlip, fEbeam, 0.0, kViolet + 1, 2.01410178, 14.0032420);
         if (fDeutLoc && fDeutLoc->GetN() > 1) { fDeutLoc->SetLineWidth(3); fDeutLoc->SetLineStyle(2); fDeutLoc->Draw("L same"); }
      }
      // An empty TGraph makes TGraphPainter print "illegal number of points (0)" on every
      // redraw; with the locus disabled that was two errors per repaint and it buried real output.
      if (fRefP && fRefP->GetN() > 1) fRefP->Draw("L");
      if (fRefD && fRefD->GetN() > 1) fRefD->Draw("L");
      if (fNew && fNew->GetN() > 1) fNew->Draw("L");
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
      int npulse = 1;
      t->SetBranchAddress("sqrtdedx", &x);
      t->SetBranchAddress("brho", &y);
      t->SetBranchAddress("polar", &pol);
      // ★ SINGLE PULSE, to match apply_gate_C15p.C. The IC window was chosen on the single-pulse
      // spectrum and pid/ic_C15p.json records singlePulse: true, so a drawer that applies the
      // window WITHOUT this condition shows a plane 8.2 % larger than the one the gate is later
      // applied to (361340 against 333849 tracks) -- the extra being pile-up. The gate would then
      // be drawn against a background that is not there at application time. apply_gate and
      // gate_events already had to be reconciled over exactly this; the drawer was still the odd
      // one out.
      const bool hasNP = t->GetBranch("npulse") != nullptr;
      if (hasNP)
         t->SetBranchAddress("npulse", &npulse);
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
         // fIcLo < 0 disables the IC cut. Without this, the -1 default compares every real IC
         // value against a window of [-1,-1] and silently discards the whole sample, which reads
         // as "this beam has no such particles" rather than as an unset gate.
         const bool icCutOn = (fIcLo >= 0 && fIcHi > fIcLo);
         if (icCutOn && hasIC && (ic < fIcLo || ic > fIcHi))
            continue;
         // A track whose run had no IC (ic = -1) cannot satisfy a window. Dropping it is right,
         // but it must be counted rather than vanish.
         if (icCutOn && hasIC && ic < 0)
            continue;
         // Same condition the applier uses, and only when the IC window is actually on: with the
         // window disabled the drawer is deliberately showing the raw plane.
         if (icCutOn && hasNP && npulse != 1)
            continue;
         ++nCut;
         fX.push_back(x);
         fY.push_back(y);
         fPol.push_back(pol);
         fH->Fill(x, y);
         // mark the tracks sitting on the (d,p) g.s. locus. Whether the stored polar is theta_lab
         // or 180 - theta_lab is a CONVENTION, so it is an argument rather than an assumption.
         double b = dpBrhoC15p(fFlip ? 180.0 - pol : pol, fEbeam, 0.0, fM3, fM4);
         if (b > 0 && std::fabs(y - b) / b < fLocusTol) {
            if (!fOnLocus) fOnLocus = new TGraph();
            fOnLocus->SetPoint(fOnLocus->GetN(), x, y);
         }
      }
      printf("=== %s: %ld tracks, %ld inside IC [%.0f,%.0f]%s (%.1f%%) ===\n", cache.Data(), nAll, nCut,
             fIcLo, fIcHi, (fIcLo >= 0 && fIcHi > fIcLo && hasNP) ? " and single-pulse" : "",
             nAll ? 100.0 * nCut / nAll : 0.0);
      if (fIcLo >= 0 && fIcHi > fIcLo && !hasNP)
         printf("\033[1;33mno npulse branch -- pile-up NOT removed, so this plane is wider than the "
                "one apply_gate_C15p.C will select on\033[0m\n");
      f->Close();
   }

   /// Overlay a reference polygon. Shown, never enforced -- these are the EXISTING proton and
   /// deuteron gates, drawn so it is obvious which band is which before the first vertex.
   TCutG *LoadRef(TString ref, int colour, const char *what)
   {
      if (ref.IsNull()) return nullptr;
      // Read with C stdio. Two cling limitations rule out the obvious C++ idioms here:
      //   * std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
      //     fails to COMPILE ("no matching conversion for functional-style cast"), taking the
      //     whole macro down so the GUI never appears;
      //   * the ostringstream/rdbuf replacement compiles but fails to LINK at run time with
      //     "symbol '__clang_call_terminate' unresolved", killing the drawer after the plane has
      //     already been read -- which looks like a GUI crash rather than an I/O problem.
      // fread has neither failure mode.
      FILE *fp = fopen(ref.Data(), "rb");
      if (fp == nullptr) { printf("no %s at %s (skipped)\n", what, ref.Data()); return nullptr; }
      std::string s;
      char buf[4096];
      size_t got = 0;
      while ((got = fread(buf, 1, sizeof(buf), fp)) > 0)
         s.append(buf, got);
      fclose(fp);
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
      // Title follows the species the drawer was OPENED with, never a hardcoded one. The old
      // fixed "14C(p,t)12C -- triton gate" text survived the port from a1954 and said triton
      // whatever Z/A was passed, which is the same class of mislabelling that already caused a
      // proton gate to be saved under the deuteron name.
      main->SetWindowName(Form("C15p PID  --  draw the %s gate (Z=%d A=%d)",
                               (fZ == 1 && fA == 1)   ? "proton"
                               : (fZ == 1 && fA == 2) ? "deuteron"
                               : (fZ == 1 && fA == 3) ? "triton"
                               : (fZ == 2 && fA == 4) ? "alpha"
                                                      : "selected",
                               fZ, fA));
      auto *bar = new TGHorizontalFrame(main);
      const char *lbl[] = {"Draw new gate", "Evaluate", "Locus check", "Locus dots on/off",
                           "(p,d) band on/off", "Save JSON", "Save PNG", "Redraw", "Quit"};
      const char *slot[] = {"DrawGate()", "Evaluate()", "LocusCheck()", "ToggleLocus()",
                            "ToggleDeuteron()", "Save()", "SavePNG()", "Redraw()", "Quit()"};
      for (int i = 0; i < 9; ++i) {
         auto *b = new TGTextButton(bar, lbl[i]);
         b->Connect("Clicked()", "C15pGateDraw", this, slot[i]);
         bar->AddFrame(b, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      }
      bar->AddFrame(new TGLabel(bar, "name:"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 12, 2, 4, 4));
      // Default the name to the output file's stem: a gate whose "name" says one species while
      // its filename says another is how triton_14C_proton got written.
      TString stem = gSystem->BaseName(fOut.Data());
      stem.ReplaceAll(".json", "");
      fName = new TGTextEntry(bar, stem.Data());
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
      fEThLo->Connect("ValueSet(Long_t)", "C15pGateDraw", this, "ApplyAngle()");
      bar->AddFrame(fEThLo, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
      bar->AddFrame(new TGLabel(bar, "to"), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 2, 4, 4));
      fEThHi = new TGNumberEntry(bar, 180.0, 5, -1, TGNumberFormat::kNESRealOne,
                                 TGNumberFormat::kNEANonNegative,
                                 TGNumberFormat::kNELLimitMinMax, 0.0, 180.0);
      fEThHi->Resize(70, 22);
      fEThHi->Connect("ValueSet(Long_t)", "C15pGateDraw", this, "ApplyAngle()");
      bar->AddFrame(fEThHi, new TGLayoutHints(kLHintsLeft, 2, 4, 4, 4));
      {
         auto *ba = new TGTextButton(bar, "Apply angle");
         ba->Connect("Clicked()", "C15pGateDraw", this, "ApplyAngle()");
         bar->AddFrame(ba, new TGLayoutHints(kLHintsLeft, 5, 4, 4, 4));
      }
      main->AddFrame(bar, new TGLayoutHints(kLHintsTop | kLHintsExpandX));

      // ---- second row: display limits and binning -------------------------------------------
      // These change only how the points are drawn. The gate is stored as a polygon in physical
      // units, so re-binning or re-ranging can never change what it selects.
      auto *bar2 = new TGHorizontalFrame(main);
      auto mkNum = [&](TGNumberEntry *&e, double v, double lo, double hi, int w, const char *lbl,
                       TGNumberFormat::EStyle st) {
         bar2->AddFrame(new TGLabel(bar2, lbl), new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 8, 2, 4, 4));
         e = new TGNumberEntry(bar2, v, 6, -1, st, TGNumberFormat::kNEAAnyNumber,
                               TGNumberFormat::kNELLimitMinMax, lo, hi);
         e->Resize(w, 22);
         e->Connect("ValueSet(Long_t)", "C15pGateDraw", this, "ApplyAxes()");
         bar2->AddFrame(e, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 4));
      };
      mkNum(fEXlo, fXlo, 0, 200, 70, "sqrt(dE/dx):", TGNumberFormat::kNESRealTwo);
      mkNum(fEXhi, fXhi, 0, 200, 70, "to",           TGNumberFormat::kNESRealTwo);
      mkNum(fENbx, fNbx, 1, 2000, 60, "bins",        TGNumberFormat::kNESInteger);
      mkNum(fEYlo, fYlo, 0, 20, 70, "  Brho:",       TGNumberFormat::kNESRealThree);
      mkNum(fEYhi, fYhi, 0, 20, 70, "to",            TGNumberFormat::kNESRealThree);
      mkNum(fENby, fNby, 1, 2000, 60, "bins",        TGNumberFormat::kNESInteger);
      {
         auto *bx = new TGTextButton(bar2, "Apply axes");
         bx->Connect("Clicked()", "C15pGateDraw", this, "ApplyAxes()");
         bar2->AddFrame(bx, new TGLayoutHints(kLHintsLeft, 8, 4, 4, 4));
         auto *br = new TGTextButton(bar2, "Reset axes");
         br->Connect("Clicked()", "C15pGateDraw", this, "ResetAxes()");
         bar2->AddFrame(br, new TGLayoutHints(kLHintsLeft, 4, 4, 4, 4));
      }
      main->AddFrame(bar2, new TGLayoutHints(kLHintsTop | kLHintsExpandX));
      fLabel = new TGLabel(main, "  Draw around the band high and right of the proton gate, then Locus check. "
                                 "Violet dashed = the (p,d) deuteron band: keep the gate off it.");
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
   TGraph *fDeutLoc = nullptr;
   bool fShowDeut = true;
   TGNumberEntry *fEXlo = nullptr, *fEXhi = nullptr, *fEYlo = nullptr, *fEYhi = nullptr;
   TGNumberEntry *fENbx = nullptr, *fENby = nullptr;
   double fXlo = 0.0, fXhi = 30.0, fYlo = 0.0, fYhi = 3.0;
   int fNbx = 300, fNby = 300;
   TGLabel *fLabel = nullptr;
   TGTextEntry *fName = nullptr;
   int fZ = 1, fA = 1;   ///< species written into the gate JSON
   double fM3 = 1.00782503, fM4 = 16.0147413; ///< ejectile/residual for THIS species' locus
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
void gate_draw_C15p(TString outJson = "pid/proton_C15p.json",
                       TString cache = "pid/points_C15p.root",
                       TString refP = "", TString refD = "",
                       double xMax = 60.0, double yMax = 2.0, double icLo = -1, double icHi = -1,
                       bool showLocus = false, double eBeam = 0, bool flipPolar = true,
                       /// Species written into the gate. Defaults to the PROTON, the (d,p) ejectile.
                       ///   deuteron: Z 1 A 2      triton: Z 1 A 3      alpha: Z 2 A 4
                       int Z = 1, int A = 1)
{
   new C15pGateDraw(outJson, cache, refP, refD, xMax, yMax, icLo, icHi, showLocus, eBeam, flipPolar, Z, A);
}
