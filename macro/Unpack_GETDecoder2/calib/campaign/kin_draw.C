// Dec 2014 alpha+alpha kinematics drawer -- ONE window.
//
// Plane: lab angle to the beam vs energy from range, with the elastic loci overlaid.
// v_D is a LIVE control: every track is stored as NEAR/FAR endpoints in (x, y, tb), so
// changing the drift velocity rebuilds every angle and every range. That is the point --
// the open question is whether the measured band can be put on the alpha+alpha locus by
// any v_D, and this lets that be judged directly.
//
// One window by construction: the locus-residual histogram goes into a sub-pad of the same
// canvas, not a second TCanvas.

#include <fstream>
#include <vector>

TTree *kT = nullptr;
TCanvas *kC = nullptr;
TH2F *kH = nullptr;
TH1F *kRes = nullptr;
TGNumberEntry *eVD, *eEB, *eP, *eOpLo, *eOpHi, *eNx, *eNy;
TGComboBox *cTgt;
TGCheckButton *kCurv, *kResOn;
TGLabel *kStat;
std::vector<double> *tD = nullptr, *tE = nullptr;
const double kT0 = -34.0;

void kinLoadTable(const char *p)
{
   if (!tD) { tD = new std::vector<double>(); tE = new std::vector<double>(); }
   std::ifstream f(p); std::string s;
   while (std::getline(f, s)) {
      if (s.empty() || s[0] == '#') continue;
      std::stringstream ss(s); double a, b; ss >> a >> b;
      tD->push_back(a); tE->push_back(b);
   }
   printf("  range table %zu points, R(7.80)=%.0f mm @150 torr\n", tD->size(), tD->back());
}
double kinEofR(double R, double press)
{
   double d = tD->back() - R * press / 150.0;
   if (d <= 0) return tE->front();
   if (d >= tD->back()) return 0.0;
   Int_t i = TMath::BinarySearch((Long64_t)tD->size(), &(*tD)[0], d);
   if (i < 0) i = 0;
   if (i >= (Int_t)tD->size() - 1) return tE->back();
   double f = (d - (*tD)[i]) / ((*tD)[i+1] - (*tD)[i]);
   return (*tE)[i] + f * ((*tE)[i+1] - (*tE)[i]);
}
double kinEscat(double E, double th, double m, double M)
{
   double s = TMath::Sin(th*TMath::DegToRad()), c = TMath::Cos(th*TMath::DegToRad());
   double disc = M*M - m*m*s*s;
   if (disc < 0) return -1;
   return E * TMath::Power((m*c + TMath::Sqrt(disc)) / (m+M), 2);
}
double kinErec(double E, double th, double m, double M)
{
   if (th > 90) return -1;
   double c = TMath::Cos(th*TMath::DegToRad());
   return 4*m*M/((m+M)*(m+M)) * E * c * c;
}
double kinNum(TGNumberEntry *e, double d) { return e ? e->GetNumber() : d; }
int    kinTgt() { return cTgt ? cTgt->GetSelected() : 0; }
double kinMass() { int t = kinTgt(); return t == 1 ? 12.0 : (t == 2 ? 16.0 : 4.0); }
const char *kinTname() { int t = kinTgt(); return t == 1 ? "^{12}C" : (t == 2 ? "^{16}O" : "#alpha"); }

void kinRedraw()
{
   const double vD = kinNum(eVD,2.142), EB = kinNum(eEB,7.80), pr = kinNum(eP,299.5);
   const double oLo = kinNum(eOpLo,80), oHi = kinNum(eOpHi,100);
   const Int_t nx = (Int_t)kinNum(eNx,180), ny = (Int_t)kinNum(eNy,150);
   const bool wantRes = kResOn && kResOn->IsOn();

   if (kH) { delete kH; kH = nullptr; }
   kH = new TH2F("kH", Form("v_{D} = %.3f cm/#mus, %.0f torr, E_{b} = %.2f MeV, opening %.0f-%.0f#circ;"
                            "lab angle to the beam  #theta [deg];energy from range [MeV]",
                            vD, pr, EB, oLo, oHi), nx, 0, 180, ny, 0, 9);
   kH->SetStats(0);
   if (kRes) { delete kRes; kRes = nullptr; }
   kRes = new TH1F("kRes", Form("residual to the #alpha+%s locus;E_{meas} - E_{kin} [MeV];tracks",
                                kinTname()), 100, -4, 6);
   kRes->SetStats(0);

   Int_t run; Float_t b[6], a1[6], a2[6];
   kT->SetBranchAddress("run",&run); kT->SetBranchAddress("b",b);
   kT->SetBranchAddress("a1",a1);    kT->SetBranchAddress("a2",a2);
   const double K = vD * 160e-3 * 10.0;
   const double M = kinMass();
   Long64_t n = kT->GetEntries(), kept = 0;
   for (Long64_t i = 0; i < n; ++i) {
      kT->GetEntry(i);
      TVector3 bn(b[0],b[1],K*(b[2]-kT0)), bf(b[3],b[4],K*(b[5]-kT0));
      TVector3 bd = (bn-bf).Unit();                 // entrance -> vertex = direction of travel
      TVector3 u[2]; double L[2];
      const Float_t *A[2] = {a1,a2};
      for (int k = 0; k < 2; ++k) {
         TVector3 nr(A[k][0],A[k][1],K*(A[k][2]-kT0)), fr(A[k][3],A[k][4],K*(A[k][5]-kT0));
         u[k] = (fr-nr).Unit(); L[k] = (fr-nr).Mag();
      }
      double op = TMath::ACos(TMath::Max(-1.,TMath::Min(1.,u[0].Dot(u[1])))) * TMath::RadToDeg();
      if (op < oLo || op > oHi) continue;
      ++kept;
      for (int k = 0; k < 2; ++k) {
         double th = TMath::ACos(TMath::Max(-1.,TMath::Min(1.,u[k].Dot(bd)))) * TMath::RadToDeg();
         double E = kinEofR(L[k], pr);
         kH->Fill(th, E);
         double ek = kinEscat(EB, th, 4.0, M);
         if (ek > 0) kRes->Fill(E - ek);
      }
   }

   kC->cd(); kC->Clear();
   if (wantRes) { kC->Divide(1,2); kC->cd(1); gPad->SetPad(0,0.32,1,1); }
   gPad->SetRightMargin(0.13); gPad->SetLogz();
   kH->Draw("COLZ");
   if (!kCurv || kCurv->IsOn()) {
      const int NP = 400;
      TGraph *gs = new TGraph(), *gr = new TGraph();
      for (int i = 0; i < NP; ++i) {
         double th = 0.5 + i*(179.0/NP);
         double es = kinEscat(EB,th,4.0,M); if (es > 0) gs->SetPoint(gs->GetN(), th, es);
         double er = kinErec(EB,th,4.0,M);  if (er > 0) gr->SetPoint(gr->GetN(), th, er);
      }
      gs->SetLineColor(kWhite); gs->SetLineWidth(3); gs->Draw("L SAME");
      gr->SetLineColor(kOrange+7); gr->SetLineWidth(2); gr->SetLineStyle(2); gr->Draw("L SAME");
      static const double kFrac[3] = {0.6,0.35,0.15};
      for (int fi = 0; fi < 3; ++fi) { double f = kFrac[fi];
         TGraph *g2 = new TGraph();
         for (int i = 0; i < NP; ++i) {
            double th = 0.5 + i*(179.0/NP), es = kinEscat(EB*f,th,4.0,M);
            if (es > 0) g2->SetPoint(g2->GetN(), th, es);
         }
         g2->SetLineColor(kGray+2); g2->SetLineWidth(1); g2->Draw("L SAME");
      }
      TLegend *lg = new TLegend(0.58,0.74,0.86,0.88);
      lg->SetFillStyle(0); lg->SetBorderSize(0); lg->SetTextColor(kWhite); lg->SetTextSize(0.030);
      lg->AddEntry(gs, Form("#alpha+%s scattered", kinTname()), "l");
      lg->AddEntry(gr, Form("%s recoil", kinTname()), "l");
      lg->Draw();
   }
   if (wantRes) {
      kC->cd(2); gPad->SetPad(0,0,1,0.32); gPad->SetTopMargin(0.02); gPad->SetBottomMargin(0.28);
      kRes->SetFillColorAlpha(kAzure+1,0.55); kRes->SetLineColor(kAzure+2);
      kRes->GetXaxis()->SetTitleSize(0.10); kRes->GetXaxis()->SetLabelSize(0.09);
      kRes->GetYaxis()->SetLabelSize(0.09);
      kRes->Draw("HIST");
      TLine *l0 = new TLine(0,0,0,kRes->GetMaximum()*1.05);
      l0->SetLineColor(kRed); l0->SetLineWidth(2); l0->Draw();
   }
   kC->Update();
   if (kStat) kStat->SetText(Form("%lld events | locus residual mean %.2f rms %.2f MeV",
                                  kept, kRes->GetMean(), kRes->GetRMS()));
}
void kinSavePNG() { kC->SaveAs("kin_draw.png"); if (kStat) kStat->SetText("wrote kin_draw.png"); }
void kinQuit()    { gApplication->Terminate(0); }

// ---------------------------------------------------------------------------------------
// GUARD. If the interpreter re-executes this macro (which it does when a cling error is hit
// after the macro has run), every re-entry would map another window. This makes every call
// after the first a no-op, so at most ONE window can ever exist per process.
bool kinAlreadyRan = false;
// ---------------------------------------------------------------------------------------

void kin_draw(const char *file = "kin_points.root",
              const char *table = "/home/yassid/dec2014_kin/eloss_alpha_heco2_11p4.txt",
              double vD = 2.142, double press = 299.5, double Ebeam = 11.40)
{
   // File lock, not a global: the respawn creates a FRESH PROCESS each time, so a
   // process-local flag is reset every cycle and cannot stop it.
   if (!gSystem->AccessPathName("/tmp/kin_draw.lock")) {
      printf("kin_draw: /tmp/kin_draw.lock exists -- refusing to start. rm it to allow.\n");
      return;
   }
   { FILE *lk = fopen("/tmp/kin_draw.lock", "w"); if (lk) { fprintf(lk, "%d\n", gSystem->GetPid()); fclose(lk); } }
   kinAlreadyRan = true;
   kinLoadTable(table);
   TFile *f = TFile::Open(file);
   if (!f || f->IsZombie()) { printf("ERROR cannot open %s\n", file); return; }
   kT = (TTree *)f->Get("kin");
   if (!kT) { printf("ERROR no 'kin' tree\n"); return; }
   printf("  %lld events\n", kT->GetEntries());
   printf("CP1 palette\n"); fflush(stdout);
   gStyle->SetPalette(kBird);

   if (gROOT->IsBatch()) {                       // headless check: render once, write a PNG
      kC = new TCanvas("kC","kin",1000,700);
      kinRedraw(); kC->SaveAs("kin_draw.png");
      printf("batch: wrote kin_draw.png\n"); return;
   }

   printf("CP2 before TGMainFrame\n"); fflush(stdout);
   TGMainFrame *mf = new TGMainFrame(gClient->GetRoot(), 1180, 780);
   printf("CP3 TGMainFrame ok\n"); fflush(stdout);
   mf->SetWindowName("Dec2014 alpha kinematics");
   TGHorizontalFrame *hf = new TGHorizontalFrame(mf);
   TRootEmbeddedCanvas *ec = new TRootEmbeddedCanvas("ec", hf, 870, 720);
   hf->AddFrame(ec, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY,2,2,2,2));
   kC = ec->GetCanvas();
   printf("CP4 canvas ok\n"); fflush(stdout);

   TGVerticalFrame *vf = new TGVerticalFrame(hf, 280);
   auto addNum = [&](const char *lab, double val, double lo, double hi, TGNumberEntry *&e,
                     TGNumberFormat::EStyle st = TGNumberFormat::kNESRealThree) {
      TGHorizontalFrame *r = new TGHorizontalFrame(vf);
      r->AddFrame(new TGLabel(r, lab), new TGLayoutHints(kLHintsLeft|kLHintsCenterY,2,4,3,3));
      e = new TGNumberEntry(r, val, 8, -1, st, TGNumberFormat::kNEAAnyNumber,
                           TGNumberFormat::kNELLimitMinMax, lo, hi);
      r->AddFrame(e, new TGLayoutHints(kLHintsRight,2,2,3,3));
      vf->AddFrame(r, new TGLayoutHints(kLHintsExpandX));
   };
   addNum("v_D  [cm/us]", vD, 0.3, 5.0, eVD);
   addNum("E_beam [MeV]", Ebeam, 0.1, 30.0, eEB);   // 11.40 measured, NOT the log's 7.80
   addNum("pressure [torr]", press, 10, 900, eP, TGNumberFormat::kNESRealOne);
   addNum("opening min", 80.0, 0, 180, eOpLo, TGNumberFormat::kNESRealOne);
   addNum("opening max", 100.0, 0, 180, eOpHi, TGNumberFormat::kNESRealOne);
   addNum("bins theta", 180, 20, 720, eNx, TGNumberFormat::kNESInteger);
   addNum("bins E", 150, 20, 720, eNy, TGNumberFormat::kNESInteger);

   TGHorizontalFrame *rt = new TGHorizontalFrame(vf);
   rt->AddFrame(new TGLabel(rt,"target"), new TGLayoutHints(kLHintsLeft|kLHintsCenterY,2,4,3,3));
   cTgt = new TGComboBox(rt, 100);
   cTgt->AddEntry("alpha",0); cTgt->AddEntry("12C",1); cTgt->AddEntry("16O",2);
   cTgt->Select(0); cTgt->Resize(110,22);
   rt->AddFrame(cTgt, new TGLayoutHints(kLHintsRight,2,2,3,3));
   vf->AddFrame(rt, new TGLayoutHints(kLHintsExpandX));

   kCurv = new TGCheckButton(vf,"kinematic curves"); kCurv->SetOn();
   vf->AddFrame(kCurv, new TGLayoutHints(kLHintsLeft,4,2,6,2));
   kResOn = new TGCheckButton(vf,"locus residual panel");
   vf->AddFrame(kResOn, new TGLayoutHints(kLHintsLeft,4,2,2,8));

   auto addBtn = [&](const char *lab, const char *slot) {
      TGTextButton *b = new TGTextButton(vf, lab);
      b->Connect("Clicked()","",0,slot);
      vf->AddFrame(b, new TGLayoutHints(kLHintsExpandX,4,4,2,2));
   };
   addBtn("Redraw",   "kinRedraw()");
   addBtn("Save PNG", "kinSavePNG()");
   addBtn("Quit",     "kinQuit()");
   kStat = new TGLabel(vf, "ready");
   vf->AddFrame(kStat, new TGLayoutHints(kLHintsExpandX,4,4,10,4));

   hf->AddFrame(vf, new TGLayoutHints(kLHintsRight|kLHintsExpandY,2,2,2,2));
   mf->AddFrame(hf, new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));
   printf("CP5 before map\n"); fflush(stdout);
   mf->MapSubwindows(); mf->Resize(mf->GetDefaultSize()); mf->MapWindow();
   printf("CP6 mapped\n"); fflush(stdout);
   kinRedraw();

   // Enter the event loop HERE, inside the macro. Do not pass -e 'gApplication->Run()' on
   // the command line: that form dies with "symbol '__clang_call_terminate' unresolved"
   // AFTER the macro has run, and if stdin is being held open the interpreter restarts and
   // re-executes the macro, mapping a new window every cycle. Blocking here means ROOT never
   // reaches the prompt, so stdin can simply be /dev/null.
   printf("CP7 entering event loop\n"); fflush(stdout);
   gApplication->Run();
}
