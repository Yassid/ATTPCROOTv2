/// @file eval_gate_pt_Be12.C
/// @brief Evaluate a saved 12Be(p,t) triton gate: what it keeps, what it steals from the proton and
///        deuteron gates, whether the selection lands on the (p,t) kinematic locus, and WHICH BEAM
///        it came from (the IC panel -- the 12Be peak is the small one at 625-750, the dominant
///        ~1900 component is a contaminant).
///
///   root -l 'pid/eval_gate_pt_Be12.C'              // live, zoomable
///   root -b -q 'pid/eval_gate_pt_Be12.C'           // just write the png
///
/// NOTE: never call c->Draw() on the TCanvas itself -- it draws the canvas INTO the current pad and
/// blanks every sub-pad. The canvas is already on screen once it is constructed.
#include "gate_draw_pt_Be12.C"

static TCutG *EG_load(TString f, const char *n, int col)
{
   std::ifstream in(f.Data()); if (!in) return nullptr;
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto pos = s.find('[', s.find("vertices")); if (pos == std::string::npos) return nullptr;
   std::vector<double> v; const char *p = s.c_str() + pos, *e = s.c_str() + s.size(); int d = 0;
   while (p < e) { if (*p=='[') d++; if (*p==']') { d--; if (d<=0) break; }
      char *np=nullptr; double x=strtod(p,&np); if (np!=p) { v.push_back(x); p=np; } else ++p; }
   if (v.size() < 6) return nullptr;
   auto *c = new TCutG(n, v.size()/2);
   for (size_t i=0,k=0; i+1<v.size(); i+=2,++k) c->SetPoint(k, v[i], v[i+1]);
   c->SetLineColor(col); c->SetLineWidth(3); return c;
}

void eval_gate_pt_Be12(TString gateFile = "pid/triton_12Be.json", double icLo = 500, double icHi = 800,
                       double eBeam = 155.0)
{
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TString d = "pid/";
   TCutG *T = EG_load(gateFile, "Tg", kBlue+1);
   if (!T) { printf("\033[1;31mno gate at %s\033[0m\n", gateFile.Data()); return; }
   TCutG *P = EG_load(d+"proton_12Be.json", "Pg", kRed+1);
   TCutG *D = EG_load(d+"deuteron_12Be.json", "Dg", kGreen+2);

   TFile *f = TFile::Open(d+"plots/points_Be12_mp30.root");
   TTree *t = (TTree*)f->Get("pts");
   Float_t sd,br,pol,arc,ic;
   t->SetBranchAddress("sqrtdedx",&sd); t->SetBranchAddress("brho",&br);
   t->SetBranchAddress("polar",&pol);   t->SetBranchAddress("arclen",&arc); t->SetBranchAddress("ic",&ic);

   auto *hPlane = new TH2F("hPlane","IC-gated PID plane + gates (blue = triton);#sqrt{dE/dx};B#rho [T m]",300,0,30,300,0,3);
   auto *hLoc   = new TH2F("hLoc","B#rho vs polar INSIDE the gate;#theta_{polar} [deg];B#rho [T m]",180,0,180,150,0,3);
   auto *hZoom  = new TH2F("hZoom","zoom 120-180 deg, INSIDE the gate;#theta_{polar} [deg];B#rho [T m]",60,120,180,120,0.4,2.8);
   auto *hPol   = new TH1F("hPol","polar of the selected tracks;#theta_{polar} [deg];tracks",180,0,180);
   auto *hIcAll = new TH1F("hIcAll","IC amplitude;IC max;tracks",300,0,3000);
   auto *hIcSel = new TH1F("hIcSel","IC of the tracks in the triton gate (NO IC cut applied);IC max;tracks",300,0,3000);

   long nIC=0,in=0,ovP=0,ovD=0,onloc=0,selAllIC=0;
   for (Long64_t i=0;i<t->GetEntries();++i) {
      t->GetEntry(i);
      if (ic >= 0) hIcAll->Fill(ic);
      // the gate's OWN IC spectrum, deliberately WITHOUT the IC cut: this is what says whether the
      // band is 12Be at all, and it is the panel that was missing.
      bool inGate = T->IsInside(sd,br);
      if (inGate && ic >= 0) { hIcSel->Fill(ic); ++selAllIC; }
      if (ic < icLo || ic > icHi) continue;
      ++nIC; hPlane->Fill(sd,br);
      if (!inGate) continue;
      ++in; hLoc->Fill(pol,br); hZoom->Fill(pol,br); hPol->Fill(pol);
      if (P && P->IsInside(sd,br)) ++ovP;
      if (D && D->IsInside(sd,br)) ++ovD;
      double b = ptBrhoAt(180.0-pol, eBeam);
      if (b>0 && std::fabs(br-b)/b < 0.15) ++onloc;
   }
   printf("\n\033[1;33m=== %s ===\033[0m\n", gateFile.Data());
   printf("IC-gated tracks [%g,%g] : %ld\n", icLo, icHi, nIC);
   printf("inside the triton gate  : %ld  (%.1f%% of IC-gated)\n", in, nIC?100.*in/nIC:0);
   printf("  also in proton gate   : %ld\n  also in deuteron gate : %ld\n", ovP, ovD);
   printf("  on the (p,t) g.s.locus: %ld  (%.1f%%)\n", onloc, in?100.*onloc/in:0);
   printf("in the gate, ANY beam   : %ld  -> the IC cut keeps %.1f%% of them\n",
          selAllIC, selAllIC?100.*in/selAllIC:0);

   auto *c = new TCanvas("cEval","(p,t) gate evaluation",1600,950);
   c->Divide(3,2);
   c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.13); hPlane->Draw("colz");
   if(P)P->Draw("l"); if(D)D->Draw("l"); T->Draw("l");
   c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.13); hLoc->Draw("colz");
   auto *g1 = ptLocusGraph(true,eBeam,0.0,kRed+1);      g1->Draw("L same");
   auto *g2 = ptLocusGraph(true,eBeam,3.368,kMagenta+1); g2->SetLineStyle(7); g2->Draw("L same");
   auto *lg = new TLegend(0.13,0.74,0.55,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(g1,"(p,t) g.s.","l"); lg->AddEntry(g2,"^{10}Be 3.368 (2^{+})","l"); lg->Draw();
   c->cd(3); gPad->SetLogz(); gPad->SetRightMargin(0.13); hZoom->Draw("colz");
   auto *g3 = ptLocusGraph(true,eBeam,0.0,kRed+1);       g3->Draw("L same");
   auto *g4 = ptLocusGraph(true,eBeam,3.368,kMagenta+1); g4->SetLineStyle(7); g4->Draw("L same");
   c->cd(4); gPad->SetLogy(); hIcAll->SetLineColor(kGray+2); hIcAll->Draw("hist");
   auto *bx = new TBox(icLo,0.5,icHi,hIcAll->GetMaximum()); bx->SetFillColorAlpha(kGreen-9,0.35);
   bx->Draw(); hIcAll->Draw("hist same");
   auto *l2 = new TLegend(0.40,0.76,0.88,0.88); l2->SetBorderSize(0); l2->SetFillStyle(0);
   l2->AddEntry(hIcAll,"all tracks","l"); l2->AddEntry(bx,Form("gate %g-%g",icLo,icHi),"f"); l2->Draw();
   c->cd(5); gPad->SetLogy(); hIcSel->SetLineColor(kBlue+1); hIcSel->SetLineWidth(2); hIcSel->Draw("hist");
   auto *bx2 = new TBox(icLo,0.5,icHi,hIcSel->GetMaximum()); bx2->SetFillColorAlpha(kGreen-9,0.35);
   bx2->Draw(); hIcSel->Draw("hist same");
   c->cd(6); hPol->SetLineColor(kBlue+1); hPol->SetLineWidth(2); hPol->Draw("hist");
   c->Modified(); c->Update();
   c->SaveAs("/home/yassid/a1954_Be12_pt_gate_eval.png");
}
