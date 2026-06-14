/// @file pid_gate_edit.C
/// @brief Interactive PID gate editor (TGMainFrame GUI). Shows the Spyral PID
/// landscape (sqrt(dE/dx) vs Brho) with the CURRENT gate polygon overlaid, and lets
/// you draw a NEW gate by clicking on the plot, count how many tracks fall inside it,
/// and save it to JSON in the AtCut2D format that fitGenfitter_a1975.C / AtGenfitter
/// SetPIDGate read. 2D TCanvas + TGMainFrame (works over X11/WSLg; no Eve/OpenGL).
///
///   root -l 'pp/pid_gate_edit.C'                 // runs 106-125, current deuteron gate
///   root -l 'pp/pid_gate_edit.C(106,189,"pid/deuteron_band.json","pid/deuteron_tight.json")'
///
/// Workflow: [Draw new gate] -> click polygon vertices on the plot, DOUBLE-CLICK to
/// close -> [Count inside] -> [Save JSON]. The new (green) gate is what gets saved.

#include <vector>

class PIDGateEditor : public TObject {
public:
   PIDGateEditor(int runLo,int runHi,TString dir,TString suffix,TString gateJson,TString outJson)
      : fOut(outJson)
   {
      gSystem->Load("libAtReconstruction.so");
      fH=new TH2F("hpid_e","^{16}C+p Spyral PID  (draw a gate to refine);sqrt(dE/dx) [arb];B#rho [T#upointm]",300,0,30,250,0,2.5);
      long n=0;
      for(int r=runLo;r<=runHi;++r){
         TString gf=Form("%srun_%04d%s.root",dir.Data(),r,suffix.Data());
         if(gSystem->AccessPathName(gf)) continue;
         TFile*f=TFile::Open(gf); auto*t=(TTree*)f->Get("cbmsim");
         TClonesArray*pe=nullptr; t->SetBranchAddress("AtPIDEvent",&pe);
         for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
            if(pe->GetEntries()==0) continue; auto*ev=(AtPIDEvent*)pe->At(0); if(!ev) continue;
            for(auto&s:ev->GetSpyral()){ if(!s.valid) continue;
               fH->Fill(s.sqrtdEdx,s.brho); fX.push_back(s.sqrtdEdx); fY.push_back(s.brho); ++n; } }
         f->Close();
      }
      std::cout<<"PIDGateEditor: loaded "<<n<<" PID points (runs "<<runLo<<"-"<<runHi<<")\n";
      // current gate
      auto pid=AtTools::AtParticleID::LoadJSON(gateJson.Data());
      if(pid.GetCut().IsValid()){ const auto&v=pid.GetCut().GetVertices(); int np=v.size();
         fCur=new TPolyLine(np+1); for(int i=0;i<np;++i) fCur->SetPoint(i,v[i].first,v[i].second);
         fCur->SetPoint(np,v[0].first,v[0].second); fCur->SetLineColor(kRed+1); fCur->SetLineWidth(3); }
      MakeGui();
      Redraw();
   }

   void Redraw(){
      fCanvas->cd(); fCanvas->SetLogz(); fCanvas->SetRightMargin(0.13);
      fH->Draw("colz");
      if(fCur) fCur->Draw("L");
      if(fNew) fNew->Draw("L");
      auto*tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.03);
      tx->SetTextColor(kRed+1);   tx->DrawLatex(0.16,0.87,"current gate");
      tx->SetTextColor(kGreen+2); tx->DrawLatex(0.16,0.83, fNew?"new gate (will be saved)":"new gate: none yet");
      fCanvas->Modified(); fCanvas->Update(); gSystem->ProcessEvents();
   }

   void DrawGate(){
      fCanvas->cd();
      std::cout<<"\n>>> Click polygon vertices on the PID plot. DOUBLE-CLICK to close it.\n";
      TCutG*c=(TCutG*)fCanvas->WaitPrimitive("CUTG","CutG");
      if(!c){ std::cout<<"no cut drawn\n"; return; }
      fNew=(TCutG*)c->Clone("newgate"); fNew->SetLineColor(kGreen+2); fNew->SetLineWidth(3);
      delete c; // remove the transient CUTG; we keep our clone
      CountInside();
      Redraw();
   }

   void CountInside(){
      const TCutG*g=fNew?fNew:nullptr; if(!g){ std::cout<<"draw a new gate first\n"; return; }
      long in=0; for(size_t i=0;i<fX.size();++i) if(g->IsInside(fX[i],fY[i])) ++in;
      fLabel->SetText(Form("new gate: %ld / %zu PID points inside (%.1f%%)",in,fX.size(),100.0*in/std::max((size_t)1,fX.size())));
      std::cout<<"inside new gate: "<<in<<" / "<<fX.size()<<"\n";
   }

   void Save(){
      if(!fNew){ std::cout<<"no new gate drawn — nothing to save\n"; return; }
      TString name=fName->GetText();
      FILE*f=fopen(fOut.Data(),"w"); if(!f){ std::cout<<"cannot write "<<fOut<<"\n"; return; }
      fprintf(f,"{\n    \"name\": \"%s\",\n    \"xaxis\": \"sqrtdedx\",\n    \"yaxis\": \"brho\",\n    \"vertices\": [\n",name.Data());
      int np=fNew->GetN();
      for(int i=0;i<np;++i){ double x,y; fNew->GetPoint(i,x,y);
         fprintf(f,"        [%.3f, %.3f]%s\n",x,y,(i==np-1)?"":","); }
      fprintf(f,"    ]\n}\n"); fclose(f);
      std::cout<<"saved new gate ("<<np<<" vertices) -> "<<fOut<<"\n";
      fLabel->SetText(Form("saved %d-vertex gate -> %s",np,fOut.Data()));
   }
   void SavePNG(){ TString p="/tmp/pid_gate_edit.png"; fCanvas->SaveAs(p); std::cout<<"saved "<<p<<"\n"; }

private:
   void MakeGui(){
      auto*main=new TGMainFrame(gClient->GetRoot(),1100,900); main->SetWindowName("PID gate editor");
      auto*bar=new TGHorizontalFrame(main);
      auto*bD=new TGTextButton(bar,"  Draw new gate  "); bD->Connect("Clicked()","PIDGateEditor",this,"DrawGate()");
      bar->AddFrame(bD,new TGLayoutHints(kLHintsLeft,6,4,4,4));
      auto*bC=new TGTextButton(bar,"  Count inside  "); bC->Connect("Clicked()","PIDGateEditor",this,"CountInside()");
      bar->AddFrame(bC,new TGLayoutHints(kLHintsLeft,4,4,4,4));
      auto*bS=new TGTextButton(bar,"  Save JSON  "); bS->Connect("Clicked()","PIDGateEditor",this,"Save()");
      bar->AddFrame(bS,new TGLayoutHints(kLHintsLeft,4,4,4,4));
      auto*bP=new TGTextButton(bar,"  Save PNG  "); bP->Connect("Clicked()","PIDGateEditor",this,"SavePNG()");
      bar->AddFrame(bP,new TGLayoutHints(kLHintsLeft,4,4,4,4));
      bar->AddFrame(new TGLabel(bar,"name:"),new TGLayoutHints(kLHintsLeft|kLHintsCenterY,12,2,4,4));
      fName=new TGTextEntry(bar,"deuteron_band"); fName->Resize(140,22);
      bar->AddFrame(fName,new TGLayoutHints(kLHintsLeft,2,4,4,4));
      main->AddFrame(bar,new TGLayoutHints(kLHintsTop|kLHintsExpandX));
      fLabel=new TGLabel(main,"  draw a polygon on the deuteron band; double-click to close, then Count / Save  ");
      main->AddFrame(fLabel,new TGLayoutHints(kLHintsTop|kLHintsExpandX,6,6,2,2));
      auto*ec=new TRootEmbeddedCanvas("ec_pidedit",main,1080,820);
      main->AddFrame(ec,new TGLayoutHints(kLHintsExpandX|kLHintsExpandY));
      fCanvas=ec->GetCanvas();
      main->MapSubwindows(); main->Resize(main->GetDefaultSize()); main->MapWindow();
   }
   std::vector<float> fX,fY;
   TH2F*fH{nullptr}; TPolyLine*fCur{nullptr}; TCutG*fNew{nullptr};
   TCanvas*fCanvas{nullptr}; TGTextEntry*fName{nullptr}; TGLabel*fLabel{nullptr}; TString fOut;
   ClassDef(PIDGateEditor,0);
};

void pid_gate_edit(int runLo=106,int runHi=125, TString gateJson="pid/deuteron_band.json",
                   TString outJson="pid/deuteron_tight.json", TString dir="/mnt/f/a1975/reco_pd/",
                   TString suffix="_genfitter_pd")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   new PIDGateEditor(runLo,runHi,dir,suffix,gateJson,outJson);
}
