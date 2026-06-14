/// @file pid_pd.C
/// @brief Plot the Spyral PID landscape (sqrt(dE/dx) vs Brho) from the a1975 genfit
/// output and overlay the deuteron gate polygon, to see band separation and any
/// triton/proton leakage into the deuteron gate. Reads AtPIDEvent (computed for ALL
/// tracks) from <run>_genfitter_pd.root.
///
///   root -l -b -q 'pp/pid_pd.C'                       // first 8 runs, default gate
///   root -l -b -q 'pp/pid_pd.C(106,189)'              // run range

void pid_pd(int runLo=106, int runHi=113, TString dir="/mnt/f/a1975/reco_pd/",
            TString suffix="_genfitter_pd", TString gateJson="pid/deuteron_band.json",
            TString out="/tmp/pid_pd.png", double icMin=950, double icMax=1350,
            TString fribDir="/mnt/f/a1975/reco/", int icTbLo=1000, int icTbHi=1350)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   bool useIC = (icMin>=0);

   auto*h=new TH2F("hpid",Form("^{16}C+p Spyral PID%s;sqrt(dE/dx)  [arb];B#rho  [T#upointm]",
                   useIC?" (IC beam-gated)":""),300,0,30,250,0,2.5);
   long ntot=0;
   for(int r=runLo;r<=runHi;++r){
      TString gf=Form("%srun_%04d%s.root",dir.Data(),r,suffix.Data());
      TString ff=Form("%srun_%04d_FRIB.root",fribDir.Data(),r);
      if(gSystem->AccessPathName(gf)) continue;
      bool haveFrib = useIC && !gSystem->AccessPathName(ff);
      TFile*f=TFile::Open(gf); auto*t=(TTree*)f->Get("cbmsim");
      TClonesArray*pe=nullptr; t->SetBranchAddress("AtPIDEvent",&pe);
      TFile*fc=haveFrib?TFile::Open(ff):nullptr; TTree*tc=haveFrib?(TTree*)fc->Get("cbmsim"):nullptr;
      TClonesArray*re=nullptr; if(tc) tc->SetBranchAddress("AtRawEvent",&re);
      Long64_t N = tc ? std::min(t->GetEntries(),tc->GetEntries()) : t->GetEntries();
      for(Long64_t i=0;i<N;++i){
         if(tc){ tc->GetEntry(i); double ic=-1; // ion-chamber 16C beam gate
            if(re->GetEntries()>0){ auto*raw=(AtRawEvent*)re->At(0);
               if(raw&&!raw->GetGenTraces().empty()){ auto&adc=raw->GetGenTraces()[0]->GetADC();
                  double mx=-1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();++b) mx=std::max(mx,adc[b]); ic=mx; } }
            if(ic<icMin||ic>icMax) continue; }
         t->GetEntry(i);
         if(pe->GetEntries()==0) continue; auto*ev=(AtPIDEvent*)pe->At(0); if(!ev) continue;
         for(auto&s:ev->GetSpyral()){ if(!s.valid) continue; h->Fill(s.sqrtdEdx,s.brho); ++ntot; } }
      f->Close(); if(fc) fc->Close();
   }
   printf("filled %ld valid PID points from runs %d-%d%s\n",ntot,runLo,runHi,useIC?" (IC-gated)":"");

   // deuteron gate polygon from the JSON (parsed by AtParticleID)
   TPolyLine*gate=nullptr;
   auto pid=AtTools::AtParticleID::LoadJSON(gateJson.Data());
   if(pid.GetCut().IsValid()){
      const auto&v=pid.GetCut().GetVertices(); int np=v.size();
      gate=new TPolyLine(np+1);
      for(int i=0;i<np;++i) gate->SetPoint(i,v[i].first,v[i].second);
      gate->SetPoint(np,v[0].first,v[0].second);
      gate->SetLineColor(kRed+1); gate->SetLineWidth(3);
   } else printf("WARN: could not load gate %s\n",gateJson.Data());

   auto*c=new TCanvas("c_pid","PID",1000,800);
   c->SetLogz(); c->SetRightMargin(0.13);
   h->Draw("colz");
   if(gate) gate->Draw("L");
   auto*tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.03); tx->SetTextColor(kRed+1);
   tx->DrawLatex(0.16,0.86,"deuteron gate");
   tx->SetTextColor(kBlack);
   tx->DrawLatex(0.16,0.82,"bands: p (low dE/dx) -> d -> t (high dE/dx)");
   c->SaveAs(out);
   printf("wrote %s\n",out.Data());
}
