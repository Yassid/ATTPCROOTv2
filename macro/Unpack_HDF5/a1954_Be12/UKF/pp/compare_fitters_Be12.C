/// @file compare_fitters_Be12.C
/// @brief UKF vs GENFIT on the doubly-gated 12Be proton sample. KE-vs-theta (2D + profile),
///        chi2/ndf, and the per-track UKF/GENFIT KE ratio vs theta (entry+trackID matched,
///        both fitters ran the SAME gated input). Reads <run>_ukf.root / <run>_genfit.root.
///
///   root -b -q 'pp/compare_fitters_Be12.C("run_0142,run_0143,...","/home/yassid/a1954_Be12_fit/")'
void compare_fitters_Be12(TString runsCSV, TString fitDir = "/home/yassid/a1954_Be12_fit/")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TString plotDir = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1954_Be12/UKF/pp/plots/";
   gSystem->mkdir(plotDir.Data(), kTRUE);

   const double THMAX = 110.0, KEMAX = 30.0; // data live here; caps kill unphysical outliers
   TH2F *hU = new TH2F("hU", "UKF  KE vs #theta_{lab};#theta_{lab} [deg];KE [MeV]", 220, 0, THMAX, 240, 0, KEMAX);
   TH2F *hG = new TH2F("hG", "GENFIT  KE vs #theta_{lab};#theta_{lab} [deg];KE [MeV]", 220, 0, THMAX, 240, 0, KEMAX);
   TProfile *pU = new TProfile("pU", "KE vs #theta profile;#theta_{lab} [deg];<KE> [MeV]", 55, 0, THMAX, 0, KEMAX);
   TProfile *pG = new TProfile("pG", "", 55, 0, THMAX, 0, KEMAX);
   TH1D *cU = new TH1D("cU", "#chi^{2}/ndf;#chi^{2}/ndf;tracks", 100, 0, 20);
   TH1D *cG = new TH1D("cG", "", 100, 0, 20);
   TProfile *pRatio = new TProfile("pRatio", "UKF/GENFIT KE ratio vs #theta;#theta_{lab} [deg];KE_{UKF}/KE_{GENFIT}", 22, 0, THMAX, 0.0, 3.0);

   auto getMap = [](TTree *t, TClonesArray *te, Long64_t i) {
      std::map<int, std::pair<double,double>> m; // trackID -> (KE, theta)
      t->GetEntry(i); if (te->GetEntries()==0) return m;
      auto *ev=(AtTrackingEvent*)te->At(0); if(!ev) return m;
      for (auto &ft : ev->GetFittedTracks()){ if(!ft)continue; auto&k=ft->GetKinematics();
         if(k.kineticEnergy<=0||k.kineticEnergy>1000)continue;
         m[ft->GetTrackID()]={k.kineticEnergy,k.theta*TMath::RadToDeg()}; }
      return m;
   };

   TObjArray *runs = runsCSV.Tokenize(",");
   long nU=0,nG=0,nMatch=0;
   for (int ri=0; ri<runs->GetEntries(); ++ri){
      TString run=((TObjString*)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString uf=fitDir+run+"_ukf.root", gf=fitDir+run+"_genfit.root";
      if(gSystem->AccessPathName(uf)||gSystem->AccessPathName(gf))continue;
      TFile *fU=TFile::Open(uf); TTree*tU=(TTree*)fU->Get("cbmsim"); TClonesArray*eU=nullptr; tU->SetBranchAddress("AtTrackingEvent",&eU);
      TFile *fG=TFile::Open(gf); TTree*tG=(TTree*)fG->Get("cbmsim"); TClonesArray*eG=nullptr; tG->SetBranchAddress("AtTrackingEvent",&eG);
      Long64_t N=std::min(tU->GetEntries(),tG->GetEntries());
      for(Long64_t i=0;i<N;i++){
         auto mU=getMap(tU,eU,i); auto mG=getMap(tG,eG,i);
         for(auto&kv:mU){ hU->Fill(kv.second.second,kv.second.first); pU->Fill(kv.second.second,kv.second.first); ++nU; }
         for(auto&kv:mG){ hG->Fill(kv.second.second,kv.second.first); pG->Fill(kv.second.second,kv.second.first); ++nG; }
         // chi2 (reopen fitted tracks for metadata)
         tU->GetEntry(i); if(eU->GetEntries()){auto*ev=(AtTrackingEvent*)eU->At(0); for(auto&ft:ev->GetFittedTracks()){if(!ft)continue; double nd=ft->GetTrackMetadata()->GetNdf(),c=ft->GetTrackMetadata()->GetChi2(); if(nd>0)cU->Fill(c/nd);}}
         tG->GetEntry(i); if(eG->GetEntries()){auto*ev=(AtTrackingEvent*)eG->At(0); for(auto&ft:ev->GetFittedTracks()){if(!ft)continue; double nd=ft->GetTrackMetadata()->GetNdf(),c=ft->GetTrackMetadata()->GetChi2(); if(nd>0)cG->Fill(c/nd);}}
         for(auto&kv:mU){ auto it=mG.find(kv.first); if(it!=mG.end()&&it->second.first>0){ pRatio->Fill(kv.second.second, kv.second.first/it->second.first); ++nMatch; } }
      }
      fU->Close(); fG->Close();
   }
   printf("UKF tracks=%ld  GENFIT tracks=%ld  matched=%ld\n",nU,nG,nMatch);

   TCanvas *c=new TCanvas("c","cmp",1600,1000); c->Divide(2,2);
   c->cd(1); hU->Draw("colz");
   c->cd(2); hG->Draw("colz");
   c->cd(3); pU->SetLineColor(kBlue); pU->SetLineWidth(2); pG->SetLineColor(kRed); pG->SetLineWidth(2);
   pU->SetTitle("<KE> vs #theta   (blue=UKF, red=GENFIT)"); pU->SetMinimum(0); pU->SetMaximum(KEMAX); pU->Draw(); pG->Draw("same");
   { auto*l=new TLegend(0.55,0.72,0.88,0.88); l->AddEntry(pU,"UKF","l"); l->AddEntry(pG,"GENFIT","l"); l->Draw(); }
   c->cd(4); pRatio->SetLineColor(kBlack); pRatio->SetLineWidth(2); pRatio->SetMinimum(0.8); pRatio->SetMaximum(1.8); pRatio->Draw();
   { auto*l=new TLine(0,1,THMAX,1); l->SetLineStyle(2); l->SetLineColor(kGray+1); l->Draw(); }
   TString png=plotDir+"compare_fitters_Be12.png"; c->SaveAs(png); printf("saved %s\n",png.Data());

   TCanvas *c2=new TCanvas("c2","chi2",700,600); c2->SetLogy();
   cU->SetLineColor(kBlue); cU->SetLineWidth(2); cG->SetLineColor(kRed); cG->SetLineWidth(2);
   cU->SetTitle("#chi^{2}/ndf   (blue=UKF, red=GENFIT)"); cU->Draw(); cG->Draw("same");
   { auto*l=new TLegend(0.6,0.72,0.88,0.88); l->AddEntry(cU,Form("UKF (med~%.2f)",cU->GetMean()),"l"); l->AddEntry(cG,Form("GENFIT (med~%.2f)",cG->GetMean()),"l"); l->Draw(); }
   TString png2=plotDir+"compare_fitters_chi2_Be12.png"; c2->SaveAs(png2); printf("saved %s\n",png2.Data());
}
