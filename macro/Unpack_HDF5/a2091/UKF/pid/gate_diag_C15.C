/// @file gate_diag_C15.C
/// @brief Diagnostics to SET the gates from data: avg IC trace, IC amplitude spectrum
///        (max vs -min), and the raw PID plane (sqrt(dEdx) vs Brho) with the current
///        proton polygon overlaid. Local slim + FRIB cache.
///
///   root -b -q 'pid/gate_diag_C15.C("run_0138,run_0058,run_0060,run_0148,run_0150")'
static TCutG *LoadCut(const char *path, const char *name)
{
   std::ifstream in(path); if (!in) return nullptr;
   std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   auto pos = s.find('[', s.find("vertices")); if (pos == std::string::npos) return nullptr;
   std::vector<double> nums; const char *p = s.c_str() + pos, *end = s.c_str() + s.size(); int depth = 0;
   while (p < end) { if (*p=='[') depth++; if (*p==']'){depth--; if(depth<=0){++p;break;}}
      char *np=nullptr; double v=strtod(p,&np); if(np!=p){nums.push_back(v);p=np;} else ++p; }
   auto *cut = new TCutG(name, nums.size()/2);
   for (size_t i=0,k=0;i+1<nums.size();i+=2,++k) cut->SetPoint(k,nums[i],nums[i+1]);
   cut->SetLineColor(kRed); cut->SetLineWidth(2); return cut;
}

void gate_diag_C15(TString runsCSV = "run_0138,run_0058,run_0060,run_0148,run_0150",
                    TString inDir = "/home/yassid/a2091_C15_reco_slim/", Int_t icTbLo = 1050, Int_t icTbHi = 1250,
                    double bField = 2.85,
                    TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_HDF5/a2091/"
                                         "UKF/pid/proton_band_C15.json")
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TString dir = getenv("VMCWORKDIR");
   TString plotDir = dir + "/macro/Unpack_HDF5/a2091/UKF/pid/plots/";
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");

   TProfile *hTr = new TProfile("hTr", "avg IC trace (trace[0]);time bucket;ADC", 2048, 0, 2048);
   TH1D *hMax = new TH1D("hMax", "IC amp = MAX in window;amplitude;events", 500, -200, 2500);
   TH1D *hMin = new TH1D("hMin", "IC amp = -MIN in window (pulse depth);amplitude;events", 500, -200, 2500);
   TH2F *hPID = new TH2F("hPID", "raw PID all tracks (red=current proton gate);#sqrt{dEdx};B#rho [T m]", 400, 0, 45, 400, 0, 1.2);

   TObjArray *runs = runsCSV.Tokenize(",");
   for (int ri=0; ri<runs->GetEntries(); ++ri) {
      TString run = ((TObjString*)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff=inDir+run+"_FRIB.root", rf=inDir+run+"_slim.root";
      if (gSystem->AccessPathName(ff)||gSystem->AccessPathName(rf)) continue;
      TFile *fF=TFile::Open(ff); TTree *tF=(TTree*)fF->Get("cbmsim");
      TClonesArray*ra=nullptr; tF->SetBranchAddress("AtRawEvent",&ra);
      for (Long64_t i=0;i<tF->GetEntries();i++){ tF->GetEntry(i); if(ra->GetEntries()==0)continue;
         auto*r=(AtRawEvent*)ra->At(0); if(!r||r->GetGenTraces().empty())continue;
         auto&adc=r->GetGenTraces()[0]->GetADC();
         for(size_t b=0;b<adc.size();b++) hTr->Fill(b,adc[b]);
         double mx=-1e9,mn=1e9; for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();b++){mx=std::max(mx,(double)adc[b]);mn=std::min(mn,(double)adc[b]);}
         hMax->Fill(mx); hMin->Fill(-mn); }
      fF->Close();
      TFile *fR=TFile::Open(rf); TTree *tR=(TTree*)fR->Get("cbmsim");
      TClonesArray*pe=nullptr; tR->SetBranchAddress("AtPatternEvent",&pe);
      for (Long64_t i=0;i<tR->GetEntries();i++){ tR->GetEntry(i); if(pe->GetEntries()==0)continue;
         auto*p=(AtPatternEvent*)pe->At(0); if(!p)continue;
         for(auto&trk:p->GetTrackCand()){AtTrack&tr=const_cast<AtTrack&>(trk); auto r=spy.Estimate(tr); if(!r.valid)continue; hPID->Fill(r.sqrtdEdx,r.brho);} }
      fR->Close();
   }
   TCanvas *c=new TCanvas("c","diag",1600,1000); c->Divide(2,2);
   c->cd(1); hTr->GetXaxis()->SetRangeUser(1000,1300); hTr->Draw();
   c->cd(2); c->cd(2)->SetLogy(); hMin->SetLineColor(kBlue); hMin->Draw(); hMax->SetLineColor(kRed); hMax->Draw("same");
   TLegend*L=new TLegend(0.5,0.75,0.88,0.88); L->AddEntry(hMin,"-min (pulse depth)","l"); L->AddEntry(hMax,"max (overshoot)","l"); L->Draw();
   c->cd(3); c->cd(3)->SetLogz(); hPID->Draw("colz"); if(pcut)pcut->Draw("L same");
   c->cd(4); c->cd(4)->SetLogy(); hMin->Draw();
   TString png=plotDir+"gate_diag_C15.png"; c->SaveAs(png); printf("saved %s\n",png.Data());
}
