/// @file show_gates_C15.C
/// @brief Show the two gates on the data BEFORE running the full analysis:
///        (L) IC max spectrum with the [icLo,icHi] cut; (R) IC-gated PID plane with the
///        proton polygon overlaid. Local slim + FRIB cache, a few runs for speed.
///
///   root -b -q 'pid/show_gates_C15.C("run_0138,run_0058,run_0060,run_0061,run_0064")'
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
   cut->SetLineColor(kRed); cut->SetLineWidth(3); return cut;
}

void show_gates_C15(TString runsCSV = "run_0138,run_0058,run_0060,run_0061,run_0064",
                     TString inDir = "/home/yassid/a2091_C15_reco_slim/", double icLo = 500, double icHi = 900,
                     Int_t icTbLo = 1050, Int_t icTbHi = 1250, double bField = 2.85,
                     TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_HDF5/a2091/"
                                          "UKF/pid/proton_15C.json")
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TString plotDir = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a2091/UKF/pid/plots/";
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");
   if (!pcut) printf("WARN: proton gate not loaded\n");

   TH1D *hIC = new TH1D("hIC", "IC amplitude (max, TB 1050-1250);amplitude;events", 500, 0, 2500);
   TH2F *hPID = new TH2F("hPID", "15C-gated PID (red = proton_15C gate);#sqrt{dEdx};B#rho [T m]", 300, 0, 45, 300, 0, 1.2);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nIC = 0, nAll = 0;
   for (int ri=0; ri<runs->GetEntries(); ++ri) {
      TString run = ((TObjString*)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff=inDir+run+"_FRIB.root", rf=inDir+run+"_slim.root";
      if (gSystem->AccessPathName(ff)||gSystem->AccessPathName(rf)) continue;
      std::map<int,double> icByID;
      TFile *fF=TFile::Open(ff); TTree *tF= fF?(TTree*)fF->Get("cbmsim"):nullptr;
      if(!tF||tF->GetEntries()==0){ if(fF)fF->Close(); continue; }
      TClonesArray*ra=nullptr; tF->SetBranchAddress("AtRawEvent",&ra);
      for(Long64_t i=0;i<tF->GetEntries();i++){ tF->GetEntry(i); if(ra->GetEntries()==0)continue;
         auto*r=(AtRawEvent*)ra->At(0); if(!r||r->GetGenTraces().empty())continue;
         auto&adc=r->GetGenTraces()[0]->GetADC(); double mx=-1e9;
         for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();b++)mx=std::max(mx,(double)adc[b]);
         icByID[i]=mx; hIC->Fill(mx); ++nAll; }
      fF->Close();
      TFile *fR=TFile::Open(rf); TTree *tR= fR?(TTree*)fR->Get("cbmsim"):nullptr;
      if(!tR||tR->GetEntries()==0){ if(fR)fR->Close(); continue; }
      TClonesArray*pe=nullptr; tR->SetBranchAddress("AtPatternEvent",&pe);
      for(Long64_t i=0;i<tR->GetEntries();i++){ tR->GetEntry(i); if(pe->GetEntries()==0)continue;
         auto it=icByID.find((int)i); if(it==icByID.end())continue;
         if(it->second<icLo||it->second>icHi)continue; ++nIC;
         auto*p=(AtPatternEvent*)pe->At(0); if(!p)continue;
         for(auto&trk:p->GetTrackCand()){AtTrack&tr=const_cast<AtTrack&>(trk); auto r=spy.Estimate(tr); if(!r.valid)continue; hPID->Fill(r.sqrtdEdx,r.brho);} }
      fR->Close();
   }
   printf("events=%ld  15C-IC-gated=%ld (%.1f%%)\n", nAll, nIC, nAll?100.0*nIC/nAll:0);

   TCanvas *c=new TCanvas("c","gates",1500,600); c->Divide(2,1);
   c->cd(1); c->cd(1)->SetLogy(); hIC->Draw();
   TLine *l1=new TLine(icLo,0.5,icLo,hIC->GetMaximum()); TLine *l2=new TLine(icHi,0.5,icHi,hIC->GetMaximum());
   l1->SetLineColor(kRed); l2->SetLineColor(kRed); l1->SetLineWidth(2); l2->SetLineWidth(2); l1->Draw(); l2->Draw();
   c->cd(2); c->cd(2)->SetLogz(); hPID->Draw("colz"); if(pcut)pcut->Draw("L same");
   TString png=plotDir+"show_gates_C15.png"; c->SaveAs(png); printf("saved %s\n",png.Data());
}
