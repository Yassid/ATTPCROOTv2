/// @file ic_subgate_C14.C
/// @brief (2) sub-divide the ~650 IC peak + (3) backward cut theta>90. Tests whether the
///        wrong-slope contamination separates by IC amplitude (a second beam species under
///        the ~650 peak). All proton-gated, backward (theta_lab>thMin) unless noted.
///
/// Panels: (1) fine IC spectrum [400,1100] w/ split lines
///         (2) IC(max) vs Brho, backward proton tracks (correlation -> two species?)
///         (3) Brho vs theta (theta>thMin), IC in [icLo,icSplit]  (LOW sub-band)
///         (4) Brho vs theta (theta>thMin), IC in [icSplit,icHi]  (HIGH sub-band)
///
///   root -b -q 'pid/ic_subgate_C14.C("run_0055,run_0058,run_0060,run_0061,run_0062,run_0063,run_0064,run_0065")'
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
   return cut;
}

void ic_subgate_C14(TString runsCSV, TString inDir = "/home/yassid/a1954_C14_reco_hdb_slim/", double icLo = 500,
                     double icHi = 900, double icSplit = 690, double thMin = 90, Int_t icTbLo = 1050,
                     Int_t icTbHi = 1250, double bField = 2.85,
                     TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1954/"
                                          "UKF/pid/proton_14C.json")
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TString plotDir = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");

   TH1D *hIC = new TH1D("hIC", "IC amplitude zoom (lines: 500/split/900);ic amplitude;events", 350, 400, 1100);
   TH2F *hICb = new TH2F("hICb", "IC(max) vs B#rho  (backward protons);ic amplitude;B#rho [T m]", 200, 400, 1100, 200, 0, 1.2);
   TH2F *bLo = new TH2F("bLo", Form("B#rho vs #theta  IC[%.0f,%.0f]  #theta>%.0f;#theta_{lab} [deg];B#rho [T m]", icLo, icSplit, thMin), 130, 80, 180, 300, 0, 1.5);
   TH2F *bHi = new TH2F("bHi", Form("B#rho vs #theta  IC[%.0f,%.0f]  #theta>%.0f;#theta_{lab} [deg];B#rho [T m]", icSplit, icHi, thMin), 130, 80, 180, 300, 0, 1.5);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nLo=0,nHi=0;
   for (int ri=0; ri<runs->GetEntries(); ++ri){
      TString run=((TObjString*)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff=inDir+run+"_FRIB.root", rf=inDir+run+"_slim.root";
      if(gSystem->AccessPathName(ff)||gSystem->AccessPathName(rf))continue;
      std::map<int,double> icByID;
      TFile*fF=TFile::Open(ff); TTree*tF=fF?(TTree*)fF->Get("cbmsim"):nullptr;
      if(!tF||tF->GetEntries()==0){if(fF)fF->Close();continue;}
      TClonesArray*ra=nullptr; tF->SetBranchAddress("AtRawEvent",&ra);
      for(Long64_t i=0;i<tF->GetEntries();i++){ tF->GetEntry(i); if(ra->GetEntries()==0)continue;
         auto*r=(AtRawEvent*)ra->At(0); if(!r||r->GetGenTraces().empty())continue;
         auto&adc=r->GetGenTraces()[0]->GetADC(); double mx=-1e9;
         for(int b=icTbLo;b<icTbHi&&b<(int)adc.size();b++)mx=std::max(mx,(double)adc[b]);
         icByID[(int)i]=mx; hIC->Fill(mx); }
      fF->Close();
      TFile*fR=TFile::Open(rf); TTree*tR=fR?(TTree*)fR->Get("cbmsim"):nullptr;
      if(!tR||tR->GetEntries()==0){if(fR)fR->Close();continue;}
      TClonesArray*pe=nullptr; tR->SetBranchAddress("AtPatternEvent",&pe);
      for(Long64_t i=0;i<tR->GetEntries();i++){ auto it=icByID.find((int)i); if(it==icByID.end())continue;
         double ic=it->second; if(ic<icLo||ic>icHi)continue;
         tR->GetEntry(i); if(pe->GetEntries()==0)continue; auto*p=(AtPatternEvent*)pe->At(0); if(!p)continue;
         for(auto&trk:p->GetTrackCand()){AtTrack&tr=const_cast<AtTrack&>(trk); auto r=spy.Estimate(tr);
            if(!r.valid||!pcut->IsInside(r.sqrtdEdx,r.brho))continue;
            double th=r.polar*TMath::RadToDeg();
            if(th>thMin) hICb->Fill(ic,r.brho);
            if(th<=thMin)continue;                         // backward cut (item 3)
            if(ic<icSplit){bLo->Fill(th,r.brho);++nLo;} else {bHi->Fill(th,r.brho);++nHi;}
         } }
      fR->Close();
   }
   printf("backward proton tracks: IC[%.0f,%.0f]=%ld  IC[%.0f,%.0f]=%ld\n",icLo,icSplit,nLo,icSplit,icHi,nHi);

   TCanvas*c=new TCanvas("c","icsub",1600,1000); c->Divide(2,2);
   c->cd(1); c->cd(1)->SetLogy(); hIC->Draw();
   for(double x:{icLo,icSplit,icHi}){auto*l=new TLine(x,0.5,x,hIC->GetMaximum()); l->SetLineColor(kRed); l->SetLineWidth(x==icSplit?1:2); if(x==icSplit)l->SetLineStyle(2); l->Draw();}
   c->cd(2); c->cd(2)->SetLogz(); hICb->Draw("colz"); { auto*l=new TLine(icSplit,0,icSplit,1.2); l->SetLineColor(kRed);l->SetLineStyle(2);l->Draw(); }
   c->cd(3); bLo->Draw("colz");
   c->cd(4); bHi->Draw("colz");
   TString png=plotDir+"ic_subgate_C14.png"; c->SaveAs(png); printf("saved %s\n",png.Data());
}
