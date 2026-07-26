/// @file ic_multiplicity_C14.C
/// @brief IC pile-up rejection via pulse-counting. Counts pulses (peaks) in the IC trace;
///        a single beam particle = 1 pulse, pile-up = >1. Splits the proton-gated Brho-vs-theta
///        by IC multiplicity to see if the wrong-slope contamination is pile-up.
///
/// Panels: (1) #IC-peaks distribution (events whose max is in the 14C window)
///         (2) Bp-vs-theta, mult==1 AND single peak in [icLo,icHi]  (CLEAN 14C)
///         (3) Bp-vs-theta, mult>=2 with a peak in [icLo,icHi]      (PILE-UP)
///         (4) Bp-vs-theta, current gate (max in [icLo,icHi], any mult)  (reference/dirty)
///
///   root -b -q 'pid/ic_multiplicity_C14.C("run_0055,run_0058,run_0060,run_0061,run_0064")'
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

// Count pulses in adc[tbLo,tbHi]: contiguous runs above `thr` (hysteresis at thr/2) = one pulse.
// Returns the per-pulse peak heights.
static std::vector<double> FindPulses(const std::vector<Double_t> &adc, double thr, int tbLo, int tbHi)
{
   std::vector<double> pk;
   int b = tbLo, N = (int)adc.size();
   if (tbHi > N) tbHi = N;
   while (b < tbHi) {
      if (adc[b] > thr) {
         double mx = adc[b]; int b2 = b;
         while (b2 < tbHi && adc[b2] > thr * 0.5) { if (adc[b2] > mx) mx = adc[b2]; ++b2; }
         pk.push_back(mx);
         b = b2;
      } else ++b;
   }
   return pk;
}

void ic_multiplicity_C14(TString runsCSV, TString inDir = "/home/yassid/a1954_C14_reco_hdb_slim/", double icLo = 500,
                          double icHi = 900, double peakThr = 200, Int_t tbLo = 800, Int_t tbHi = 1500,
                          double bField = 2.85,
                          TString protonGate = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/"
                                               "a1954/UKF/pid/proton_14C.json")
{
   gSystem->Load("libAtTools.so"); gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetPalette(kBird); gStyle->SetNumberContours(255);
   TString plotDir = TString(getenv("VMCWORKDIR")) + "/macro/Unpack_HDF5/a1954/UKF/pid/plots/";
   AtTools::AtSpyralPID spy; spy.SetBField(bField);
   TCutG *pcut = LoadCut(protonGate.Data(), "pcut");

   TH1D *hN = new TH1D("hN", "IC pulse multiplicity (max in 14C window);# IC pulses;events", 8, 0.5, 8.5);
   TH2F *b1 = new TH2F("b1", "B#rho vs #theta  MULT==1 (clean 14C);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 300, 0, 1.5);
   TH2F *b2 = new TH2F("b2", "B#rho vs #theta  MULT#geq2 (pile-up);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 300, 0, 1.5);
   TH2F *b0 = new TH2F("b0", "B#rho vs #theta  current gate (any mult);#theta_{lab} [deg];B#rho [T m]", 180, 0, 180, 300, 0, 1.5);

   TObjArray *runs = runsCSV.Tokenize(",");
   long nMaxGate = 0, nClean = 0, nPile = 0;
   for (int ri=0; ri<runs->GetEntries(); ++ri){
      TString run=((TObjString*)runs->At(ri))->GetString().Strip(TString::kBoth);
      TString ff=inDir+run+"_FRIB.root", rf=inDir+run+"_slim.root";
      if(gSystem->AccessPathName(ff)||gSystem->AccessPathName(rf))continue;
      // per-entry: pulse heights + classification
      std::map<int,int> classOf; // 0=not-max-gated, 1=clean single, 2=pileup
      TFile*fF=TFile::Open(ff); TTree*tF=fF?(TTree*)fF->Get("cbmsim"):nullptr;
      if(!tF||tF->GetEntries()==0){if(fF)fF->Close();continue;}
      TClonesArray*ra=nullptr; tF->SetBranchAddress("AtRawEvent",&ra);
      for(Long64_t i=0;i<tF->GetEntries();i++){ tF->GetEntry(i); if(ra->GetEntries()==0)continue;
         auto*r=(AtRawEvent*)ra->At(0); if(!r||r->GetGenTraces().empty())continue;
         auto&adc=r->GetGenTraces()[0]->GetADC();
         double mx=-1e9; for(int b=1050;b<1250&&b<(int)adc.size();b++)mx=std::max(mx,(double)adc[b]);
         bool maxGate=(mx>=icLo&&mx<=icHi);
         if(!maxGate){classOf[(int)i]=0;continue;}
         auto pk=FindPulses(adc,peakThr,tbLo,tbHi);
         int nInWin=0; for(double h:pk) if(h>=icLo&&h<=icHi)++nInWin;
         hN->Fill(std::min((int)pk.size(),8));
         if(pk.size()==1) classOf[(int)i]=1;        // exactly one pulse total -> clean
         else classOf[(int)i]=2;                    // more than one pulse -> pile-up
      }
      fF->Close();
      TFile*fR=TFile::Open(rf); TTree*tR=fR?(TTree*)fR->Get("cbmsim"):nullptr;
      if(!tR||tR->GetEntries()==0){if(fR)fR->Close();continue;}
      TClonesArray*pe=nullptr; tR->SetBranchAddress("AtPatternEvent",&pe);
      for(Long64_t i=0;i<tR->GetEntries();i++){ auto ci=classOf.find((int)i); if(ci==classOf.end()||ci->second==0)continue;
         tR->GetEntry(i); if(pe->GetEntries()==0)continue; auto*p=(AtPatternEvent*)pe->At(0); if(!p)continue;
         for(auto&trk:p->GetTrackCand()){AtTrack&tr=const_cast<AtTrack&>(trk); auto r=spy.Estimate(tr);
            if(!r.valid||!pcut->IsInside(r.sqrtdEdx,r.brho))continue;
            double th=r.polar*TMath::RadToDeg();
            b0->Fill(th,r.brho); ++nMaxGate;
            if(ci->second==1){b1->Fill(th,r.brho);++nClean;}
            else {b2->Fill(th,r.brho);++nPile;}
         } }
      fR->Close();
   }
   printf("proton tracks: maxGate=%ld  clean(mult1)=%ld (%.1f%%)  pileup=%ld (%.1f%%)\n", nMaxGate, nClean,
          nMaxGate?100.0*nClean/nMaxGate:0, nPile, nMaxGate?100.0*nPile/nMaxGate:0);

   TCanvas*c=new TCanvas("c","icmult",1600,1000); c->Divide(2,2);
   c->cd(1); c->cd(1)->SetLogy(); hN->Draw();
   c->cd(2); b1->Draw("colz");
   c->cd(3); b2->Draw("colz");
   c->cd(4); b0->Draw("colz");
   TString png=plotDir+"ic_multiplicity_C14.png"; c->SaveAs(png); printf("saved %s\n",png.Data());
}
