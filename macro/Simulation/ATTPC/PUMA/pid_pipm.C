/// @file pid_pipm.C
/// @brief Quantify the pi+ / pi- PID separation. pi+ and pi- have identical mass
///        (identical dE/dx), so they can ONLY be separated by CHARGE SIGN = the
///        curvature sense in the 4 T field. The PID variable is signed magnetic
///        rigidity  q*Brho = sign * p/0.2998 [T.m]:  pi+ -> +, pi- -> -.
///        Reports charge-ID accuracy, mis-ID rate, and a separation power, and
///        draws the signed-Brho distributions for truth pi+ vs pi-.
/// Run: root -b -q 'pid_pipm.C("data/output_digi_pi8_ring.root","data/attpcsim.root","UKF")'
double medOf(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }
double iqrS(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }

void pid_pipm(TString fitFile = "data/output_digi_pi8_ring.root", TString simFile = "data/attpcsim.root",
              TString fitter = "UKF", double Bfield = 4.0,
              TString outdir = "/Users/quantumlab/fair_install/puma_slides/figs", TString tag = "")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62,"xyz"); gStyle->SetTitleFont(62,"xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const double m = 139.57039; // pion mass MeV
   const double kMatchTol = 3.0;

   TFile fD(fitFile); auto *tD=(TTree*)fD.Get("cbmsim");
   TFile fS(simFile); auto *tS=(TTree*)fS.Get("cbmsim");
   if(!tD||!tS){printf("missing tree\n");return;}
   auto *teArr=new TClonesArray("AtTrackingEvent");
   auto *patArr=new TClonesArray("AtPatternEvent");
   tD->SetBranchAddress(Form("AtTrackingEvent%s",fitter.Data()),&teArr);
   tD->SetBranchAddress("AtPatternEvent",&patArr);
   auto *mcPts=new TClonesArray("AtMCPoint"); auto *mcTrks=new TClonesArray("AtMCTrack");
   tS->SetBranchAddress("AtTpcPoint",&mcPts); tS->SetBranchAddress("MCTrack",&mcTrks);

   std::vector<double> brhoPlus, brhoMinus; // signed q*Brho by TRUTH species
   int nPlus=0,nMinus=0, plusRight=0, minusRight=0;
   Long64_t nE=std::min(tD->GetEntries(),tS->GetEntries());
   for(Long64_t e=0;e<nE;e++){ tD->GetEntry(e); tS->GetEntry(e);
      if(patArr->GetEntries()==0||teArr->GetEntries()==0) continue;
      auto *pat=(AtPatternEvent*)patArr->At(0); auto trackArr=pat->GetTrackCand();
      int nMC=mcPts->GetEntries(); std::vector<double> mcX(nMC),mcY(nMC); std::vector<int> mcPdg(nMC);
      for(int k=0;k<nMC;k++){ auto*mp=(AtMCPoint*)mcPts->At(k); mcX[k]=mp->GetX()*10; mcY[k]=mp->GetY()*10;
         int tid=mp->GetTrackID(); auto*mt=(tid>=0&&tid<mcTrks->GetEntries())?(AtMCTrack*)mcTrks->At(tid):nullptr;
         mcPdg[k]=mt?mt->GetPdgCode():0; }
      auto *te=(AtTrackingEvent*)teArr->At(0);
      for(const auto &ft: te->GetFittedTracks()){
         const auto &kin=ft->GetKinematics(0); double KE=kin.kineticEnergy; if(!(KE>0)) continue;
         double p=std::sqrt(KE*KE+2*KE*m)/1000.0; // GeV/c
         double brho=p/0.2997925;                 // T.m (|q|=1)
         // reco charge sign
         const auto &pinfo=ft->GetParticleInfo(0); int fitSign=0;
         if(pinfo.charge!=0) fitSign=pinfo.charge>0?1:-1;
         else { TString id=pinfo.idPDG; if(id.Contains("+"))fitSign=1; else if(id.Contains("-"))fitSign=-1; }
         if(fitSign==0) continue;
         // truth sign via hit->MC majority
         int tid=ft->GetTrackID(); if(tid<0||tid>=(int)trackArr.size()) continue;
         std::map<int,int> votes;
         for(const auto &h: trackArr[tid].GetHitArray()){ const auto&pp=h->GetPosition(); double best=kMatchTol*kMatchTol; int bp=0;
            for(int k=0;k<nMC;k++){ double d2=(pp.X()-mcX[k])*(pp.X()-mcX[k])+(pp.Y()-mcY[k])*(pp.Y()-mcY[k]); if(d2<best){best=d2;bp=mcPdg[k];}}
            if(bp!=0)votes[bp]++; }
         int truthPdg=0,bv=0; for(auto&kv:votes) if(kv.second>bv){bv=kv.second;truthPdg=kv.first;}
         if(truthPdg==0) continue;
         int truthSign=truthPdg>0?1:-1;
         double signedBrho=fitSign*brho;
         if(truthSign>0){ nPlus++; brhoPlus.push_back(signedBrho); if(fitSign==truthSign)plusRight++; }
         else           { nMinus++; brhoMinus.push_back(signedBrho); if(fitSign==truthSign)minusRight++; }
      }
   }
   // separation power on |Brho| peaks: S = |mean+ - mean-| / (sig+ + sig-)
   double mp2=medOf(brhoPlus), mm=medOf(brhoMinus), sp=iqrS(brhoPlus), sm=iqrS(brhoMinus);
   double sep=(sp+sm>0)?std::abs(mp2-mm)/(sp+sm):0;
   printf("\n==== pi+/pi- PID  (%s, B=%.1f T) ====\n", fitter.Data(), Bfield);
   printf("  truth pi+: %d   charge-ID %.1f%%   median q*Brho=%+.3f T.m  sigma=%.3f\n",nPlus,nPlus?100.*plusRight/nPlus:0,mp2,sp);
   printf("  truth pi-: %d   charge-ID %.1f%%   median q*Brho=%+.3f T.m  sigma=%.3f\n",nMinus,nMinus?100.*minusRight/nMinus:0,mm,sm);
   printf("  overall charge-ID: %.1f%%   |  pi+/pi- SEPARATION POWER S=%.2f\n",
          (nPlus+nMinus)?100.*(plusRight+minusRight)/(nPlus+nMinus):0, sep);

   // draw
   double lim=1.2*std::max(std::abs(mp2)+3*sp,std::abs(mm)+3*sm); if(lim<=0)lim=1;
   auto *hP=new TH1F("hP","",50,-lim,lim); auto *hM=new TH1F("hM","",50,-lim,lim);
   for(double v:brhoPlus)hP->Fill(v); for(double v:brhoMinus)hM->Fill(v);
   auto *c=new TCanvas("c","",760,560); c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
   hP->SetLineColor(kRed+1); hP->SetLineWidth(3); hM->SetLineColor(kAzure+2); hM->SetLineWidth(3); hM->SetLineStyle(2);
   hP->SetTitle(Form("#pi^{+}/#pi^{-} PID via signed rigidity (%s);q_{reco}#times B#rho  [T m];tracks",fitter.Data()));
   hP->SetMaximum(1.25*std::max(hP->GetMaximum(),hM->GetMaximum())); hP->Draw("hist"); hM->Draw("hist same");
   auto *l0=new TLine(0,0,0,hP->GetMaximum()); l0->SetLineColor(kGray+2); l0->SetLineStyle(3); l0->Draw();
   auto *lg=new TLegend(0.62,0.74,0.88,0.9); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(hP,Form("truth #pi^{+} (%.0f%%)",nPlus?100.*plusRight/nPlus:0),"l");
   lg->AddEntry(hM,Form("truth #pi^{-} (%.0f%%)",nMinus?100.*minusRight/nMinus:0),"l");
   lg->Draw();
   auto *tx=new TLatex(); tx->SetNDC(); tx->SetTextFont(62); tx->SetTextSize(0.045);
   tx->DrawLatex(0.15,0.84,Form("charge-ID %.1f%%",(nPlus+nMinus)?100.*(plusRight+minusRight)/(nPlus+nMinus):0));
   c->SaveAs(outdir+"/pid_pipm"+tag+".png");
   printf("PIDPIPM_DONE\n");
}
