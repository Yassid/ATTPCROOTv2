/// @file diag_shape_2p.C
/// @brief Momentum-shape 2x2: rows = beam momentum {375, 200 MeV/c}, each overlaying
///        pad-centre (grey) vs resistive sub-pad (blue) PRA rigidity dp/p. Tests the
///        pad-quantization hypothesis two ways at once: (a) sub-pad centroiding should
///        wash the comb into a Gaussian at fixed p; (b) lowering p to 200 (sagitta
///        ~2.3mm ~ 1 pad, vs 1.2mm sub-pad at 375) should already be more Gaussian even
///        for pad-centre hits. Prints skew/exKurt for all four.
/// Run: root -b -q diag_shape_2p.C
double skw(const std::vector<double>&v){ if(v.size()<3)return 0; double m=0;for(double x:v)m+=x;m/=v.size();
   double s2=0,s3=0;for(double x:v){double d=x-m;s2+=d*d;s3+=d*d*d;}s2/=v.size();s3/=v.size();return s2>0?s3/std::pow(s2,1.5):0;}
double krt(const std::vector<double>&v){ if(v.size()<4)return 0; double m=0;for(double x:v)m+=x;m/=v.size();
   double s2=0,s4=0;for(double x:v){double d=x-m;s2+=d*d;s4+=d*d*d*d;}s2/=v.size();s4/=v.size();return s2>0?s4/(s2*s2)-3:0;}
double iqrp(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }

std::vector<double> collectDP(TString file,double B,double p0){
   std::vector<double> d; TFile f(file); TTree*t=(TTree*)f.Get("cbmsim");
   if(!t){printf("  MISSING %s\n",file.Data()); return d;}
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent",&pat);
   for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e); if(!pat->GetEntries())continue;
      for(auto&tr:((AtPatternEvent*)pat->At(0))->GetTrackCand()){ double r=tr.GetGeoRadius();
         if(r>0&&r<1e5) d.push_back(100*(0.299792458*B*r-p0)/p0); } }
   return d;
}

void panel(std::vector<double>&pad,std::vector<double>&sub,const char*title){
   gPad->SetGrid();
   auto*hp=new TH1F(Form("hp%s",title),Form("%s;#Deltap/p [%%];tracks (norm)",title),120,-60,60);
   auto*hs=new TH1F(Form("hs%s",title),"",120,-60,60);
   for(double x:pad)hp->Fill(x); for(double x:sub)hs->Fill(x);
   if(hp->Integral()>0)hp->Scale(1./hp->Integral()); if(hs->Integral()>0)hs->Scale(1./hs->Integral());
   hp->SetLineColor(kGray+2); hp->SetLineWidth(2); hs->SetLineColor(kAzure+1); hs->SetLineWidth(2);
   hp->SetMaximum(1.3*std::max(hp->GetMaximum(),hs->GetMaximum()));
   hp->Draw("hist"); hs->Draw("hist same");
   TLatex L; L.SetNDC(); L.SetTextSize(0.045);
   L.SetTextColor(kGray+2); L.DrawLatex(0.14,0.85,Form("pad: skew %+.1f kurt %+.0f IQR%.0f%%",skw(pad),krt(pad),iqrp(pad)));
   L.SetTextColor(kAzure+1); L.DrawLatex(0.14,0.79,Form("sub: skew %+.1f kurt %+.0f IQR%.0f%%",skw(sub),krt(sub),iqrp(sub)));
}

void diag_shape_2p(double B=4.0){
   gSystem->Load("libAtReconstruction.so"); gStyle->SetOptStat(0);
   TString D="/mnt/f/puma_sweep/";
   auto pad375=collectDP(D+"output_digi_pi_pid.root",     B,374.9);
   auto sub375=collectDP(D+"output_digi_pi_prfring.root", B,374.9);
   auto pad200=collectDP(D+"output_digi_pi200_pad.root",  B,200.0);
   auto sub200=collectDP(D+"output_digi_pi200_prfring.root",B,200.0);
   printf("\n==== momentum-shape vs beam p and readout (skew,exKurt,IQR) ====\n");
   auto rep=[&](const char*n,std::vector<double>&v){ printf("  %-16s n=%5zu  skew=%+6.1f  exKurt=%+8.1f  IQR=%.1f%%\n",n,v.size(),skw(v),krt(v),iqrp(v)); };
   rep("375 pad-centre",pad375); rep("375 sub-pad",sub375);
   rep("200 pad-centre",pad200); rep("200 sub-pad",sub200);

   auto*c=new TCanvas("s2","shape vs p",1300,900); c->Divide(2,2);
   c->cd(1); panel(pad375,sub375,"p=375 MeV/c (sagitta~1.2mm, sub-pad)");
   c->cd(2); { auto*l=new TLegend(0.4,0.4,0.88,0.6); auto*a=new TH1F("la","",1,0,1),*b=new TH1F("lb","",1,0,1);
      a->SetLineColor(kGray+2);a->SetLineWidth(2); b->SetLineColor(kAzure+1);b->SetLineWidth(2);
      l->AddEntry(a,"pad-centre hits","l"); l->AddEntry(b,"resistive sub-pad (PRF+ring)","l");
      l->SetHeader("Gaussian => skew~0, exKurt~0"); l->Draw();
      TLatex T;T.SetNDC();T.SetTextSize(0.04);T.DrawLatex(0.1,0.85,"comb = pad-quantized radius");
      T.DrawLatex(0.1,0.30,"lower p (row 2) = bigger sagitta");
      T.DrawLatex(0.1,0.24,"= better sampled by the pads"); }
   c->cd(3); panel(pad200,sub200,"p=200 MeV/c (sagitta~2.3mm ~ 1 pad)");
   c->cd(4); { gPad->SetGrid(); // radius overlay at 200 to show quantization directly
      auto*Hr=new TH1F("HR2p",";PRA radius R [mm];tracks(norm)",140,50,500),*Hs=new TH1F("HR2s","",140,50,500);
      TFile f(D+"output_digi_pi200_pad.root"); TTree*t=(TTree*)f.Get("cbmsim");
      if(t){ TClonesArray*pt=new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent",&pt);
         for(Long64_t e=0;e<t->GetEntries();++e){t->GetEntry(e); if(!pt->GetEntries())continue;
            for(auto&tr:((AtPatternEvent*)pt->At(0))->GetTrackCand()){double r=tr.GetGeoRadius(); if(r>0&&r<1e5)Hr->Fill(r);} } }
      TFile f2(D+"output_digi_pi200_prfring.root"); TTree*t2=(TTree*)f2.Get("cbmsim");
      if(t2){ TClonesArray*pt=new TClonesArray("AtPatternEvent"); t2->SetBranchAddress("AtPatternEvent",&pt);
         for(Long64_t e=0;e<t2->GetEntries();++e){t2->GetEntry(e); if(!pt->GetEntries())continue;
            for(auto&tr:((AtPatternEvent*)pt->At(0))->GetTrackCand()){double r=tr.GetGeoRadius(); if(r>0&&r<1e5)Hs->Fill(r);} } }
      if(Hr->Integral()>0)Hr->Scale(1./Hr->Integral()); if(Hs->Integral()>0)Hs->Scale(1./Hs->Integral());
      Hr->SetLineColor(kGray+2);Hr->SetLineWidth(2);Hs->SetLineColor(kAzure+1);Hs->SetLineWidth(2);
      Hr->SetMaximum(1.3*std::max(Hr->GetMaximum(),Hs->GetMaximum()));
      Hr->SetTitle("PRA radius R at p=200"); Hr->Draw("hist"); Hs->Draw("hist same");
      TLine*ln=new TLine(166.8,0,166.8,Hr->GetMaximum()); ln->SetLineStyle(2);ln->SetLineColor(kRed);ln->Draw(); }
   c->SaveAs("./data/diag_shape_2p.png");
   printf("wrote ./data/diag_shape_2p.png\n\n");
}
