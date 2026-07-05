/// @file diag_shape_compare.C
/// @brief Direct overlay: does resistive sub-pad centroiding turn the pad-quantized
///        momentum COMB into a smooth Gaussian? Overlays the PRA rigidity (=0.2998*B*R)
///        and radius distributions for two digitizations of the SAME monoenergetic pi
///        sim: pad-centre hits (no charge sharing) vs charge-sharing PRF + ring
///        centroiding (sub-pad). If the comb+tail -> unimodal bump, the "19%" comb was
///        a pad-quantization artifact and the resistive readout fixes the SHAPE.
/// Run: root -b -q diag_shape_compare.C
double skew(const std::vector<double>&v){ if(v.size()<3)return 0; double m=0;for(double x:v)m+=x;m/=v.size();
   double s2=0,s3=0;for(double x:v){double d=x-m;s2+=d*d;s3+=d*d*d;}s2/=v.size();s3/=v.size();return s2>0?s3/std::pow(s2,1.5):0;}
double kurt(const std::vector<double>&v){ if(v.size()<4)return 0; double m=0;for(double x:v)m+=x;m/=v.size();
   double s2=0,s4=0;for(double x:v){double d=x-m;s2+=d*d;s4+=d*d*d*d;}s2/=v.size();s4/=v.size();return s2>0?s4/(s2*s2)-3:0;}

void collect(TString file,double B,std::vector<double>&P,std::vector<double>&R){
   TFile f(file); TTree*t=(TTree*)f.Get("cbmsim"); if(!t){printf("missing %s\n",file.Data());return;}
   TClonesArray*pat=new TClonesArray("AtPatternEvent"); t->SetBranchAddress("AtPatternEvent",&pat);
   for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e); if(!pat->GetEntries())continue;
      for(auto&tr:((AtPatternEvent*)pat->At(0))->GetTrackCand()){ double r=tr.GetGeoRadius();
         if(r>0&&r<1e5){ R.push_back(r); P.push_back(0.299792458*B*r);} } }
}

void diag_shape_compare(TString rawFile="/mnt/f/puma_sweep/output_digi_pi_pid.root",
                        TString prfFile="/mnt/f/puma_sweep/output_digi_pi_prfring.root",
                        double B=4.0,double p0=374.9){
   gSystem->Load("libAtReconstruction.so"); gStyle->SetOptStat(0);
   std::vector<double> Praw,Rraw,Pprf,Rprf;
   collect(rawFile,B,Praw,Rraw); collect(prfFile,B,Pprf,Rprf);
   auto dp=[&](std::vector<double>&P){ std::vector<double>d; for(double x:P)d.push_back(100*(x-p0)/p0); return d; };
   auto draw=[&](std::vector<double>&P,const char*nm,int col,int fbins){
      auto*h=new TH1F(Form("h%s",nm),Form(";#Deltap/p [%%];tracks (norm)"),fbins,-60,60);
      for(double x:dp(P)) h->Fill(x); if(h->Integral()>0)h->Scale(1.0/h->Integral());
      h->SetLineColor(col); h->SetLineWidth(2); return h; };
   printf("\n==== momentum-shape: pad-centre vs resistive sub-pad (monoenergetic pi) ====\n");
   { auto d=dp(Praw); printf("  pad-centre : n=%zu skew=%+.2f exKurt=%+.1f\n",Praw.size(),skew(d),kurt(d)); }
   { auto d=dp(Pprf); printf("  sub-pad PRF: n=%zu skew=%+.2f exKurt=%+.1f\n",Pprf.size(),skew(d),kurt(d)); }

   auto*c=new TCanvas("cmp","shape compare",1300,560); c->Divide(2,1);
   c->cd(1); gPad->SetGrid();
   auto*hr=draw(Praw,"raw",kGray+2,120); auto*hp=draw(Pprf,"prf",kAzure+1,120);
   double mx=std::max(hr->GetMaximum(),hp->GetMaximum()); hr->SetMaximum(1.25*mx);
   hr->SetTitle("#Deltap/p shape: pad-centre (grey) vs sub-pad PRF (blue)");
   hr->Draw("hist"); hp->Draw("hist same");
   { auto*l=new TLegend(0.6,0.72,0.88,0.88); l->AddEntry(hr,"pad-centre (comb)","l");
     l->AddEntry(hp,"resistive sub-pad","l"); l->Draw(); }
   c->cd(2); gPad->SetGrid();
   auto*Hr=new TH1F("HRr",";PRA radius R [mm];tracks (norm)",120,100,700);
   auto*Hp=new TH1F("HRp","",120,100,700);
   for(double r:Rraw)Hr->Fill(r); for(double r:Rprf)Hp->Fill(r);
   if(Hr->Integral()>0)Hr->Scale(1./Hr->Integral()); if(Hp->Integral()>0)Hp->Scale(1./Hp->Integral());
   Hr->SetLineColor(kGray+2); Hr->SetLineWidth(2); Hp->SetLineColor(kAzure+1); Hp->SetLineWidth(2);
   Hr->SetMaximum(1.25*std::max(Hr->GetMaximum(),Hp->GetMaximum()));
   Hr->SetTitle("PRA radius R: pad-centre vs sub-pad"); Hr->Draw("hist"); Hp->Draw("hist same");
   { TLine*ln=new TLine(p0/(0.299792458*B),0,p0/(0.299792458*B),Hr->GetMaximum()); ln->SetLineStyle(2); ln->SetLineColor(kRed); ln->Draw(); }
   c->SaveAs("./data/diag_shape_compare.png");
   printf("wrote ./data/diag_shape_compare.png\n\n");
}
