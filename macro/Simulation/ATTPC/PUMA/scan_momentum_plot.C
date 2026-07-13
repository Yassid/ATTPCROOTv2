/// @file scan_momentum_plot.C
/// @brief Plot PUMA pion resolution vs momentum from the scan fit files
///        (max PSA + ring clustering). Momentum resolution uses the analytic
///        truth p0 = sqrt(E^2 - m^2) from the beam energy (no sim file needed),
///        so it is robust to the per-point sim being overwritten.
/// Run: root -b -q 'scan_momentum_plot.C'
double medOf(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }
double iqrS(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }

// returns {pMedFracBias%, pSigIQR%, thetaSig_deg, nfit} for one fitter branch
std::vector<double> extract(TString file, TString branch, double p0MeV)
{
   const double m=139.57039;
   TFile f(file); auto*t=(TTree*)f.Get("cbmsim"); if(!t) return {0,0,0,0};
   auto*arr=new TClonesArray("AtTrackingEvent"); t->SetBranchAddress(branch,&arr);
   std::vector<double> dp, th;
   for(Long64_t e=0;e<t->GetEntries();e++){ t->GetEntry(e); if(arr->GetEntries()==0)continue;
      auto*te=(AtTrackingEvent*)arr->At(0);
      for(const auto&ft: te->GetFittedTracks()){ const auto&kin=ft->GetKinematics(0); double KE=kin.kineticEnergy;
         if(!(KE>0))continue; double p=std::sqrt(KE*KE+2*KE*m);
         dp.push_back(100.*(p-p0MeV)/p0MeV); th.push_back(kin.theta*180./M_PI-90.); }
   }
   return {medOf(dp), iqrS(dp), iqrS(th), (double)dp.size()};
}

void scan_momentum_plot(TString prefix="data/scan_p", TString tag="", TString outdir="/Users/quantumlab/fair_install/puma_slides/figs")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62,"xyz"); gStyle->SetTitleFont(62,"xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   const double m=0.13957;
   // p[MeV], E[GeV], file
   struct Pt{double p; TString f;};
   std::vector<Pt> pts={{150,prefix+"0.150.root"},{250,prefix+"0.250.root"},{375,prefix+"0.375.root"},
                        {500,prefix+"0.500.root"},{700,prefix+"0.700.root"},{900,prefix+"0.900.root"}};
   auto gU=new TGraph, gG=new TGraph, gUt=new TGraph, gGt=new TGraph, gUb=new TGraph;
   printf("\n  p[MeV]   UKF: bias%%  sig%%  th   |  GENFIT: bias%%  sig%%  th   (nfit)\n");
   for(size_t i=0;i<pts.size();i++){ double p0=pts[i].p;
      auto u=extract(pts[i].f,"AtTrackingEventUKF",p0);
      auto g=extract(pts[i].f,"AtTrackingEventGenfit",p0);
      printf("  %5.0f    UKF %+6.1f %5.1f %4.1f   | GF %+6.1f %5.1f %4.1f   (%.0f)\n",
             p0,u[0],u[1],u[2],g[0],g[1],g[2],u[3]);
      gU->SetPoint(i,p0,u[1]); gG->SetPoint(i,p0,g[1]);
      gUt->SetPoint(i,p0,u[2]); gGt->SetPoint(i,p0,g[2]); gUb->SetPoint(i,p0,u[0]);
   }
   auto*c=new TCanvas("c","",1150,500); c->Divide(2,1);
   c->cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGrid();
   gU->SetTitle("Pion momentum resolution vs |p|;truth |p| [MeV/c];#sigma_{IQR}(p)/p  [%]");
   gU->SetMarkerStyle(20); gU->SetMarkerSize(1.4); gU->SetMarkerColor(kAzure+2); gU->SetLineColor(kAzure+2); gU->SetLineWidth(3);
   gU->GetYaxis()->SetRangeUser(0,60); gU->GetXaxis()->SetLimits(100,950); gU->Draw("ALP");
   gG->SetMarkerStyle(24); gG->SetMarkerSize(1.4); gG->SetMarkerColor(kRed+1); gG->SetLineColor(kRed+1); gG->SetLineWidth(3); gG->SetLineStyle(2); gG->Draw("LP same");
   { auto*lg=new TLegend(0.18,0.74,0.5,0.89); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
     lg->AddEntry(gU,"UKF","lp"); lg->AddEntry(gG,"GENFIT","lp"); lg->Draw(); }
   c->cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGrid();
   gUt->SetTitle("Pion polar-angle resolution vs |p|;truth |p| [MeV/c];#sigma_{IQR}(#theta)  [deg]");
   gUt->SetMarkerStyle(20); gUt->SetMarkerSize(1.4); gUt->SetMarkerColor(kAzure+2); gUt->SetLineColor(kAzure+2); gUt->SetLineWidth(3);
   gUt->GetYaxis()->SetRangeUser(0,12); gUt->GetXaxis()->SetLimits(100,950); gUt->Draw("ALP");
   gGt->SetMarkerStyle(24); gGt->SetMarkerSize(1.4); gGt->SetMarkerColor(kRed+1); gGt->SetLineColor(kRed+1); gGt->SetLineWidth(3); gGt->SetLineStyle(2); gGt->Draw("LP same");
   { auto*lg=new TLegend(0.18,0.74,0.5,0.89); lg->SetTextFont(62); lg->SetBorderSize(0); lg->SetFillStyle(0);
     lg->AddEntry(gUt,"UKF","lp"); lg->AddEntry(gGt,"GENFIT","lp"); lg->Draw(); }
   c->SaveAs(outdir+"/momentum_scan"+tag+".png");
   printf("SCANPLOT_DONE%s\n", tag.Data());
}
