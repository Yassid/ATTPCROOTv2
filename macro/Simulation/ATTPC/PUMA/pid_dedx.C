/// @file pid_dedx.C
/// @brief Fix K/pi mass discrimination with a Bethe-Bloch dE/dx vs rigidity PID.
///
/// The multi-hypothesis Kalman fit picks K vs pi by reduced-chi2, which is ~coin-flip
/// (~55%): at p~375 MeV/c the K and pi helices are nearly identical over PUMA's short
/// arc, so the kinematic fit carries almost no mass information. dE/dx does: at fixed
/// RIGIDITY, K (heavier, betagamma~0.76, on the Bethe-Bloch rise) ionizes ~2.5x more
/// than pi (betagamma~2.7, min-ionizing). So classify on (Brho, dE/dx), NOT the fit.
///
/// Momentum MUST be mass-independent here -> use rigidity p = 0.2998*B*R from the PRA
/// circle radius, never the mass-dependent fit KE. dE/dx = truncated mean of per-hit
/// dQ/dl (drops the Landau tail). Expected dE/dx per species from Bethe-Bloch with ONE
/// calibration constant k fit from the pi sample (a detector-response constant, not
/// per-species training). Classify by nearest expected dE/dx in log space.
///
/// Reads matched-momentum labeled samples (branch-8 pi, branch-10 K at p~375, made by
/// pid_samples.sh with species="both"). Reports dE/dx-PID vs chi2-PID confusion.
/// Run: root -b -q pid_dedx.C
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

// Bethe-Bloch shape (arb. units; absolute scale absorbed by calibration k). Mean dE/dx
// ~ (1/beta^2)[0.5 ln(2 me beta^2 gamma^2 Tmax / I^2) - beta^2], me=0.511 MeV, I~40 eV.
double bbShape(double p_MeV, double m_MeV)
{
   const double me = 0.510999, I = 40e-6; // I in MeV
   double bg = p_MeV / m_MeV;             // beta*gamma
   double b2 = bg * bg / (1.0 + bg * bg); // beta^2
   double g2 = 1.0 + bg * bg;             // gamma^2
   double Tmax = 2 * me * bg * bg / (1 + 2 * std::sqrt(g2) * me / m_MeV + (me / m_MeV) * (me / m_MeV));
   double lnArg = 2 * me * b2 * g2 * Tmax / (I * I);
   double val = (1.0 / b2) * (0.5 * std::log(lnArg) - b2);
   return val > 0 ? val : 1e-6;
}

// truncated-mean dE/dx: per-hit dQ/dl along the phi-ordered track, drop top 30%.
double truncDeDx(const AtTrack &tr)
{
   auto cen = tr.GetGeoCenter();
   const auto &hits = tr.GetHitArray();
   if (hits.size() < 4) return -1;
   std::vector<std::tuple<double,double,double,double,double>> s; // phi,x,y,z,q
   for (const auto &h : hits) { const auto &p = h->GetPosition();
      s.emplace_back(std::atan2(p.Y()-cen.second, p.X()-cen.first), p.X(), p.Y(), p.Z(), h->GetCharge()); }
   std::sort(s.begin(), s.end());
   std::vector<double> dedx;
   for (size_t i = 1; i < s.size(); ++i) {
      double dx=std::get<1>(s[i])-std::get<1>(s[i-1]), dy=std::get<2>(s[i])-std::get<2>(s[i-1]), dz=std::get<3>(s[i])-std::get<3>(s[i-1]);
      double dl = std::sqrt(dx*dx+dy*dy+dz*dz);
      if (dl > 0.05) dedx.push_back(std::get<4>(s[i]) / dl); // charge of hit i over step to it
   }
   if (dedx.size() < 3) return -1;
   std::sort(dedx.begin(), dedx.end());
   size_t keep = (size_t)std::ceil(0.70 * dedx.size()); // drop top 30% (Landau tail)
   double sum = 0; for (size_t i = 0; i < keep; ++i) sum += dedx[i];
   return sum / keep;
}

struct Trk { double p, dedx; int truth, chi2pid; }; // truth/chi2pid: +1=K, -1=pi, 0=unknown

void loadSample(TString digiFile, TString simFile, double B, std::vector<Trk> &out)
{
   TFile fD(digiFile); TTree *tD = (TTree*)fD.Get("cbmsim");
   TFile fS(simFile);  TTree *tS = (TTree*)fS.Get("cbmsim");
   if (!tD || !tS) { printf("  MISSING %s or %s\n", digiFile.Data(), simFile.Data()); return; }
   TClonesArray *pat = new TClonesArray("AtPatternEvent"); tD->SetBranchAddress("AtPatternEvent", &pat);
   TClonesArray *ukf = new TClonesArray("AtTrackingEvent"); tD->SetBranchAddress("AtTrackingEventUKF", &ukf);
   TClonesArray *mcPts = new TClonesArray("AtMCPoint"); tS->SetBranchAddress("AtTpcPoint", &mcPts);
   TClonesArray *mcTrks = new TClonesArray("AtMCTrack"); tS->SetBranchAddress("MCTrack", &mcTrks);
   const double kTol2 = 9.0; // 3mm xy match
   Long64_t n = std::min(tD->GetEntries(), tS->GetEntries());
   for (Long64_t e = 0; e < n; ++e) {
      tD->GetEntry(e); tS->GetEntry(e);
      if (!pat->GetEntries()) continue;
      int nMC = mcPts->GetEntries();
      std::vector<double> mcX(nMC), mcY(nMC); std::vector<int> mcPdg(nMC);
      for (int k = 0; k < nMC; ++k) { auto *mp=(AtMCPoint*)mcPts->At(k); mcX[k]=mp->GetX()*10; mcY[k]=mp->GetY()*10;
         int t=mp->GetTrackID(); auto *mt=(t>=0&&t<mcTrks->GetEntries())?(AtMCTrack*)mcTrks->At(t):nullptr; mcPdg[k]=mt?mt->GetPdgCode():0; }
      // chi2-PID: map fitted-track (by TrackID) -> chosen species
      std::map<int,int> chi2pid;
      if (ukf->GetEntries()) for (const auto &ft : ((AtTrackingEvent*)ukf->At(0))->GetFittedTracks()) {
         TString id = ft->GetParticleInfo(0).idPDG;
         int sp = id.Contains("K") ? +1 : (id.Contains("pi") ? -1 : 0);
         chi2pid[ft->GetTrackID()] = sp;
      }
      int tid = 0;
      for (auto &tr : ((AtPatternEvent*)pat->At(0))->GetTrackCand()) {
         double R = tr.GetGeoRadius();
         double p = 0.299792458 * B * R;  // MeV/c, mass-independent rigidity (|q|=1, R in mm)
         double dedx = truncDeDx(tr);
         if (!(R > 0 && R < 1e5) || !(dedx > 0) || !(p > 50 && p < 1500)) { tid++; continue; }
         // truth species: majority PDG of matched MC points
         std::map<int,int> votes;
         for (const auto &h : tr.GetHitArray()) { const auto &p3=h->GetPosition(); double best=kTol2; int bp=0;
            for (int k=0;k<nMC;++k){ double d2=(p3.X()-mcX[k])*(p3.X()-mcX[k])+(p3.Y()-mcY[k])*(p3.Y()-mcY[k]); if(d2<best){best=d2;bp=mcPdg[k];} }
            if (bp) votes[bp]++; }
         int tp=0,bv=0; for (auto &kv:votes) if (kv.second>bv){bv=kv.second;tp=kv.first;}
         int truth = (std::abs(tp)==321)?+1 : (std::abs(tp)==211?-1:0);
         Trk t; t.p=p; t.dedx=dedx; t.truth=truth;
         t.chi2pid = chi2pid.count(tid)?chi2pid[tid]:0;
         out.push_back(t); tid++;
      }
   }
}

void pid_dedx(double B = 4.0)
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const double mPi = 139.57039, mK = 493.677;
   TString dir = "/mnt/f/puma_sweep/";

   std::vector<Trk> all;
   printf("loading pi sample...\n");  loadSample(dir+"output_digi_pi_pid.root", "./data/attpcsim_pi.root", B, all);
   size_t nPi = all.size();
   printf("loading K sample...\n");   loadSample(dir+"output_digi_K_pid.root",  "./data/attpcsim_K.root",  B, all);
   printf("tracks: %zu total (%zu from pi file, %zu from K file)\n", all.size(), nPi, all.size()-nPi);

   // calibrate k from truth-pi tracks: k = median(dedx_pi) / bbShape(median p, mPi)
   std::vector<double> pdPi, ppPi;
   for (auto &t : all) if (t.truth==-1) { pdPi.push_back(t.dedx); ppPi.push_back(t.p); }
   double kCal = med(pdPi) / bbShape(med(ppPi), mPi);
   printf("calibration: median pi dE/dx=%.3g at p=%.0f -> k=%.4g\n", med(pdPi), med(ppPi), kCal);

   // classify every track by nearest expected dE/dx (log space)
   auto classify = [&](const Trk &t){
      double ePi = kCal*bbShape(t.p, mPi), eK = kCal*bbShape(t.p, mK);
      return (std::abs(std::log(t.dedx)-std::log(eK)) < std::abs(std::log(t.dedx)-std::log(ePi))) ? +1 : -1;
   };

   // confusion for dE/dx-PID and chi2-PID (restrict to truth-labeled tracks)
   int cD[2][2]={{0,0},{0,0}}, cC[2][2]={{0,0},{0,0}}; // [truthIsK][pidIsK]
   int nChi=0;
   for (auto &t : all) {
      if (t.truth==0) continue;
      int tr = (t.truth==+1)?1:0;
      cD[tr][ classify(t)==+1 ? 1:0 ]++;
      if (t.chi2pid!=0){ cC[tr][ t.chi2pid==+1?1:0 ]++; nChi++; }
   }
   auto acc=[&](int c[2][2]){ int corr=c[0][0]+c[1][1], tot=c[0][0]+c[0][1]+c[1][0]+c[1][1]; return tot?100.0*corr/tot:0; };
   auto eff=[&](int c[2][2],int k){ int tot=c[k][0]+c[k][1]; return tot?100.0*c[k][k]/tot:0; };

   printf("\n================ K/pi PID confusion (matched p~375 MeV/c) ================\n");
   printf(" dE/dx-PID (Bethe-Bloch on rigidity):  overall %.1f%%\n", acc(cD));
   printf("   true pi -> pi %.1f%%  | K %.1f%%      (n=%d)\n", eff(cD,0), 100-eff(cD,0), cD[0][0]+cD[0][1]);
   printf("   true K  -> K  %.1f%%  | pi %.1f%%      (n=%d)\n", eff(cD,1), 100-eff(cD,1), cD[1][0]+cD[1][1]);
   printf(" chi2-PID (multi-hyp Kalman baseline): overall %.1f%%  (n=%d labeled+fitted)\n", acc(cC), nChi);
   printf("   true pi -> pi %.1f%%  | K %.1f%%      (n=%d)\n", eff(cC,0), 100-eff(cC,0), cC[0][0]+cC[0][1]);
   printf("   true K  -> K  %.1f%%  | pi %.1f%%      (n=%d)\n", eff(cC,1), 100-eff(cC,1), cC[1][0]+cC[1][1]);
   printf("==========================================================================\n\n");

   // figure: (p, dE/dx) scatter colored by truth + BB expected curves
   auto *c = new TCanvas("pid","dedx PID",900,650); gPad->SetGrid();
   auto *gPi=new TGraph(), *gK=new TGraph();
   gPi->SetMarkerStyle(20); gPi->SetMarkerColor(kGreen+2); gPi->SetMarkerSize(0.5);
   gK->SetMarkerStyle(20);  gK->SetMarkerColor(kRed+1);    gK->SetMarkerSize(0.5);
   std::vector<double> dvals;
   for (auto &t:all){ if(t.truth==-1) gPi->SetPoint(gPi->GetN(),t.p,t.dedx);
      else if(t.truth==+1) gK->SetPoint(gK->GetN(),t.p,t.dedx); if(t.truth) dvals.push_back(t.dedx); }
   std::sort(dvals.begin(),dvals.end());
   double dcap = dvals.empty()?4000 : 1.15*dvals[(size_t)(0.98*dvals.size())]; // 98th pct headroom
   auto *frame=new TH2F("fr",";rigidity p=0.2998 B R [MeV/c];dE/dx (trunc. mean) [a.u./mm]",10,150,600,10,0,dcap);
   frame->Draw();
   auto *bPi=new TGraph(), *bK=new TGraph(), *bBnd=new TGraph();
   for (int i=0;i<=50;++i){ double p=200+i*7.0; bPi->SetPoint(i,p,kCal*bbShape(p,mPi)); bK->SetPoint(i,p,kCal*bbShape(p,mK));
      bBnd->SetPoint(i,p,std::sqrt(kCal*bbShape(p,mPi)*kCal*bbShape(p,mK))); }
   gPi->Draw("P same"); gK->Draw("P same");
   bPi->SetLineColor(kGreen+2); bPi->SetLineWidth(3); bPi->Draw("L same");
   bK->SetLineColor(kRed+1);    bK->SetLineWidth(3);  bK->Draw("L same");
   bBnd->SetLineColor(kBlack); bBnd->SetLineStyle(2); bBnd->SetLineWidth(2); bBnd->Draw("L same");
   auto *leg=new TLegend(0.55,0.68,0.88,0.88);
   leg->AddEntry(gPi,"true #pi","p"); leg->AddEntry(gK,"true K","p");
   leg->AddEntry(bPi,"BB expected #pi","l"); leg->AddEntry(bK,"BB expected K","l");
   leg->AddEntry(bBnd,"decision boundary","l"); leg->Draw();
   c->SaveAs("./data/pid_dedx.png");
   printf("wrote ./data/pid_dedx.png\n");
}
