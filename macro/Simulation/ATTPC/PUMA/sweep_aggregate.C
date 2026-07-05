/// @file sweep_aggregate.C
/// @brief Aggregate the PRF charge-dispersion sweep (SPSC-P-361 sec 7.4.1). For each
///        sweep point it reports sigma_p/p (IQR) AND median (bias) for three estimators
///        on the SAME digitized events:
///          UKF        — full unscented-KF pipeline (AtTrackingEventUKF)
///          genfit     — full genfit KalmanFitterRefTrack (AtTrackingEventGenfit)
///          Taubin*    — offline framework Taubin fit on per-ring charge centroids
///                       (the physics ceiling of resistive centroiding for that sigma)
///        Always prints MEDIAN so Kasa-style bias can never masquerade as resolution.
/// Run: root -b -q sweep_aggregate.C
double iqr(std::vector<double> v){ if(v.size()<4)return 0; std::sort(v.begin(),v.end()); return (v[3*v.size()/4]-v[v.size()/4])/1.349; }
double med(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(),v.end()); return v[v.size()/2]; }

struct Pt { const char *tag; double sig; int ring; const char *file; };

void sweep_aggregate()
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   const double Bfield = 4.0, p0 = 374.9, mass = 139.57039, vz0 = 75.0;

   std::vector<Pt> pts = {
      {"base ",   0.0, 0, "/mnt/f/puma_sweep/output_digi_sw_base.root"},
      {"ring0",   0.0, 1, "/mnt/f/puma_sweep/output_digi_sw_r00.root"},
      {"0.3mm",   0.3, 1, "/mnt/f/puma_sweep/output_digi_sw_r03.root"},
      {"0.5mm",   0.5, 1, "/mnt/f/puma_sweep/output_digi_sw_r05.root"},
      {"0.8mm",   0.8, 1, "/mnt/f/puma_sweep/output_digi_sw_r08.root"},
      {"1.2mm",   1.2, 1, "/mnt/f/puma_sweep/output_digi_sw_r12.root"},
      {"1.5mm",   1.5, 1, "/mnt/f/puma_sweep/output_digi_sw_r15.root"},
   };

   // ring geometry for offline Taubin ceiling
   const double kRin = 62.9, kRout = 121.1; const int kNring = 16;
   std::vector<double> redge(kNring + 1); redge[0] = kRin;
   double dA = (kRout * kRout - kRin * kRin) / kNring;
   for (int i = 0; i < kNring; ++i) redge[i + 1] = std::sqrt(redge[i] * redge[i] + dA);
   auto ringOf = [&](double r) { for (int i = 0; i < kNring; ++i) if (r >= redge[i] && r < redge[i + 1]) return i; return (r < kRin) ? -1 : kNring - 1; };
   auto taubin = [](const std::vector<std::pair<double,double>> &q) -> double {
      if (q.size() < 5) return -1;
      std::vector<AtHit> hs; hs.reserve(q.size()); std::vector<const AtHit*> hp;
      for (auto &a : q) hs.emplace_back(0, ROOT::Math::XYZPoint(a.first, a.second, 0.0), 1.0);
      for (auto &h : hs) hp.push_back(&h);
      AtPatterns::AtPatternCircle2D c; c.AtPattern::FitPattern(hp, -1.0);
      double R = c.GetRadius(); return (R > 0 && R < 1e5) ? R : -1;
   };

   struct Row { double sig; double us,um; int un; double gs,gm; int gn; double ts,tm; int tn; };
   std::vector<Row> rows;

   printf("\n============ PRF charge-dispersion sweep (branch-8 pi, B=%.1fT, |p|=%.1f MeV/c) ============\n", Bfield, p0);
   printf(" grounded in SPSC-P-361 sec 7.4.1: Telegraph sigma=sqrt(2t/RC); pad pitch ~2.3mm azim.\n");
   printf(" %-6s %-5s | %-22s | %-22s | %-22s\n", "point", "sig", "UKF (full)", "genfit (full)", "Taubin ring (ceiling)");
   printf(" %-6s %-5s | %-22s | %-22s | %-22s\n", "", "mm", "sig%%  med%%   n", "sig%%  med%%   n", "sig%%  med%%   n");
   printf(" -------------------------------------------------------------------------------------------\n");

   for (auto &P : pts) {
      TFile *f = TFile::Open(P.file);
      if (!f || f->IsZombie()) { printf(" %-6s  MISSING (%s)\n", P.tag, P.file); if(f) f->Close(); continue; }
      TTree *t = (TTree*)f->Get("cbmsim");
      if (!t) { printf(" %-6s  no tree (%s)\n", P.tag, P.file); f->Close(); continue; }
      TClonesArray *ukf = new TClonesArray("AtTrackingEvent"); t->SetBranchAddress("AtTrackingEventUKF", &ukf);
      TClonesArray *gf  = new TClonesArray("AtTrackingEvent"); t->SetBranchAddress("AtTrackingEventGenfit", &gf);
      TClonesArray *pe  = new TClonesArray("AtPatternEvent");  t->SetBranchAddress("AtPatternEvent", &pe);
      std::vector<double> du, dg, dt;
      TClonesArray *arr[2] = {ukf, gf}; std::vector<double> *dst[2] = {&du, &dg};
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         for (int fi = 0; fi < 2; ++fi) {
            if (!arr[fi]->GetEntries()) continue;
            for (const auto &ft : ((AtTrackingEvent*)arr[fi]->At(0))->GetFittedTracks()) {
               double KE = ft->GetKinematics(0).kineticEnergy; if (!(KE > 0)) continue;
               double p = std::sqrt(KE*KE + 2*KE*mass); dst[fi]->push_back((p - p0)/p0);
            }
         }
         // offline Taubin ring-centroid ceiling
         if (pe->GetEntries()) for (auto &tr : ((AtPatternEvent*)pe->At(0))->GetTrackCand()) {
            std::vector<double> sx(kNring,0), sy(kNring,0), sq(kNring,0);
            for (auto &h : tr.GetHitArray()) {
               auto pp = h->GetPosition(); double q = std::max(1e-6,(double)h->GetCharge());
               double r = std::hypot(pp.X(),pp.Y()); int ri = ringOf(r); if (ri<0) continue;
               sx[ri]+=q*pp.X(); sy[ri]+=q*pp.Y(); sq[ri]+=q;
            }
            std::vector<std::pair<double,double>> rp;
            for (int i=0;i<kNring;++i) if (sq[i]>0) rp.emplace_back(sx[i]/sq[i], sy[i]/sq[i]);
            double Rc = taubin(rp); if (Rc>0) dt.push_back((0.299792458*Bfield*Rc - p0)/p0);
         }
      }
      printf(" %-6s %-5.2f | %5.1f %+6.1f %5zu | %5.1f %+6.1f %5zu | %5.1f %+6.1f %5zu\n",
             P.tag, P.sig,
             100*iqr(du),100*med(du),du.size(),
             100*iqr(dg),100*med(dg),dg.size(),
             100*iqr(dt),100*med(dt),dt.size());
      rows.push_back({P.sig, 100*iqr(du),100*med(du),(int)du.size(),
                             100*iqr(dg),100*med(dg),(int)dg.size(),
                             100*iqr(dt),100*med(dt),(int)dt.size()});
      f->Close();
   }
   printf(" -------------------------------------------------------------------------------------------\n");
   printf(" sig=IQR/1.349 resolution; med=bias (should be near 0 for an honest fit)\n\n");

   // plot sigma_p/p vs prfSigma for the ring points (skip base at sig=0 ring=0)
   std::vector<double> vx, vu, vg, vt;
   for (auto &r : rows) { vx.push_back(r.sig); vu.push_back(r.us); vg.push_back(r.gs); vt.push_back(r.ts); }
   if (vx.size() >= 2) {
      auto *c = new TCanvas("sw","PRF sweep",900,650);
      auto mk = [&](std::vector<double>&y,int col,int mst){ auto*g=new TGraph(vx.size(),vx.data(),y.data());
         g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mst); g->SetLineWidth(2); g->SetMarkerSize(1.4); return g; };
      auto *gu = mk(vu,kBlue+1,20), *gg = mk(vg,kRed+1,21), *gt = mk(vt,kGreen+2,22);
      auto *mg = new TMultiGraph(); mg->Add(gu); mg->Add(gg); mg->Add(gt);
      mg->SetTitle("PUMA resistive readout: #sigma_{p}/p vs Telegraph charge-dispersion #sigma;PRF #sigma [mm];#sigma_{p}/p [%]");
      mg->Draw("ALP"); mg->SetMinimum(0);
      auto *leg = new TLegend(0.55,0.70,0.88,0.88);
      leg->AddEntry(gu,"UKF (full pipeline)","lp");
      leg->AddEntry(gg,"genfit (full pipeline)","lp");
      leg->AddEntry(gt,"Taubin ring-centroid (offline ceiling)","lp");
      leg->Draw();
      c->SaveAs("./data/sweep_resolution.png");
      printf(" wrote ./data/sweep_resolution.png\n");
   }
}
