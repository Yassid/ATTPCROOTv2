/// @file refit_tracks_dp.C
/// @brief Refit a LIST of tracks under every cluster-ordering mode, so the browser viewer can
///        switch between them instantly. The browser cannot run genfit; it does not need to --
///        the fits are precomputed here and embedded in the page.
///
/// Ordering modes, matching fit_one_dp.C:
///   0  z-sort inside AtGenfitter  (the production; MUST reproduce the cache)
///   1  the clusteriser's own order (SetUseClusterOrder)
///   2  greedy nearest-neighbour arc walk computed here
///
/// The fitter configuration is copied from fit_one_dp.C, which is checked track-by-track against
/// dp_kin_dv1104.root. Note SetELossHybrid: omitting it silently makes this a different fitter
/// (KE 4.517 instead of 0.385 on run_0020/9789). Every setter here matters.
///
///   root -b -q 'refit_tracks_dp.C("list.txt","out.json")'
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

static std::vector<int> rfArcWalk(const std::vector<AtHitCluster> &hc)
{
   const int n = hc.size();
   std::vector<ROOT::Math::XYZPoint> p(n);
   for (int i = 0; i < n; ++i) p[i] = hc[i].GetPosition();
   double cx=0, cy=0, cz=0;
   for (auto &q : p) { cx+=q.X(); cy+=q.Y(); cz+=q.Z(); }
   cx/=n; cy/=n; cz/=n;
   auto d2=[&](int a,int b){ return (p[a].X()-p[b].X())*(p[a].X()-p[b].X())
        +(p[a].Y()-p[b].Y())*(p[a].Y()-p[b].Y())+(p[a].Z()-p[b].Z())*(p[a].Z()-p[b].Z()); };
   int far=0; double best=-1;
   for(int i=0;i<n;++i){ double d=(p[i].X()-cx)*(p[i].X()-cx)+(p[i].Y()-cy)*(p[i].Y()-cy)+(p[i].Z()-cz)*(p[i].Z()-cz);
      if(d>best){best=d;far=i;} }
   int start=far; best=-1;
   for(int i=0;i<n;++i) if(d2(far,i)>best){best=d2(far,i);start=i;}
   std::vector<char> used(n,0); std::vector<int> out; int cur=start;
   used[cur]=1; out.push_back(cur);
   for(int k=1;k<n;++k){ int nx=-1; double bd=1e30;
      for(int i=0;i<n;++i) if(!used[i]&&d2(cur,i)<bd){bd=d2(cur,i);nx=i;}
      if(nx<0) break; used[nx]=1; out.push_back(nx); cur=nx; }
   return out;
}

void refit_tracks_dp(TString listFile = "/mnt/f/a1975/caches/.explorer/refit_list.txt",
                     TString outJson  = "/mnt/f/a1975/caches/.explorer/refits_dp.json",
                     TString gfDir    = "/mnt/f/a1975/gf_dp_cateloss/",
                     TString geoName  = "ATTPC_D300torr_v2_geomanager.root",
                     TString eLossTable = "proton_D2_300torr.txt")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = gSystem->Getenv("VMCWORKDIR");
   { TFile *gf = TFile::Open(dir + "/geometry/" + geoName); gf->Get("FAIRGeom"); }
   TString elossPath = dir + "/macro/Unpack_HDF5/a1975/D2_UKF/" + eLossTable;

   auto fitter = std::make_unique<EventFit::AtGenfitter>(-2.85, 2212, 1.00782503207, 1,
                                                         elossPath.Data(), kFALSE, 2, 5);
   fitter->SetCatimaMaterial(kTRUE, kTRUE);
   fitter->SetCatimaELoss(kTRUE, kFALSE);
   fitter->SetELossHybrid(kTRUE, 6.61e-5);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(4.0);
   fitter->SetSeedFromSpyral(kFALSE);
   fitter->SetMatEffectsFallback(kFALSE);
   fitter->SetThetaWindow(10.0, 170.0);
   fitter->SetBackwardSeedFix(kTRUE);
   fitter->SetBackExtrapToAxis(kTRUE);
   fitter->Init();

   // group the requests by run so each 400 MB file is opened once
   std::map<TString, std::vector<std::pair<Long64_t,int>>> want;
   { std::ifstream in(listFile.Data()); TString r; Long64_t e; int tid;
     std::string rs;
     while (in >> rs >> e >> tid) want[TString(rs.c_str())].push_back({e,tid}); }
   long total=0; for (auto &kv : want) total += kv.second.size();
   printf("  %ld tracks over %zu runs\n", total, want.size());

   std::ofstream o(outJson.Data());
   o << "{";
   bool firstTrk = true; long done = 0;
   for (auto &kv : want) {
      TString fn = gfDir + "run_" + kv.first + "_multifit_genfitter_p.root";
      TFile *f = TFile::Open(fn);
      if (!f || f->IsZombie()) { printf("  [skip] %s\n", fn.Data()); continue; }
      auto *t = (TTree *)f->Get("cbmsim");
      TClonesArray *te = nullptr; t->SetBranchAddress("AtTrackingEvent", &te);
      for (auto &pr : kv.second) {
         t->GetEntry(pr.first);
         auto *ev = (AtTrackingEvent *)te->At(0);
         if (!ev) continue;
         auto trks = ev->GetTrackArray();
         AtTrack *src = nullptr;
         for (auto &tr : trks) if (tr.GetTrackID() == pr.second) src = &tr;
         if (!src) continue;
         std::vector<AtHitCluster> orig = *src->GetHitClusterArray();
         o << (firstTrk ? "" : ",") << "\"" << kv.first << ":" << pr.first << ":" << pr.second << "\":[";
         firstTrk = false;
         for (int mode = 0; mode < 3; ++mode) {
            AtTrack track = *src;
            auto *hc = track.GetHitClusterArray();
            if (mode == 2) { auto ord = rfArcWalk(orig);
               std::vector<AtHitCluster> re; re.reserve(ord.size());
               for (int i : ord) re.push_back(orig[i]);
               *hc = re; }
            else *hc = orig;
            fitter->SetUseClusterOrder(mode != 0);
            AtPatternEvent pe; pe.AddTrack(track);
            AtTrackingEvent tev;
            fitter->FitEvent(&tev, &pe, nullptr, nullptr, nullptr);
            o << (mode ? "," : "");
            if (tev.GetFittedTracks().empty()) { o << "null"; continue; }
            auto *ft = tev.GetFittedTracks().front().get();
            auto &k = ft->GetKinematicsXtr();
            double ndf = ft->GetTrackMetadata()->GetNdf(), c2 = ft->GetTrackMetadata()->GetChi2();
            auto &sp = ft->GetSmoothedPositions();
            o << "{\"ke\":" << Form("%.4f", k.kineticEnergy)
              << ",\"theta\":" << Form("%.2f", k.theta * TMath::RadToDeg())
              << ",\"chi2ndf\":" << Form("%.3f", ndf > 0 ? c2 / ndf : 1e9)
              << ",\"ndf\":" << Form("%.0f", ndf)
              << ",\"nfit\":" << sp.size() << ",\"fit\":[";
            for (size_t i = 0; i < sp.size(); ++i)     // z back to the CLUSTER frame
               o << (i ? "," : "") << "[" << Form("%.1f,%.1f,%.1f", sp[i].X(), sp[i].Y(), 1000.0 - sp[i].Z()) << "]";
            o << "]}";
         }
         o << "]";
         if (++done % 25 == 0) printf("  %ld / %ld\n", done, total);
      }
      f->Close();
   }
   o << "}";
   o.close();
   printf("  wrote %s  (%ld tracks x 3 modes)\n", outJson.Data(), done);
}
