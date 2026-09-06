/// @file dump_tracks_dp.C
/// @brief Dump the HIT CLOUD of selected 16C(d,p)17C tracks to JSON, for the browser track viewer.
///
/// The aggregate explorer answers "what does this population look like"; this answers "is THIS
/// track real". It exists because the backward sub-MeV selection (theta_lab > 90, KE < 1 MeV)
/// is short-track territory -- a 1 MeV proton ranges 264 mm and a 0.5 MeV one only 86 mm -- and
/// no histogram can tell a genuine short track from a fragment of a curler or a noise cluster.
///
/// SELECTION is applied on the SAME quantities the cache uses (GetKinematicsXtr, chi2/ndf from
/// the metadata, the IC file entry-aligned with the fit tree), so a track shown here is a track
/// in dp_kin_dv1104.root, not an approximation of one.
///
/// DUPLICATE ENTRIES: the genfitter trees repeat ~22% of events, bit-identical, in the NEXT
/// entry. The copy carries an EMPTY AtPIDEvent, so requiring a valid Spyral record -- which the
/// selection needs anyway for brho/dEdx -- skips the copies for free. Without that the viewer
/// would show the same track twice and it would look like two independent confirmations.
///
///   root -b -q 'dump_tracks_dp.C("out.json", 240, 90, 180, 0, 1.0)'
void dump_tracks_dp(TString outJson = "tracks_dp_backward.json", int nWant = 240,
                    double thMin = 90.0, double thMax = 180.0,
                    double keMin = 0.0, double keMax = 1.0,
                    double icLo = 900, double icHi = 1300, double chi2Max = 1e8,
                    /// Qualifying-track count, if already known (39242-track cache says 3150 for
                    /// the backward sub-MeV window). Passing it SKIPS the counting pass, which is
                    /// half the I/O -- and the I/O here is 19 GB over drvfs, so that is minutes.
                    long knownTotal = 0,
                    /// Take every runStride-th run. The sample still spans 0016-0103 at a fraction
                    /// of the reading. 1 = every run.
                    int runStride = 1,
                    TString gfDir = "/mnt/f/a1975/gf_dp_cateloss/",
                    TString icDir = "/mnt/f/a1975/ic_d2/", double Ebeam = 184.25)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   const double u = 931.49401;
   const double m1 = 16.0147013*u, m2 = 2.0135532*u, m3 = 1.00727646688*u, m4 = 17.0225864*u;
   auto om=[](double x,double y,double z){return std::sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);};
   auto kine=[&](double th,double ke)->std::pair<double,double>{
      double Et1=Ebeam+m1, Et3=ke+m3;
      double s=m1*m1+m2*m2+2*m2*Et1, uu=m2*m2+m3*m3-2*m2*Et3;
      double arg=(std::cos(th)*om(s,m1*m1,m2*m2)*om(uu,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-uu))/(2*m2*m2)+s+uu-m2*m2;
      if(arg<0) return {std::nan(""),std::nan("")};
      double m4x=std::sqrt(arg);
      double t=m2*m2+m4x*m4x-2*m2*(Et1+m2-Et3);
      double tcm=TMath::Pi()-std::acos((s*s+s*(2*t-m1*m1-m2*m2-m3*m3-m4x*m4x)+(m1*m1-m2*m2)*(m3*m3-m4x*m4x))/
                                       (om(s,m1*m1,m2*m2)*om(s,m3*m3,m4x*m4x)));
      return {m4x-m4, tcm*TMath::RadToDeg()};
   };

   const char *runs[] = {"0016","0017","0018","0019","0020","0021","0022","0023","0026","0027","0031","0032",
      "0034","0036","0037","0038","0039","0040","0041","0042","0043","0044","0046","0048","0057","0058","0076",
      "0077","0078","0079","0080","0082","0083","0084","0085","0086","0087","0088","0089","0091","0092","0095",
      "0096","0097","0098","0102","0103"};

   // FIRST PASS counts the qualifying tracks so the sample can be spread evenly over all 47 runs.
   // Taking the first nWant would sample only the earliest runs, and run-to-run differences are
   // exactly the thing a by-eye check should be exposed to.
   long total = knownTotal;
   std::vector<std::string> out;
   for (int pass = (knownTotal > 0 ? 1 : 0); pass < 2; ++pass) {
      long seen = 0; int stride = 1;
      if (pass == 1) { stride = (total/runStride > nWant) ? (int)(total/runStride/nWant) : 1;
                       printf("  qualifying tracks: %ld -> taking every %d\n", total, stride); }
      int ri = -1;
      for (const char *rn : runs) {
         if (++ri % runStride) continue;
         TString gf = gfDir + TString("run_") + rn + "_multifit_genfitter_p.root";
         if (gSystem->AccessPathName(gf)) continue;
         TFile *fu = TFile::Open(gf); if(!fu||fu->IsZombie()) continue;
         TTree *tu = (TTree*)fu->Get("cbmsim");
         TClonesArray *te=nullptr,*pe=nullptr;
         tu->SetBranchAddress("AtTrackingEvent",&te); tu->SetBranchAddress("AtPIDEvent",&pe);
         std::vector<float> icv;
         TString icf = icDir + TString("run_") + rn + "_IC.root";
         if(!gSystem->AccessPathName(icf)){ TFile*fi=TFile::Open(icf);
            if(fi&&!fi->IsZombie()){ TTree*ti=(TTree*)fi->Get("ic"); float v=0; ti->SetBranchAddress("icmax",&v);
               for(Long64_t j=0;j<ti->GetEntries();++j){ti->GetEntry(j); icv.push_back(v);} fi->Close(); } }
         for (Long64_t i=0;i<tu->GetEntries();++i){
            tu->GetEntry(i);
            if(!pe||pe->GetEntries()==0||!te||te->GetEntries()==0) continue;   // skips the duplicates
            auto*pidev=(AtPIDEvent*)pe->At(0); auto*ev=(AtTrackingEvent*)te->At(0);
            if(!pidev||!ev) continue;
            float ic = (i<(Long64_t)icv.size())?icv[i]:-1.f;
            if(icLo>0 && (ic<icLo||ic>icHi)) continue;
            auto trkArr = ev->GetTrackArray();
            std::map<int,AtFittedTrack*> fmap;
            for(auto&ft:ev->GetFittedTracks()) if(ft) fmap[ft->GetTrackID()]=ft.get();
            for(auto&sr:pidev->GetSpyral()){
               if(!sr.valid) continue;
               auto it=fmap.find(sr.trackID); if(it==fmap.end()) continue;
               auto*ft=it->second; auto&k=ft->GetKinematicsXtr();
               double ndf=ft->GetTrackMetadata()->GetNdf(), c2=ft->GetTrackMetadata()->GetChi2();
               double c2n = ndf>0 ? c2/ndf : 1e9;
               double ke=k.kineticEnergy, thDeg=k.theta*TMath::RadToDeg();
               if(ke<=keMin||ke>keMax||c2n>chi2Max) continue;
               if(thDeg<thMin||thDeg>thMax) continue;
               ++seen;
               if(pass==0){ ++total; continue; }
               if(seen % stride) continue;
               if((int)out.size()>=nWant) break;
               // the hit cloud of the PATTERN track this fit came from
               AtTrack *pt=nullptr;
               for(auto&tr:trkArr) if(tr.GetTrackID()==sr.trackID){ pt=&tr; break; }
               if(!pt) continue;
               auto*hc=pt->GetHitClusterArray();
               if(!hc||hc->size()<3) continue;
               auto v=ft->GetVertex();
               auto [ex,tcm]=kine(k.theta,ke);
               // THE FIT ITSELF. AtGenfitter stores the fitted state at every measurement point
               // (AtGenfitter.cxx:670, "fitted trajectory for display/QA"), so this is the actual
               // Kalman trajectory, not a helix re-derived from the summary quantities -- which
               // would hide exactly the failure the viewer exists to catch, a fit that does not
               // follow its own hits. Its LENGTH is also diagnostic: fewer fit points than
               // clusters means the fit dropped measurements.
               // *** THE Z FRAMES ARE MIRRORED AND THIS IS NOT OPTIONAL ***
               // GetSmoothedPositions comes back in the FITTER's frame, where the pad plane is at
               // z = 1000 (AtGenfitter's hardcoded SetZPadPlane(1000.0)), while AtHitCluster z runs
               // the other way. MEASURED on 2641 fit points matched to clusters by (x,y) ALONE, so
               // z could not bias the match:  z_fit + z_cluster = 1000.00.
               // Plotting them on one axis without this flip put the trajectory 300 mm off its own
               // hits: median hit-fit rms was 306 mm for fits with chi2/ndf < 2, and 5.1 mm after.
               // x and y need no transformation -- they agreed to 0.1 mm all along, which is what
               // makes the mistake look like a physics effect instead of a frame error.
               const double kZPadFit = 1000.0;
               auto &sp = ft->GetSmoothedPositions();
               std::string fitpts="[";
               for(size_t c=0;c<sp.size();++c)
                  fitpts += TString::Format("%s[%.1f,%.1f,%.1f]", c?",":"",
                                            sp[c].X(), sp[c].Y(), kZPadFit - sp[c].Z()).Data();
               fitpts += "]";
               std::string pts="[";
               for(size_t c=0;c<hc->size();++c){ auto p=(*hc)[c].GetPosition();
                  pts += TString::Format("%s[%.1f,%.1f,%.1f,%.0f]", c?",":"", p.X(),p.Y(),p.Z(),(*hc)[c].GetCharge()).Data(); }
               pts += "]";
               out.push_back(TString::Format(
                  "{\"run\":\"%s\",\"entry\":%lld,\"tid\":%d,\"ke\":%.3f,\"theta\":%.2f,\"phi\":%.1f,"
                  "\"ex\":%.3f,\"thcm\":%.1f,\"brho\":%.4f,\"dedx\":%.1f,\"chi2ndf\":%.3f,\"ncl\":%d,"
                  "\"ic\":%.0f,\"vx\":%.1f,\"vy\":%.1f,\"vz\":%.1f,\"nfit\":%d,\"fit\":%s,\"hits\":%s}",
                  rn, i, sr.trackID, ke, thDeg, k.phi*TMath::RadToDeg(),
                  std::isnan(ex)?-999:ex, std::isnan(tcm)?-999:tcm, sr.brho, sr.dEdx, c2n,
                  (int)hc->size(), ic, v.X(), v.Y(), v.Z(), (int)sp.size(), fitpts.c_str(),
                  pts.c_str()).Data());
            }
            if(pass==1 && (int)out.size()>=nWant) break;
         }
         fu->Close();
         if(pass==1 && (int)out.size()>=nWant) break;
      }
   }
   std::ofstream f(outJson.Data());
   f << "{\"sel\":\"theta_lab " << thMin << "-" << thMax << " deg, KE " << keMin << "-" << keMax
     << " MeV, IC [" << icLo << "," << icHi << "], chi2/ndf < " << chi2Max << "\","
     << "\"ebeam\":" << Ebeam << ",\"ntotal\":" << total << ",\"tracks\":[";
   for(size_t i=0;i<out.size();++i) f << (i?",":"") << out[i];
   f << "]}";
   f.close();
   printf("  wrote %s : %zu tracks of %ld qualifying\n", outJson.Data(), out.size(), total);
}
