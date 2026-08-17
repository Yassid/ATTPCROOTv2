/// @file mkexp_pp.C
/// @brief Convert an a1975 (p,p') cache (pp/cache_pp_run.C layout) into the column names the
///        a2091 browser explorer expects, applying the IC beam gate on the way.
///
/// Our cache calls the vertex column `vz` and carries `ic`/`run`; make_explorer_html.C wants
/// `vertexz` and has no notion of an ion chamber, so the gate has to be baked in here.
///
///   root -b -q 'pp/mkexp_pp.C("in.root","out.root")'

/// keCol/thCol pick WHICH energy column to build the page from. The (d,t) cache carries two per
/// track since AtGenfitter started writing them separately: "ke"/"theta" are back-extrapolated to
/// the beam axis (plus CATIMA over the vertex gap), "kefit"/"thetafit" are what the fit returned
/// at the first measurement point. Pointing two slots at the same cache with different columns is
/// how the two are compared without deriving a second file.
void mkexp_pp(TString in, TString out, double chi2max = 1e9, double icMin = 950, double icMax = 1350,
              TString keCol = "ke", TString thCol = "theta")
{
   TFile *f = TFile::Open(in);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", in.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree `pk` in %s\n", in.Data()); return; }
   // vertex and ion-chamber columns are optional: the D2 caches name the vertex `vertexz`, and the
   // D2 reco has no _FRIB.root at all so there is no `ic`. A flat vertexz makes the page hide its
   // z panel, which is the honest outcome when the column was never filled.
   float ke, theta, vz = 0, chi2ndf, ic = -1;
   if (!t->GetBranch(keCol) || !t->GetBranch(thCol)) {
      printf("mkexp_pp: no `%s`/`%s` in %s\n", keCol.Data(), thCol.Data(), in.Data()); return;
   }
   t->SetBranchAddress(keCol,&ke); t->SetBranchAddress(thCol,&theta);
   t->SetBranchAddress("chi2ndf",&chi2ndf);
   if (t->GetBranch("vz")) t->SetBranchAddress("vz",&vz);
   else if (t->GetBranch("vertexz")) t->SetBranchAddress("vertexz",&vz);
   const bool hasIC = t->GetBranch("ic") != nullptr;
   if (hasIC) t->SetBranchAddress("ic",&ic);
   if (!hasIC && icMin > 0) { printf("mkexp_pp: no `ic` branch -- IC gate disabled\n"); icMin = -1; }

   TFile o(out, "RECREATE");
   float oke, oth, ovz, oc2;
   TTree *pk = new TTree("pk", "16C(p,p') kinematics for the browser explorer");
   pk->Branch("ke",&oke); pk->Branch("theta",&oth);
   pk->Branch("vertexz",&ovz); pk->Branch("chi2ndf",&oc2);

   // COLLAPSED FITS ARE DROPPED UNCONDITIONALLY, and this cannot be left to chi2max.
   // ex_dt_a1975 writes chi2ndf = 1e9 EXACTLY when ndf <= 0, i.e. when the Kalman filter kept
   // essentially no measurement. Such a track still carries kinematics and reaches the ntuple
   // like any other. Every caller passes chi2max = 1e9 to mean "no chi2 cut", and
   // `chi2ndf > chi2max` is then 1e9 > 1e9, which is FALSE -- so the sentinel sailed through
   // and the explorer pages were built on it. Measured on dt_kin_maton.root (material effects
   // on, genfit's own Highland): 22433 of 37127 tracks, 60.4%, were collapsed fits being shown
   // as data. dt_kin_dv1104.root (material effects off) is only 0.3%, which is why this went
   // unnoticed for so long -- it is invisible until material effects are turned on.
   const double kCollapsed = 1e9;
   Long64_t nCollapsed = 0;

   Long64_t n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (chi2ndf >= kCollapsed) { ++nCollapsed; continue; }
      if (ke <= 0 || chi2ndf > chi2max) continue;
      if (icMin > 0 && (ic < icMin || ic > icMax)) continue;
      oke = ke; oth = theta; ovz = vz; oc2 = chi2ndf;
      pk->Fill(); ++n;
   }
   o.cd(); pk->Write(); o.Close(); f->Close();
   // report the collapsed count rather than swallowing it: a cache that is mostly collapsed
   // fits is a finding about the production, not a detail of this converter
   printf("mkexp_pp: %lld -> %lld tracks (IC[%.0f,%.0f], %lld collapsed ndf<=0 dropped = %.1f%%) -> %s\n",
          t->GetEntries(), n, icMin, icMax, nCollapsed,
          t->GetEntries() ? 100.0 * nCollapsed / t->GetEntries() : 0.0, out.Data());
}
