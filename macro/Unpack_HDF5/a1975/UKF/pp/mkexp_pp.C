/// @file mkexp_pp.C
/// @brief Convert an a1975 (p,p') cache (pp/cache_pp_run.C layout) into the column names the
///        a2091 browser explorer expects, applying the IC beam gate on the way.
///
/// Our cache calls the vertex column `vz` and carries `ic`/`run`; make_explorer_html.C wants
/// `vertexz` and has no notion of an ion chamber, so the gate has to be baked in here.
///
///   root -b -q 'pp/mkexp_pp.C("in.root","out.root")'

void mkexp_pp(TString in, TString out, double chi2max = 1e9, double icMin = 950, double icMax = 1350)
{
   TFile *f = TFile::Open(in);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", in.Data()); return; }
   TTree *t = (TTree *)f->Get("pk");
   if (!t) { printf("no tree `pk` in %s\n", in.Data()); return; }
   // vertex and ion-chamber columns are optional: the D2 caches name the vertex `vertexz`, and the
   // D2 reco has no _FRIB.root at all so there is no `ic`. A flat vertexz makes the page hide its
   // z panel, which is the honest outcome when the column was never filled.
   float ke, theta, vz = 0, chi2ndf, ic = -1;
   t->SetBranchAddress("ke",&ke); t->SetBranchAddress("theta",&theta);
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

   Long64_t n = 0;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ke <= 0 || chi2ndf > chi2max) continue;
      if (icMin > 0 && (ic < icMin || ic > icMax)) continue;
      oke = ke; oth = theta; ovz = vz; oc2 = chi2ndf;
      pk->Fill(); ++n;
   }
   o.cd(); pk->Write(); o.Close(); f->Close();
   printf("mkexp_pp: %lld -> %lld tracks (IC[%.0f,%.0f]) -> %s\n",
          t->GetEntries(), n, icMin, icMax, out.Data());
}
