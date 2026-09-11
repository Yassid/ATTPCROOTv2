// CSV -> TTree for the Dec2014 alpha kinematics drawer.
// Tracks are stored as NEAR/FAR endpoints in (x, y, tb) so the drawer can rebuild the
// 3D geometry for ANY drift velocity:  z = vD * (tb - t0) * 160 ns * 10  [mm].
void kin_points(const char *in = "kin_points.csv", const char *out = "kin_points.root")
{
   std::ifstream f(in);
   if (!f.is_open()) { printf("ERROR cannot open %s\n", in); return; }
   std::string line; std::getline(f, line);
   Int_t run; Float_t b[6], a1[6], a2[6];
   TFile fo(out, "RECREATE");
   TTree *t = new TTree("kin", "Dec2014 alpha track endpoints");
   t->Branch("run", &run, "run/I");
   t->Branch("b", b, "b[6]/F"); t->Branch("a1", a1, "a1[6]/F"); t->Branch("a2", a2, "a2[6]/F");
   Long64_t n = 0;
   while (std::getline(f, line)) {
      std::stringstream ss(line); std::string tok; std::vector<double> v;
      while (std::getline(ss, tok, ',')) v.push_back(atof(tok.c_str()));
      if (v.size() < 19) continue;
      run = (Int_t)v[0];
      for (int i = 0; i < 6; ++i) { b[i] = v[1+i]; a1[i] = v[7+i]; a2[i] = v[13+i]; }
      t->Fill(); ++n;
   }
   t->Write(); fo.Close();
   printf("wrote %s : %lld events\n", out, n);
}
