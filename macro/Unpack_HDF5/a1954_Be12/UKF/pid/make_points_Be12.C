/// @file make_points_Be12.C
/// @brief Turn the PID csv (pid/dump_pid_Be12.C output) into the "pts" tree the buttoned gate
///        drawers expect, so gate_draw_pt_Be12.C can load it in a second instead of re-reading
///        40 GB of reco. Same branch names/types as a1975 pid_plane_dt.C.
///
///   root -b -q 'make_points_Be12.C("pid_Be12_data_mp30.csv","plots/points_Be12_mp30.root")'
void make_points_Be12(TString csvName = "pid_Be12_data_mp30.csv", TString outName = "plots/points_Be12_mp30.root")
{
   TString dir = gSystem->DirName(__FILE__);
   TString csv = dir + "/" + csvName;
   if (gSystem->AccessPathName(csv)) { printf("\033[1;31mno %s\033[0m\n", csv.Data()); return; }
   TTree in("in", "in");
   in.ReadFile(csv, "", ',');   // sqrtdedx,brho,dedx,polar,arclen,npts,ic
   Float_t sd, br, de, po, al, np_, ic;
   in.SetBranchAddress("sqrtdedx", &sd); in.SetBranchAddress("brho", &br);
   in.SetBranchAddress("dedx", &de);     in.SetBranchAddress("polar", &po);
   in.SetBranchAddress("arclen", &al);   in.SetBranchAddress("npts", &np_);
   in.SetBranchAddress("ic", &ic);

   gSystem->mkdir((dir + "/plots").Data(), kTRUE);
   TFile f(dir + "/" + outName, "RECREATE");
   TTree t("pts", "PID points");
   Float_t o_sd, o_br, o_pol, o_arc, o_ic;
   t.Branch("sqrtdedx", &o_sd); t.Branch("brho", &o_br); t.Branch("polar", &o_pol);
   t.Branch("arclen", &o_arc);  t.Branch("ic", &o_ic);
   Long64_t n = in.GetEntries(), kept = 0;
   for (Long64_t i = 0; i < n; ++i) {
      in.GetEntry(i);
      // polar is ALREADY IN DEGREES in the csv -- dump_pid_Be12.C writes r.polar*RadToDeg().
      // Converting again put the maximum at 10311 "degrees"; the drawers want degrees, as here.
      o_sd = sd; o_br = br; o_pol = po; o_arc = al; o_ic = ic;
      if (!std::isfinite(o_br) || o_br <= 0 || o_br > 50) continue;  // straight tracks -> huge brho
      t.Fill(); ++kept;
   }
   t.Write(); f.Close();
   printf("\033[1;32mpts tree: %lld / %lld tracks -> %s/%s\033[0m\n", kept, n, dir.Data(), outName.Data());
}
