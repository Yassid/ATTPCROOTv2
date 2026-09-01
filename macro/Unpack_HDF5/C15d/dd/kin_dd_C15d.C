/// @file kin_dd_C15d.C
/// @brief (d,d') kinematics: select the gated deuterons from the fitted sample and map KE vs theta.
///
///   root -b -q 'dd/kin_dd_C15d.C()'
///
/// Chains every <run>_kin_d.root, keeps the tracks in the deuteron selection (joined on
/// run/event/trackID), and writes the KE-vs-theta_lab map plus a compact `dd` ntuple.
///
/// ★ WHAT THE MAP IS FOR. 15C(d,d)15C elastic is a two-body reaction with the beam energy as its
/// only free parameter, so the recoil-deuteron ridge KE(theta) MEASURES the beam energy:
///
///     T_d = 2 m_d p1^2 cos^2(theta) / [ (E1 + m_d)^2 - p1^2 cos^2(theta) ]
///
/// in inverse kinematics with a heavy beam on a light target. That is the same construction the
/// a2091 15C(p,p') analysis used on its proton ridge, and it needs no assumed beam energy to start
/// from -- which is why it is worth doing even though the number is expected from elsewhere: an
/// independent measurement is a check on the value, not a duplicate of it.
///
/// ★ theta IS NOT FOLDED. Elastic recoils in inverse kinematics are forward, but the sample is not
/// pre-restricted to the forward hemisphere: a backward population here would mean the ejectile
/// assignment or the seeding is wrong, and folding would hide it.

void kin_dd_C15d(TString fitDir = "/home/yassid/C15d_fit/", TString selFile = "pid/sel_deuteron_C15d.root",
                 TString outDir = "dd/plots/", Double_t chi2Cut = 5.0, Double_t keMax = 60.0,
                 /// The IC value is carried THROUGH to the ntuple, not cut on here, so the beam
                 /// window can be moved in the viewer instead of being frozen at cache-build time.
                 /// Every earlier IC decision was baked into a selection file, which meant changing
                 /// the window cost a full rebuild and there was no way to see what a window was
                 /// actually selecting before committing to it.
                 TString pointsFile = "pid/points_C15d.root")
{
   gSystem->mkdir(outDir, kTRUE);

   // ---- selection ---------------------------------------------------------------------------
   std::set<Long64_t> keep;
   auto key = [](int r, int e, int t) { return ((Long64_t)r << 44) | ((Long64_t)e << 12) | (Long64_t)(t & 0xFFF); };
   if (selFile.Length()) {
      if (gSystem->AccessPathName(selFile)) {
         std::cout << "\033[1;31mERROR: " << selFile << " not found (run apply_gate_C15d.C).\033[0m\n";
         return;
      }
      TFile *fs = TFile::Open(selFile);
      TTree *ts = (TTree *)fs->Get("sel");
      Int_t r, e, t;
      ts->SetBranchAddress("run", &r);
      ts->SetBranchAddress("event", &e);
      ts->SetBranchAddress("trackID", &t);
      for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
         ts->GetEntry(i);
         keep.insert(key(r, e, t));
      }
      fs->Close();
      std::cout << "  selection : " << keep.size() << " gated deuteron tracks\n";
   }

   // ---- IC join ------------------------------------------------------------------------------
   // Same (run, event, trackID) key as the selection above. npulse rides along because the IC
   // window was chosen on the single-pulse spectrum, so a viewer that lets the window move has to
   // be able to reproduce the single-pulse condition too.
   std::map<Long64_t, std::pair<float, int>> icMap;
   if (pointsFile.Length() && !gSystem->AccessPathName(pointsFile)) {
      TFile *fp = TFile::Open(pointsFile);
      TTree *tp = fp ? (TTree *)fp->Get("pts") : nullptr;
      if (tp) {
         Int_t pr, pe, pt, pn;
         Float_t pic;
         tp->SetBranchAddress("run", &pr);
         tp->SetBranchAddress("event", &pe);
         tp->SetBranchAddress("trackID", &pt);
         tp->SetBranchAddress("npulse", &pn);
         tp->SetBranchAddress("ic", &pic);
         for (Long64_t i = 0; i < tp->GetEntries(); ++i) {
            tp->GetEntry(i);
            icMap[key(pr, pe, pt)] = {pic, pn};
         }
      }
      fp->Close();
      std::cout << "  IC join   : " << icMap.size() << " tracks carry an IC value\n";
   } else {
      std::cout << "\033[1;33m  IC join   : " << pointsFile
                << " not found -- ic will be -1 and the viewer's IC control stays disabled\033[0m\n";
   }

   TChain ch("kin");
   const int nf = ch.Add(fitDir + "*_kin_d.root");
   Int_t run, event, track, ndf, fwd;
   Double_t ke, th, keX, thX, vz, c2;
   ch.SetBranchAddress("run", &run);
   ch.SetBranchAddress("event", &event);
   ch.SetBranchAddress("trackID", &track);
   ch.SetBranchAddress("ke", &ke);
   ch.SetBranchAddress("theta", &th);
   ch.SetBranchAddress("keXtr", &keX);
   ch.SetBranchAddress("thetaXtr", &thX);
   ch.SetBranchAddress("vz", &vz);
   ch.SetBranchAddress("chi2ndf", &c2);
   ch.SetBranchAddress("ndf", &ndf);
   ch.SetBranchAddress("dirFwd", &fwd);

   auto *h = new TH2D("hdd", "15C(d,d') gated deuterons;#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 240, 0, keMax);
   auto *hAll = new TH2D("hddAll", "all fitted deuterons;#theta_{lab} [deg];KE [MeV]", 180, 0, 180, 240, 0, keMax);

   TFile fo(outDir + "dd_kin_C15d.root", "RECREATE");
   TTree dd("dd", "gated (d,d') kinematics");
   dd.Branch("run", &run, "run/I");
   dd.Branch("event", &event, "event/I");
   dd.Branch("trackID", &track, "trackID/I");
   dd.Branch("ke", &keX, "ke/D");
   dd.Branch("theta", &thX, "theta/D");
   dd.Branch("vz", &vz, "vz/D");
   dd.Branch("chi2ndf", &c2, "chi2ndf/D");
   Float_t icv = -1;
   Int_t npul = 0;
   dd.Branch("ic", &icv, "ic/F");
   dd.Branch("npulse", &npul, "npulse/I");

   Long64_t nAll = 0, nSel = 0, nBack = 0;
   for (Long64_t i = 0; i < ch.GetEntries(); ++i) {
      ch.GetEntry(i);
      if (!(keX > 0) || keX > keMax * 10 || c2 > chi2Cut)
         continue;
      ++nAll;
      hAll->Fill(thX, keX);
      if (!keep.empty() && keep.find(key(run, event, track)) == keep.end())
         continue;
      ++nSel;
      if (thX > 90)
         ++nBack;
      h->Fill(thX, keX);
      {
         auto it = icMap.find(key(run, event, track));
         icv = (it == icMap.end()) ? -1.f : it->second.first;
         npul = (it == icMap.end()) ? 0 : it->second.second;
      }
      dd.Fill();
   }
   fo.cd();
   dd.Write();
   h->Write();
   hAll->Write();

   auto *c = new TCanvas("cdd", "dd", 1100, 850);
   c->SetLogz();
   c->SetRightMargin(0.13);
   h->Draw("colz");
   c->SaveAs(outDir + "dd_kin_C15d.png");
   fo.Close();

   std::cout << "\033[1;33m=== (d,d') kinematics ===\033[0m\n"
             << "  kin files  : " << nf << "\n"
             << "  fitted     : " << nAll << "  (keXtr>0, chi2/ndf<" << chi2Cut << ")\n"
             << "  in gate    : " << nSel << "\n"
             << "  backward   : " << nBack << " = " << (nSel ? 100.0 * nBack / nSel : 0.)
             << "%  (elastic recoils in inverse kinematics should be FORWARD -- a large backward "
                "fraction means the ejectile or the seeding is wrong)\n"
             << "  \033[1;32mwrote\033[0m " << outDir << "dd_kin_C15d.{root,png}\n";
}
