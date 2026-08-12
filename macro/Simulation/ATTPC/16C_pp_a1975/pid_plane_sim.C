/// @file pid_plane_sim.C
/// @brief The simulated PID plane. Nothing else -- no gate, no data, no truth matching.
///
///   root -b -q 'pid_plane_sim.C("/mnt/f/a1975_C16_pp_pid","s2001,s2002,s2003,s2004,s2005,s2006")'

void pid_plane_sim(TString dir = "/mnt/f/a1975_C16_pp_pid", TString tags = "s2001,s2002,s2003,s2004,s2005,s2006",
                   TString png = "plots/pid_plane_sim.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);

   auto *h = new TH2D("hSim", "simulated PID, 16C(p,p');#sqrt{dE/dx} [arb];B#rho [T#upointm]", 300, 0, 40, 250, 0, 1.2);
   long n = 0;

   TObjArray *ta = tags.Tokenize(",");
   for (int it = 0; it < ta->GetEntries(); ++it) {
      TString tg = ((TObjString *)ta->At(it))->GetString().Strip(TString::kBoth);
      TString fn = dir + "/" + tg + "_pid.root";
      if (gSystem->AccessPathName(fn)) {
         printf("  skip %s (missing)\n", tg.Data());
         continue;
      }
      TFile *f = TFile::Open(fn);
      TTree *t = f ? (TTree *)f->Get("cbmsim") : nullptr;
      if (!t) {
         if (f) f->Close();
         continue;
      }
      TClonesArray *pe = nullptr;
      t->SetBranchAddress("AtPIDEvent", &pe);
      long nt = 0;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (!pe || pe->GetEntriesFast() == 0)
            continue;
         auto *ev = (AtPIDEvent *)pe->At(0);
         if (!ev)
            continue;
         for (auto &sp : ev->GetSpyral()) {
            h->Fill(sp.sqrtdEdx, sp.brho);
            ++nt;
         }
      }
      printf("  %-8s %7ld entries\n", tg.Data(), nt);
      n += nt;
      f->Close();
   }
   printf("\n  %ld PID entries total\n", n);

   TCanvas *c = new TCanvas("cS", "simulated PID", 900, 700);
   gPad->SetLogz();
   gPad->SetRightMargin(0.13);
   h->Draw("colz");
   gSystem->mkdir(gSystem->DirName(png), kTRUE);
   c->SaveAs(png);
   printf("  wrote %s\n", png.Data());
}
