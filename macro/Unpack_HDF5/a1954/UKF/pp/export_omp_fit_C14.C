/// @file export_omp_fit_C14.C
/// @brief Dump everything the optical-model fit needs as plain text: the measured dsigma/dOmega,
///        the theta_cm migration matrix and the acceptance. Keeping the fit outside ROOT lets it
///        call FRESCO a few hundred times without paying ROOT startup each iteration.
void export_omp_fit_C14(TString corrCache = "plots/proton_kin_300gfx_nc_tc.root",
                        TString accDir = "/mnt/f/a1954_C14_acc_gf_nochi2/",
                        TString respFile = "", Double_t exWin = 0.6,
                        Double_t cmMin = 20.0, Double_t cmMax = 120.0, TString outDir = "plots/ompfit")
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (respFile.IsNull())
      respFile = here + "/../../../../Simulation/ATTPC/14C_pp/diagnostics/response_genfit.root";
   gSystem->mkdir(here + "/" + outDir, kTRUE);
   TFile *fc = TFile::Open(here + "/" + corrCache);
   TFile *fa = TFile::Open(accDir + "acceptance_merged_gs.root");
   TFile *frp = TFile::Open(respFile);
   if (!fc || fc->IsZombie() || !fa || fa->IsZombie() || !frp || frp->IsZombie()) { printf("missing input\n"); return; }
   TTree *t = (TTree *)fc->Get("pk");
   auto *acc = (TH1D *)fa->Get("hAcc_gs_sum");
   auto *R = (TH2D *)frp->Get("response");
   if (!t || !acc || !R) { printf("missing object\n"); return; }

   FILE *fd = fopen((here + "/" + outDir + "/data.txt").Data(), "w");
   fprintf(fd, "# theta_cm_center  dsdo(arb)  err   [selection |Ex|<%.2f, acceptance-corrected]\n", exWin);
   for (int b = 1; b <= acc->GetNbinsX(); ++b) {
      double lo = acc->GetBinLowEdge(b), wid = acc->GetBinWidth(b), ctr = acc->GetBinCenter(b);
      if (ctr < cmMin || ctr > cmMax) continue;
      auto *h = new TH1D(TString::Format("hq%d", b), "", 120, -3, 3);
      t->Draw(TString::Format("ex>>hq%d", b), TString::Format("thcm>=%g&&thcm<%g", lo, lo + wid), "goff");
      h->SetDirectory(nullptr);
      double y = h->Integral(h->FindBin(-exWin), h->FindBin(exWin));
      delete h;
      double A = acc->GetBinContent(b), s = std::sin(ctr * TMath::DegToRad());
      if (y <= 0 || A <= 0.05 || s <= 1e-3) continue;
      fprintf(fd, "%8.3f %14.6g %14.6g\n", ctr, y / A / s, std::sqrt(y) / A / s);
   }
   fclose(fd);

   FILE *fA = fopen((here + "/" + outDir + "/acceptance.txt").Data(), "w");
   for (int b = 1; b <= acc->GetNbinsX(); ++b)
      fprintf(fA, "%8.3f %10.5f\n", acc->GetBinCenter(b), acc->GetBinContent(b));
   fclose(fA);

   FILE *fR = fopen((here + "/" + outDir + "/response.txt").Data(), "w");
   fprintf(fR, "# nx ny  then rows: theta_true theta_reco P(reco|true)\n%d %d\n", R->GetNbinsX(), R->GetNbinsY());
   for (int x = 1; x <= R->GetNbinsX(); ++x)
      for (int y = 1; y <= R->GetNbinsY(); ++y)
         if (R->GetBinContent(x, y) > 0)
            fprintf(fR, "%8.3f %8.3f %12.6g\n", R->GetXaxis()->GetBinCenter(x), R->GetYaxis()->GetBinCenter(y),
                    R->GetBinContent(x, y));
   fclose(fR);
   printf("wrote %s/{data,acceptance,response}.txt\n", (here + "/" + outDir).Data());
}
