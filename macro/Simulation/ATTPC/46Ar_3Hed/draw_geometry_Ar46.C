/// @file draw_geometry_Ar46.C
/// @brief To-scale side view of the 46Ar setup: the AT-TPC drift volume and the forward telescope.
///
///   root -b -q 'draw_geometry_Ar46.C()'
///
/// EVERY BOX DRAWN HERE IS READ OUT OF THE GEOMETRY FILES, not typed from the design. The point
/// of the picture is to check that the telescope is where it is supposed to be relative to the
/// gas, so a drawing made from the same numbers the design was written with would be worthless --
/// it would agree with itself no matter what the simulation actually contains.
///
/// The two geometries are separate FairModules and are merged only inside FairRunSim, so they are
/// opened separately here and overlaid in one z-r frame, which is exactly how Geant4 sees them.
void draw_geometry_Ar46(TString outPng = "plots/geometry_Ar46.png")
{
   gStyle->SetOptStat(0);
   TString dir = gSystem->Getenv("VMCWORKDIR");
   if (dir.IsNull()) dir = "/home/yassid/fair_install/ATTPCROOTv2-OpenKF";

   struct Box { TString name; double z0, z1, r0, r1; bool slab = false;
                double xc = 0, yc = 0, dx = 0, dy = 0; };
   std::vector<Box> boxes;

   // Walk one geometry file's top volume and record every placed node as (z range, r range).
   auto scan = [&](TString file, TString tag) {
      TFile *f = TFile::Open(dir + "/geometry/" + file);
      if (!f || f->IsZombie()) { printf("  cannot open %s\n", file.Data()); return; }
      auto *gm = (TGeoManager *)f->Get(f->GetListOfKeys()->First()->GetName());
      if (!gm) { printf("  no TGeoManager in %s\n", file.Data()); return; }
      TGeoVolume *top = gm->GetTopVolume();
      // two levels: the version assembly sits under TOP
      std::function<void(TGeoNode *, double)> walk = [&](TGeoNode *nd, double zoff) {
         auto *v = nd->GetVolume();
         const double *tr = nd->GetMatrix()->GetTranslation();
         double z = zoff + tr[2];
         auto *bb = (TGeoBBox *)v->GetShape();
         bool assembly = v->IsAssembly();
         if (!assembly && bb) {
            double dz = bb->GetDZ(), dx = bb->GetDX(), dy = bb->GetDY();
            // THE TRANSVERSE OFFSET MATTERS. An array element is placed off-axis, and ignoring
            // tr[0] drew all 36 CsI elements stacked on the beam axis -- a picture that hid
            // exactly the thing the array is for (covering the silicon face). x is used as the
            // vertical here since this is an x-z section.
            double xc = tr[0], yc = tr[1];
            double rmin = 0, rmax = std::max(dx, dy);
            bool tube = TString(v->GetShape()->ClassName()) == "TGeoTube";
            if (tube) {
               auto *tb = (TGeoTube *)v->GetShape();
               rmin = tb->GetRmin(); rmax = tb->GetRmax();
            }
            // a tube is symmetric about the axis; a box is a slab centred on xc
            if (tube)
               boxes.push_back({TString(v->GetName()), z - dz, z + dz, rmin, rmax});
            else
               boxes.push_back({TString(v->GetName()), z - dz, z + dz, xc - dx, xc + dx, true,
                                xc, yc, dx, dy});
         }
         for (int i = 0; i < v->GetNdaughters(); ++i) walk(v->GetNode(i), z);
      };
      for (int i = 0; i < top->GetNdaughters(); ++i) walk(top->GetNode(i), 0.0);
      printf("  %-42s -> %zu volumes\n", (tag + " (" + file + ")").Data(), boxes.size());
      f->Close();
   };

   printf("\nreading the geometry files:\n");
   size_t n0 = 0;
   scan("ATTPC_He3CO2_300torr_geomanager.root", "AT-TPC");
   n0 = boxes.size();
   scan("Ar46_telescope_v1.0_geomanager.root", "telescope");

   if (boxes.empty()) { printf("  nothing to draw\n"); return; }

   TCanvas *c = new TCanvas("cg", "Ar46 geometry", 1500, 620);
   c->Divide(1, 2, 0.001, 0.001);

   auto drawFrame = [&](double zlo, double zhi, double rhi, const char *title) {
      TH1F *fr = gPad->DrawFrame(zlo, -rhi, zhi, rhi);
      fr->SetTitle(Form("%s;z along the beam [cm];x or r [cm]", title));
      for (size_t i = 0; i < boxes.size(); ++i) {
         const auto &b = boxes[i];
         const bool tel = (i >= n0);
         if (b.z1 < zlo || b.z0 > zhi) continue;
         int col = !tel ? kGray + 1
                        : b.name.Contains("_dE") ? kGreen + 2
                        : b.name.Contains("CsI") ? kOrange - 3 : kAzure + 2;
         // A SOLID volume (rmin = 0) spans -rmax..+rmax as ONE band. Only a hollow one (a vessel
         // shell, rmin > 0) is two bands. Drawing the solid case as a single positive band was
         // wrong and made the DSSDs look like half-detectors.
         std::vector<std::pair<double, double>> bands;
         if (b.slab)
            bands.push_back({b.r0, b.r1});   // signed limits already, do not mirror
         else if (b.r0 <= 0)
            bands.push_back({-b.r1, b.r1});
         else {
            bands.push_back({b.r0, b.r1});
            bands.push_back({-b.r1, -b.r0});
         }
         for (auto &bd : bands) {
            auto *bx = new TBox(b.z0, bd.first, b.z1, bd.second);
            bx->SetFillColorAlpha(col, tel ? 0.85 : 0.25);
            bx->SetLineColor(col); bx->SetLineWidth(tel ? 2 : 1);
            bx->Draw("l");
         }
      }
      // the beam
      auto *beam = new TArrow(zlo, 0, zhi, 0, 0.012, "|>");
      beam->SetLineColor(kRed + 1); beam->SetLineWidth(2); beam->SetLineStyle(2); beam->Draw();
   };

   c->cd(1);
   gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.07); gPad->SetRightMargin(0.02);
   drawFrame(-5, 118, 30, "full setup -- REVERSED: beam enters at the pad plane (z=0), leaves through the cathode (z=100)");
   {
      auto *t = new TLatex(); t->SetTextSize(0.055);
      t->SetTextColor(kGray + 3); t->DrawLatex(30, 22, "AT-TPC drift volume (gas)");
      t->SetTextColor(kAzure + 3); t->DrawLatex(101, -25, "telescope");
      t->SetTextColor(kRed + 1); t->SetTextSize(0.05); t->DrawLatex(-3, 3, "beam");
   }

   c->cd(2);
   gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.07); gPad->SetRightMargin(0.02);
   drawFrame(99.5, 110.5, 7, "zoom -- dE 500 um (green), E 1000 um (blue), CsI array 18x18x25 mm (orange)");

   // --- face-on view -------------------------------------------------------------------------
   // AN x-z SECTION CANNOT SHOW AN ARRAY. The CsI elements tile contiguously in x, so in a side
   // view all 36 collapse into one solid bar and the segmentation -- and whether the array
   // actually covers the silicon -- is invisible. This panel looks down the beam, which is the
   // only projection in which "does the CsI cover the DSSD" is a question you can answer.
   TCanvas *c2 = new TCanvas("cgf", "face on", 620, 620);
   gPad->SetGridx(); gPad->SetGridy(); gPad->SetLeftMargin(0.13);
   TH1F *ff = gPad->DrawFrame(-7, -7, 7, 7);
   ff->SetTitle("looking down the beam;x [cm];y [cm]");
   for (const auto &b : boxes) {
      if (!b.slab) continue;
      const bool csi = b.name.Contains("CsI");
      auto *bx = new TBox(b.xc - b.dx, b.yc - b.dy, b.xc + b.dx, b.yc + b.dy);
      bx->SetFillStyle(0);
      bx->SetLineColor(csi ? kOrange - 3 : (b.name.Contains("_dE") ? kGreen + 2 : kAzure + 2));
      bx->SetLineWidth(csi ? 1 : 3);
      bx->Draw("l");
   }
   {
      auto *t = new TLatex(); t->SetTextSize(0.033);
      t->SetTextColor(kGreen + 2); t->DrawLatex(-6.6, 6.1, "DSSD 10 x 10 cm (dE green, E blue behind)");
      t->SetTextColor(kOrange - 3); t->DrawLatex(-6.6, 5.5, "CsI array, 18 x 18 mm elements");
   }
   TString facePng = outPng; facePng.ReplaceAll(".png", "_faceon.png");
   gSystem->mkdir("plots", kTRUE);
   c2->SaveAs(facePng);

   gSystem->mkdir("plots", kTRUE);
   c->SaveAs(outPng);

   printf("\n  volumes found, in z order:\n");
   std::sort(boxes.begin(), boxes.end(), [](const Box &a, const Box &b) { return a.z0 < b.z0; });
   for (const auto &b : boxes)
      printf("    %-24s z %8.4f .. %8.4f cm   r %6.2f .. %6.2f cm   (thickness %.1f um)\n",
             b.name.Data(), b.z0, b.z1, b.r0, b.r1, (b.z1 - b.z0) * 1e4);
   printf("\n  wrote %s\n\n", outPng.Data());
}
