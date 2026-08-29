/// @file dp_qa_C14.C
/// @brief One-page reconstruction QA for a single 14C(d,p)15C sample.
///
///   root -b -q 'dp_qa_C14.C("/mnt/f/a1954_C14dp","gs_s9001_b285_attpc","/mnt/f/a1954_C14dp/sims_b285/gs_s9001_sim.root")'
///
/// Six panels, in the order you would actually check a new production:
///   A  the kinematic plane, with the two-body locus over it -- do the points sit on the curve
///   B  KE_reco/KE_true against the true lab angle -- the direct test that the fit is unbiased,
///      and in particular that backward tracks (theta_lab > 90) behave like forward ones. Before
///      the measurement-order fix in AtGenfitter this panel fell off a cliff past 90 deg.
///   C  the energy resolution and bias per lab-angle slice
///   D  acceptance against theta_cm, from MC truth
///   E  the reconstructed excitation energy
///   F  E_x against theta_cm -- a discrete state must be FLAT here. Any slope is a kinematic
///      systematic (field, drift velocity, energy loss, vertex) and it smears the integrated
///      spectrum in panel E without ever showing up as a bad fit.

#include <algorithm>
#include <vector>

static const double U = 931.49401;
static const double M1 = 14.003242 * U, M2 = 2.0141018 * U, M3 = 1.007825 * U, M4 = 15.0105993 * U;

static double qa_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static bool qa_forward(double m4, double K, double thcm_deg, double &ke, double &thlab)
{
   double E1 = K + M1, s = M1 * M1 + M2 * M2 + 2 * M2 * E1, rs = std::sqrt(s);
   if (rs < M3 + m4) return false;
   double pcm = qa_om2(s, M3 * M3, m4 * m4) / (2 * rs);
   double E3cm = std::sqrt(pcm * pcm + M3 * M3);
   double plab1 = std::sqrt(E1 * E1 - M1 * M1);
   double beta = plab1 / (E1 + M2), gamma = 1.0 / std::sqrt(1 - beta * beta);
   double th = thcm_deg * TMath::DegToRad();
   double pz = gamma * (pcm * std::cos(th) + beta * E3cm), pt = pcm * std::sin(th);
   ke = gamma * (E3cm + beta * pcm * std::cos(th)) - M3;
   thlab = std::atan2(pt, pz) * TMath::RadToDeg();
   return ke > 0;
}
static double qa_q(std::vector<double> v, double p)
{
   if (v.size() < 15) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_qa_C14(TString dir = "/mnt/f/a1954_C14dp", TString tag = "gs_s9001_b285_attpc",
               TString simFile = "", Double_t Ebeam = 155.9, TString outDir = "")
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(255);
   if (outDir.IsNull())
      outDir = TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/plots/";
   gSystem->mkdir(outDir.Data(), kTRUE);

   TFile *fe = TFile::Open(dir + "/exres_" + tag + ".root");
   if (!fe || fe->IsZombie()) { printf("\033[1;31mno exres file for %s\033[0m\n", tag.Data()); return; }
   TTree *t = (TTree *)fe->Get("res");
   double exReco, thTrue, thReco, keTrue, keReco, cmTrue, zTrue, zReco;
   t->SetBranchAddress("exReco", &exReco);
   t->SetBranchAddress("thTrue", &thTrue);
   t->SetBranchAddress("thReco", &thReco);
   t->SetBranchAddress("keTrue", &keTrue);
   t->SetBranchAddress("keReco", &keReco);
   t->SetBranchAddress("cmTrue", &cmTrue);
   t->SetBranchAddress("zTrue", &zTrue);
   t->SetBranchAddress("zReco", &zReco);

   auto *hKin = new TH2D("hKin", "A  reconstructed kinematics;#theta_{lab} [deg];proton KE [MeV]", 180, 0, 180, 240, 0,
                         60);
   auto *hRat = new TH2D("hRat", "B  energy closure;#theta_{lab} true [deg];KE_{reco} / KE_{true}", 90, 0, 180, 120,
                         0.6, 1.4);
   auto *hEx = new TH1D("hEx", "E  excitation energy;E_{x} [MeV];counts", 200, -2, 4);
   auto *hExCm = new TH2D("hExCm", "F  E_{x} vs #theta_{cm};#theta_{cm} true [deg];E_{x} [MeV]", 60, 0, 180, 120, -2, 4);

   const int NS = 12;
   std::vector<double> sl[NS];
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      hKin->Fill(thReco, keReco);
      if (keTrue > 0) hRat->Fill(thTrue, keReco / keTrue);
      hEx->Fill(exReco);
      hExCm->Fill(cmTrue, exReco);
      int b = (int)(thTrue / 15.0);
      if (b >= 0 && b < NS && keTrue > 0) sl[b].push_back(100.0 * (keReco - keTrue) / keTrue);
   }

   auto *cv = new TCanvas("dpqa", "dpqa", 1650, 1000);
   cv->Divide(3, 2);

   cv->cd(1); gPad->SetLogz(); gPad->SetLeftMargin(0.13);
   hKin->Draw("colz");
   auto *loc = new TGraph();
   for (double tcm = 1; tcm <= 179; tcm += 0.5) {
      double ke, thl;
      if (qa_forward(M4, Ebeam, tcm, ke, thl) && ke < 60) loc->SetPoint(loc->GetN(), thl, ke);
   }
   loc->SetLineColor(kRed + 1); loc->SetLineWidth(3); loc->Draw("l same");

   cv->cd(2); gPad->SetLogz(); gPad->SetLeftMargin(0.13);
   hRat->Draw("colz");
   auto *l1 = new TLine(0, 1, 180, 1); l1->SetLineColor(kRed + 1); l1->SetLineWidth(2); l1->Draw();
   auto *l90 = new TLine(90, 0.6, 90, 1.4); l90->SetLineStyle(2); l90->SetLineColor(kGray + 2); l90->Draw();
   auto *tb = new TLatex(0.62, 0.20, "backward #rightarrow"); tb->SetNDC(); tb->SetTextSize(0.04);
   tb->SetTextColor(kGray + 3); tb->Draw();

   cv->cd(3); gPad->SetLeftMargin(0.14); gPad->SetGridy();
   auto *fr3 = new TH1F("fr3", "C  bias and spread;#theta_{lab} true [deg];KE_{reco}-KE_{true} [%]", 12, 0, 180);
   fr3->GetYaxis()->SetRangeUser(-12, 12); fr3->SetLineColor(kWhite); fr3->Draw();
   auto *gB = new TGraphErrors();
   printf("\n  %s\n  %-12s %7s %10s %10s\n", tag.Data(), "theta_lab", "n", "bias [%]", "spread [%]");
   for (int b = 0; b < NS; ++b) {
      if (sl[b].size() < 15) continue;
      double m = qa_q(sl[b], .5), s = (qa_q(sl[b], .75) - qa_q(sl[b], .25)) / 1.349;
      int i = gB->GetN();
      gB->SetPoint(i, 15 * b + 7.5, m);
      gB->SetPointError(i, 7.5, s);
      printf("  %4d-%-7d %7zu %+10.2f %10.2f\n", 15 * b, 15 * (b + 1), sl[b].size(), m, s);
   }
   gB->SetMarkerStyle(20); gB->SetMarkerColor(kAzure + 2); gB->SetLineColor(kAzure + 2);
   gB->SetLineWidth(2); gB->Draw("p same");
   auto *l0 = new TLine(0, 0, 180, 0); l0->SetLineStyle(2); l0->SetLineColor(kGray + 2); l0->Draw();

   cv->cd(4); gPad->SetLeftMargin(0.14); gPad->SetGridy();
   TFile *fa = TFile::Open(dir + "/acceptance_" + tag + ".root");
   if (fa && !fa->IsZombie()) {
      TH1D *hA = nullptr;
      TIter nx(fa->GetListOfKeys());
      while (auto *k = (TKey *)nx()) if (TString(k->GetName()).BeginsWith("hAcc_")) hA = (TH1D *)fa->Get(k->GetName());
      if (hA) {
         hA->SetTitle("D  acceptance;#theta_{cm} [deg];reconstructed / generated");
         hA->GetYaxis()->SetRangeUser(0, 1.05);
         hA->SetLineColor(kAzure + 2); hA->SetLineWidth(2); hA->SetMarkerStyle(20); hA->Draw("e1");
      }
   }

   cv->cd(5); gPad->SetLeftMargin(0.14);
   hEx->SetLineColor(kAzure + 2); hEx->SetLineWidth(2); hEx->Draw("hist");
   double q25 = 0, q50 = 0, q75 = 0;
   {
      std::vector<double> v;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); if (std::fabs(exReco) < 3) v.push_back(exReco); }
      q25 = qa_q(v, .25); q50 = qa_q(v, .5); q75 = qa_q(v, .75);
      auto *tx = new TLatex(0.17, 0.82,
                            TString::Format("median %+.3f, #sigma = %.3f MeV", q50, (q75 - q25) / 1.349));
      tx->SetNDC(); tx->SetTextSize(0.042); tx->Draw();
   }

   cv->cd(6); gPad->SetLogz(); gPad->SetLeftMargin(0.14);
   hExCm->Draw("colz");
   {  // the median per theta_cm slice, which is what shows a slope the 2D cloud hides
      auto *gp = new TGraph();
      for (int b = 1; b <= hExCm->GetNbinsX(); ++b) {
         auto *py = hExCm->ProjectionY("_py", b, b);
         if (py->GetEntries() < 25) continue;
         double qy;
         double half = 0.5;
         py->GetQuantiles(1, &qy, &half);
         gp->SetPoint(gp->GetN(), hExCm->GetXaxis()->GetBinCenter(b), qy);
      }
      gp->SetMarkerStyle(20); gp->SetMarkerColor(kRed + 1); gp->SetLineColor(kRed + 1);
      gp->SetLineWidth(2); gp->Draw("pl same");
      auto *lz = new TLine(0, 0, 180, 0);
      lz->SetLineStyle(2); lz->SetLineColor(kGray + 2); lz->Draw();
   }

   cv->SaveAs(outDir + "dp_qa_" + tag + ".png");
   printf("\n  Ex: median %+.3f  sigma %.3f MeV   (entries %lld)\n", q50, (q75 - q25) / 1.349, t->GetEntries());
   printf("  wrote %sdp_qa_%s.png\n\n", outDir.Data(), tag.Data());
   fe->Close();
}
