/// @file omp_fit_result_C14.C
/// @brief Final g.s. angular distribution (theta_cm 20-120) against KD03 and the fitted potential.
/// Both DWBA curves are folded with the measured theta_cm response and acceptance before comparison
/// -- fitting the raw curve would let the optical model absorb the detector resolution.
void omp_fit_result_C14(Double_t cmMax = 120.0, TString fitDir = "plots/ompfit",
                        TString cmpDir = "plots/ompfit06", TString cmpLbl = "same data, |E_{x}|<0.6 (old)",
                        Bool_t doFold = kTRUE, TString outTag = "")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString F = here + "/../fresco/outputs/";
   auto load = [](TString p) { auto *g = new TGraph(); std::ifstream in(p.Data()); double a, x; int n = 0;
                               while (in >> a >> x) g->SetPoint(n++, a, x); return g; };
   TGraph *kd = load(F + "p14C_el_161_dsdo.dat"), *bf = load(F + "p14C_el_161_bestfit_dsdo.dat");
   // data + response + acceptance from the exported text (identical inputs to the fit)
   auto *gd = new TGraphErrors(); int nd = 0;
   { std::ifstream in((here + "/" + fitDir + "/data.txt").Data()); std::string l;
     while (std::getline(in, l)) { if (l.empty() || l[0] == '#') continue; std::istringstream is(l);
       double a, y, e; if (is >> a >> y >> e) { gd->SetPoint(nd, a, y); gd->SetPointError(nd, 0, e); ++nd; } } }
   std::map<double,double> acc; { std::ifstream in((here + "/" + fitDir + "/acceptance.txt").Data());
     double a, v; while (in >> a >> v) acc[a] = v; }
   std::vector<std::array<double,3>> R; { std::ifstream in((here + "/" + fitDir + "/response.txt").Data()); std::string l;
     while (std::getline(in, l)) { if (l.empty() || l[0] == '#') continue; std::istringstream is(l);
       double t, r, p; if (is >> t >> r >> p) R.push_back({t, r, p}); } }
   auto foldit = [&](TGraph *g) { std::map<double,double> num;
     for (auto &e : R) { double a = acc.count(e[0]) ? acc[e[0]] : 0; if (a <= 0) continue;
       num[e[1]] += e[2] * a * g->Eval(e[0]) * std::sin(e[0] * TMath::DegToRad()); }
     auto *o = new TGraph(); int n = 0;
     for (auto &kv : num) { double a = acc.count(kv.first) ? acc[kv.first] : 0, s = std::sin(kv.first * TMath::DegToRad());
       if (a > 0.05 && s > 1e-3) o->SetPoint(n++, kv.first, kv.second / a / s); }
     return o; };
   // doFold = kFALSE compares against the RAW DWBA. Keep in mind what that means: the raw curve
   // has a sharp zero-resolution minimum the apparatus cannot produce, so the data must sit above
   // it there by construction. The fold is not cosmetic -- it is the part of the prediction that
   // belongs to the detector.
   TGraph *kdF = doFold ? foldit(kd) : (TGraph *)kd->Clone("kdR");
   TGraph *bfF = doFold ? foldit(bf) : (TGraph *)bf->Clone("bfR");
   auto *gc = new TGraphErrors(); int nc2 = 0;
   if (!cmpDir.IsNull()) { std::ifstream in((here + "/" + cmpDir + "/data.txt").Data()); std::string l;
     while (std::getline(in, l)) { if (l.empty() || l[0] == '#') continue; std::istringstream is(l);
       double a, y, e; if (is >> a >> y >> e) { gc->SetPoint(nc2, a, y); gc->SetPointError(nc2, 0, e); ++nc2; } } }
   auto norm = [&](TGraph *g) { double sn = 0, sd = 0;
     for (int i = 0; i < gd->GetN(); ++i) { double a = gd->GetX()[i], y = gd->GetY()[i], e = gd->GetEY()[i];
       double f = g->Eval(a); if (f <= 0 || e <= 0) continue; sn += y * f / (e * e); sd += f * f / (e * e); }
     return sd > 0 ? sn / sd : 1.0; };
   double kK = norm(kdF), kB = norm(bfF);
   for (int i = 0; i < kdF->GetN(); ++i) kdF->SetPointY(i, kdF->GetPointY(i) * kK);
   for (int i = 0; i < bfF->GetN(); ++i) bfF->SetPointY(i, bfF->GetPointY(i) * kB);

   TCanvas *c = new TCanvas("c", "omp fit", 1300, 560); c->Divide(2, 1);
   c->cd(1); gPad->SetLogy();
   auto *fr = new TH1D("fr", TString::Format("^{14}C(p,p) g.s., #theta_{cm} 20-120#circ%s;#theta_{cm} [deg];d#sigma/d#Omega [arb]", doFold ? "" : "   (UNFOLDED DWBA)"), 1, 15, cmMax + 5);
   fr->SetMinimum(50); fr->SetMaximum(1e5); fr->Draw();
   kdF->SetLineColor(kGray + 2); kdF->SetLineWidth(3); kdF->SetLineStyle(2); kdF->Draw("L same");
   bfF->SetLineColor(kRed + 1); bfF->SetLineWidth(3); bfF->Draw("L same");
   if (nc2 > 0) { // the comparison set gets its OWN normalisation, else the window change would
                  // masquerade as a scale difference rather than the shape change it is
     double sn = 0, sd = 0;
     for (int i = 0; i < gc->GetN(); ++i) { double a = gc->GetX()[i], y = gc->GetY()[i], e = gc->GetEY()[i];
       double f = kdF->Eval(a); if (f <= 0 || e <= 0) continue; sn += y * f / (e * e); sd += f * f / (e * e); }
     double kc = sd > 0 ? sd / sn : 1.0;
     for (int i = 0; i < gc->GetN(); ++i) { gc->SetPointY(i, gc->GetY()[i] * kc); gc->SetPointError(i, 0, gc->GetEY()[i] * kc); }
     gc->SetMarkerStyle(24); gc->SetMarkerColor(kAzure + 2); gc->SetLineColor(kAzure + 2); gc->Draw("P same"); }
   gd->SetMarkerStyle(20); gd->SetMarkerColor(kBlack); gd->SetLineColor(kBlack); gd->SetLineWidth(2); gd->Draw("P same");
   auto *lg = new TLegend(0.40, 0.66, 0.89, 0.88);
   lg->AddEntry(gd, "GENFIT, #theta-corr, |E_{x}|<2.0", "lp");
   if (nc2 > 0) lg->AddEntry(gc, cmpLbl, "p");
   lg->AddEntry(kdF, doFold ? "KD03 global, folded" : "KD03 global, RAW", "l");
   lg->AddEntry(bfF, doFold ? "fitted OMP, folded" : "fitted OMP, RAW", "l");
   lg->SetTextSize(0.032); lg->Draw();
   c->cd(2);
   auto *rk = new TGraph(), *rb = new TGraph();
   for (int i = 0; i < gd->GetN(); ++i) { double a = gd->GetX()[i], y = gd->GetY()[i];
     rk->SetPoint(i, a, y / kdF->Eval(a)); rb->SetPoint(i, a, y / bfF->Eval(a)); }
   rk->SetTitle(TString::Format("data / %s DWBA;#theta_{cm} [deg];ratio", doFold ? "folded" : "RAW"));
   rk->SetMarkerStyle(24); rk->SetMarkerColor(kGray + 2); rk->SetLineColor(kGray + 2); rk->SetLineWidth(2);
   rk->GetYaxis()->SetRangeUser(0.3, 1.9); rk->GetXaxis()->SetLimits(15, cmMax + 5); rk->Draw("ALP");
   rb->SetMarkerStyle(20); rb->SetMarkerColor(kRed + 1); rb->SetLineColor(kRed + 1); rb->SetLineWidth(2); rb->Draw("LP same");
   auto *one = new TLine(15, 1, cmMax + 5, 1); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
   auto *rc = new TGraph();
   if (nc2 > 0) for (int i = 0; i < gc->GetN(); ++i) { double a = gc->GetX()[i];
     rc->SetPoint(i, a, gc->GetY()[i] / kdF->Eval(a)); }
   rc->SetMarkerStyle(24); rc->SetMarkerColor(kAzure + 2); rc->SetLineColor(kAzure + 2); rc->SetLineStyle(2);
   if (nc2 > 0) rc->Draw("LP same");
   auto *l2 = new TLegend(0.40, 0.70, 0.89, 0.89);
   l2->AddEntry(rk, "KD03 global, |E_{x}|<2.0", "lp");
   if (nc2 > 0) l2->AddEntry(rc, "KD03 global, |E_{x}|<0.6", "lp");
   l2->AddEntry(rb, "fitted OMP, |E_{x}|<2.0", "lp"); l2->SetTextSize(0.032); l2->Draw();
   c->SaveAs(here + "/plots/omp_fit_result_C14" + outTag + ".png");
   printf("data points: %d (new), %d (comparison)\n", gd->GetN(), nc2);
   printf("wrote plots/omp_fit_result_C14%s.png\n", outTag.Data());
}
