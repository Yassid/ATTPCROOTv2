/// @file omp_fit_result_C14.C
/// @brief Final g.s. angular distribution (theta_cm 20-120) against KD03 and the fitted potential.
/// Both DWBA curves are folded with the measured theta_cm response and acceptance before comparison
/// -- fitting the raw curve would let the optical model absorb the detector resolution.
void omp_fit_result_C14(Double_t cmMax = 120.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString F = here + "/../fresco/outputs/";
   auto load = [](TString p) { auto *g = new TGraph(); std::ifstream in(p.Data()); double a, x; int n = 0;
                               while (in >> a >> x) g->SetPoint(n++, a, x); return g; };
   TGraph *kd = load(F + "p14C_el_161_dsdo.dat"), *bf = load(F + "p14C_el_161_bestfit_dsdo.dat");
   // data + response + acceptance from the exported text (identical inputs to the fit)
   auto *gd = new TGraphErrors(); int nd = 0;
   { std::ifstream in((here + "/plots/ompfit/data.txt").Data()); std::string l;
     while (std::getline(in, l)) { if (l.empty() || l[0] == '#') continue; std::istringstream is(l);
       double a, y, e; if (is >> a >> y >> e) { gd->SetPoint(nd, a, y); gd->SetPointError(nd, 0, e); ++nd; } } }
   std::map<double,double> acc; { std::ifstream in((here + "/plots/ompfit/acceptance.txt").Data());
     double a, v; while (in >> a >> v) acc[a] = v; }
   std::vector<std::array<double,3>> R; { std::ifstream in((here + "/plots/ompfit/response.txt").Data()); std::string l;
     while (std::getline(in, l)) { if (l.empty() || l[0] == '#') continue; std::istringstream is(l);
       double t, r, p; if (is >> t >> r >> p) R.push_back({t, r, p}); } }
   auto foldit = [&](TGraph *g) { std::map<double,double> num;
     for (auto &e : R) { double a = acc.count(e[0]) ? acc[e[0]] : 0; if (a <= 0) continue;
       num[e[1]] += e[2] * a * g->Eval(e[0]) * std::sin(e[0] * TMath::DegToRad()); }
     auto *o = new TGraph(); int n = 0;
     for (auto &kv : num) { double a = acc.count(kv.first) ? acc[kv.first] : 0, s = std::sin(kv.first * TMath::DegToRad());
       if (a > 0.05 && s > 1e-3) o->SetPoint(n++, kv.first, kv.second / a / s); }
     return o; };
   TGraph *kdF = foldit(kd), *bfF = foldit(bf);
   auto norm = [&](TGraph *g) { double sn = 0, sd = 0;
     for (int i = 0; i < gd->GetN(); ++i) { double a = gd->GetX()[i], y = gd->GetY()[i], e = gd->GetEY()[i];
       double f = g->Eval(a); if (f <= 0 || e <= 0) continue; sn += y * f / (e * e); sd += f * f / (e * e); }
     return sd > 0 ? sn / sd : 1.0; };
   double kK = norm(kdF), kB = norm(bfF);
   for (int i = 0; i < kdF->GetN(); ++i) kdF->SetPointY(i, kdF->GetPointY(i) * kK);
   for (int i = 0; i < bfF->GetN(); ++i) bfF->SetPointY(i, bfF->GetPointY(i) * kB);

   TCanvas *c = new TCanvas("c", "omp fit", 1300, 560); c->Divide(2, 1);
   c->cd(1); gPad->SetLogy();
   auto *fr = new TH1D("fr", "^{14}C(p,p) g.s., #theta_{cm} 20-120#circ;#theta_{cm} [deg];d#sigma/d#Omega [arb]", 1, 15, cmMax + 5);
   fr->SetMinimum(50); fr->SetMaximum(1e5); fr->Draw();
   kdF->SetLineColor(kGray + 2); kdF->SetLineWidth(3); kdF->SetLineStyle(2); kdF->Draw("L same");
   bfF->SetLineColor(kRed + 1); bfF->SetLineWidth(3); bfF->Draw("L same");
   gd->SetMarkerStyle(20); gd->SetMarkerColor(kBlack); gd->SetLineColor(kBlack); gd->SetLineWidth(2); gd->Draw("P same");
   auto *lg = new TLegend(0.40, 0.70, 0.89, 0.88);
   lg->AddEntry(gd, "GENFIT, #theta-corrected, acc-corrected", "lp");
   lg->AddEntry(kdF, "KD03 global, folded", "l");
   lg->AddEntry(bfF, "fitted OMP, folded", "l");
   lg->SetTextSize(0.032); lg->Draw();
   c->cd(2);
   auto *rk = new TGraph(), *rb = new TGraph();
   for (int i = 0; i < gd->GetN(); ++i) { double a = gd->GetX()[i], y = gd->GetY()[i];
     rk->SetPoint(i, a, y / kdF->Eval(a)); rb->SetPoint(i, a, y / bfF->Eval(a)); }
   rk->SetTitle("data / folded DWBA;#theta_{cm} [deg];ratio");
   rk->SetMarkerStyle(24); rk->SetMarkerColor(kGray + 2); rk->SetLineColor(kGray + 2); rk->SetLineWidth(2);
   rk->GetYaxis()->SetRangeUser(0.3, 2.0); rk->GetXaxis()->SetLimits(15, cmMax + 5); rk->Draw("ALP");
   rb->SetMarkerStyle(20); rb->SetMarkerColor(kRed + 1); rb->SetLineColor(kRed + 1); rb->SetLineWidth(2); rb->Draw("LP same");
   auto *one = new TLine(15, 1, cmMax + 5, 1); one->SetLineStyle(2); one->SetLineColor(kGray + 2); one->Draw();
   auto *l2 = new TLegend(0.45, 0.75, 0.89, 0.89);
   l2->AddEntry(rk, "KD03 global", "lp"); l2->AddEntry(rb, "fitted OMP", "lp"); l2->Draw();
   c->SaveAs(here + "/plots/omp_fit_result_C14.png");
   printf("wrote plots/omp_fit_result_C14.png\n");
}
