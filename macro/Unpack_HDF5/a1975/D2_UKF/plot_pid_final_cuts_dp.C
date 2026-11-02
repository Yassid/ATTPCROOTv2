/// @file plot_pid_final_cuts_dp.C
/// @brief The final Spyral PID (fake-triton veto, vetoed slow deuterons re-added from the circle; plot_pid_spyral_final_dp.C)
///        with two TRACK-QUALITY cuts on tracks that stop in the gas (range_dp.C friend "r"; zfar 10-955 mm) (2026-11-01):
///          MIDPEAK   the charge profile peaks well before the candidate's end and the rest is weak -- the particle
///                    stopped and the track finder appended other hits (67 % of the "short triton-like" tracks)
///          TOO SHORT range shorter than a triton's at the same circle Brho by more than e^-0.4 (0.67x): no Z = 1 or 2
///                    species stops shorter than a triton, so this is species-blind. Triton range curve = sim median
///                    vs circle Brho, pol4 in ln-ln, valid 0.13-0.44 T m (cut applied only there).
///        Top: plane after the cuts, and the tracks the cuts remove. Bottom: ln-mass per slice after (blue) / before (grey).
///        midpeakTritonOnly (default true, 2026-11-02): MIDPEAK only for tracks tLo-tHi above their slice's deuteron mode in
///        ln-mass -- applied everywhere it also removed 17 % of the slow deuterons and a strip of protons.
///   root -l -b -q 'plot_pid_final_cuts_dp.C("plots/pid_final_cuts.png")'
void plot_pid_final_cuts_dp(TString png, double shortCut = -0.4, bool midpeakTritonOnly = true, double tLo = 0.20, double tHi = 0.60, double vetoLn = 0.25, double vetoMaxBrho = 0.45, double vetoHighSq = 30,
                            TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root", TString rr = "/mnt/f/a1975/pid_range/R_all.root",
                            TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&sqp_g>0&&br>0&&bA>0&&s1_g>0")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   TFile f(pv);
   auto *w = (TTree *)f.Get("w");
   w->AddFriend("r=r", rr);
   TString drop = Form("(log(br/bA)>=%g&&(br<%g||sqp_g>%g))", vetoLn, vetoMaxBrho, vetoHighSq);
   // ln of the sim triton median range at the circle Brho (pol4 in ln bA, fitted 2026-11-01 on range_dp.C sims)
   TString lrt = "(21.1479+39.9049*log(bA)+42.3467*pow(log(bA),2)+20.4999*pow(log(bA),3)+3.6219*pow(log(bA),4))";
   TString stop = "(r.ok&&r.cont&&r.zfar>10&&r.zfar<955)";
   w->Draw("sqp_g/s1_g", sel + "&&!" + drop + "&&br>0.12&&br<0.25&&abs(log(br/bA))<0.1", "goff");
   std::vector<double> rs(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(rs.begin(), rs.end());
   const double sc = rs[rs.size() / 2];
   TString B = Form("(%s?bA:br)", drop.Data()), SQ = Form("(%s?s1_g*%g:sqp_g)", drop.Data(), sc);
   w->Draw("sqp_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2], xhi = x[x.size() * 995 / 1000] * 1.1;
   // per-track flags computed in C++ (as a single TTree formula the per-slice deuteron mode exceeds TTreeFormula's
   // operator limit): q.bad = stopping && ((MIDPEAK [&& in the triton region]) || too short); q.mid, q.shrt the parts.
   // Triton region = tLo..tHi above the deuteron mode of the track's 0.05 T m slice of the plotted Brho, in ln-mass.
   const Long64_t N = w->GetEntries();
   w->SetEstimate(N + 1);
   w->Draw(Form("%s:log(%s/%g*%s):(%s):r.midpeak:(bA>0.13&&bA<0.44&&log(r.len)-%s<%g):(%s)", B.Data(), SQ.Data(), md, B.Data(),
               stop.Data(), lrt.Data(), shortCut, sel.Data()), "", "goff");
   std::vector<double> vB(w->GetVal(0), w->GetVal(0) + N), vL(w->GetVal(1), w->GetVal(1) + N), vStop(w->GetVal(2), w->GetVal(2) + N),
      vMid(w->GetVal(3), w->GetVal(3) + N), vShort(w->GetVal(4), w->GetVal(4) + N), vSel(w->GetVal(5), w->GetVal(5) + N);
   std::vector<double> dmode(20, 0);
   for (int b = 1; b < 20; ++b) {
      TH1D hm("hm", "", 600, -4, 2);
      for (Long64_t i = 0; i < N; ++i) if (vSel[i] && vB[i] > 0.05 * b && vB[i] <= 0.05 * (b + 1)) hm.Fill(vL[i]);
      hm.Smooth(2);
      if (hm.GetEntries() > 50) dmode[b] = hm.GetBinCenter(hm.GetMaximumBin());
   }
   auto *q = new TTree("q", "cut flags");
   q->SetDirectory(nullptr);
   int qbad, qmid, qshrt;
   q->Branch("bad", &qbad); q->Branch("mid", &qmid); q->Branch("shrt", &qshrt);
   for (Long64_t i = 0; i < N; ++i) {
      const int b = std::min(19, std::max(0, (int)(vB[i] / 0.05)));
      const double rel = vL[i] - dmode[b];
      const bool inT = rel > tLo && rel < tHi;
      qmid = vStop[i] && vMid[i] && (!midpeakTritonOnly || inT);
      qshrt = vStop[i] && vShort[i] && !qmid;
      qbad = qmid || qshrt;
      q->Fill();
   }
   w->AddFriend(q);
   TString bad = "q.bad", badMid = "q.mid", badShort = "q.shrt";
   const long nAll = w->GetEntries(sel);
   printf("tracks %ld; removed by MIDPEAK %lld, by TOO SHORT %lld (total %lld = %.1f %%)\n", nAll, w->GetEntries(sel + "&&" + badMid),
          w->GetEntries(sel + "&&" + badShort), w->GetEntries(sel + "&&" + bad), 100. * w->GetEntries(sel + "&&" + bad) / nAll);
   auto *c = new TCanvas("c", "c", 2000, 1300);
   auto *top = new TPad("top", "", 0, 0.45, 1, 1); top->Draw();
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.45); bot->Draw();
   top->Divide(2, 1);
   const char *ttl[2] = {midpeakTritonOnly ? "final PID after the cuts (MIDPEAK in the triton region, too short)" : "final PID after the quality cuts (MIDPEAK, too short)", "removed by the quality cuts"};
   TString cuts[2] = {sel + "&&!" + bad, sel + "&&" + bad};
   for (int p = 0; p < 2; ++p) {
      top->cd(p + 1); gPad->SetLogz(); gPad->SetRightMargin(0.12);
      auto *h = new TH2D(Form("p%d", p), Form("%s;#sqrt{dE/dx} (Spyral scale);B#rho [T m]", ttl[p]), 300, 0, xhi, 300, 0, 1.0);
      w->Draw(Form("%s:%s>>p%d", B.Data(), SQ.Data(), p), cuts[p], "goff");
      h->Draw("colz");
   }
   const double bl[] = {0.10, 0.15, 0.25, 0.35, 0.45, 0.60, 0.80};
   bot->Divide(6, 1);
   TString L = Form("log(%s/%g*%s)", SQ.Data(), md, B.Data());
   TH1D hr("hr", "", 600, -4, 2);
   w->Draw(L + ">>hr", sel + Form("&&%s>0.15&&%s<0.25", B.Data(), B.Data()), "goff");
   hr.Smooth(2);
   const double dref = hr.GetBinCenter(hr.GetMaximumBin());
   for (int b = 0; b < 6; ++b) {
      bot->cd(b + 1); gPad->SetLogy(); gPad->SetGridx();
      TString sl = Form("&&%s>%g&&%s<%g", B.Data(), bl[b], B.Data(), bl[b + 1]);
      TH1D tmp(Form("t%d", b), "", 600, -4, 2);
      w->Draw(L + Form(">>t%d", b), sel + sl, "goff");
      tmp.Smooth(2);
      int i0 = tmp.FindBin(dref - 0.3), i1 = tmp.FindBin(dref + 0.3), bd = i0;
      for (int i = i0; i <= i1; ++i) if (tmp.GetBinContent(i) > tmp.GetBinContent(bd)) bd = i;
      const double d = tmp.GetBinCenter(bd);
      TH1D *h[2];
      for (int k = 0; k < 2; ++k) {
         h[k] = new TH1D(Form("h%d%d", b, k), Form("B#rho %.2f-%.2f T m;ln mass - ln m_{d};tracks", bl[b], bl[b + 1]), 160, -1.2, 1.2);
         w->Draw(Form("%s-(%g)>>h%d%d", L.Data(), d, b, k), (k ? sel + "&&!" + bad : sel) + sl, "goff");
         h[k]->SetLineWidth(2); h[k]->SetLineColor(k ? kBlue + 2 : kGray + 1);
      }
      h[0]->SetMinimum(0.8); h[0]->SetMaximum(3 * h[0]->GetMaximum());
      h[0]->Draw("hist"); h[1]->Draw("hist same");
      for (double v : {std::log(938.272 / 1875.613), std::log(2808.921 / 1875.613), std::log(3727.379 / 1875.613)}) {
         auto *l = new TLine(v, 0.8, v, h[0]->GetMaximum()); l->SetLineStyle(3); l->Draw();
      }
      // triton-window counts before/after
      long n0 = w->GetEntries(sel + sl + Form("&&%s-(%g)>0.25&&%s-(%g)<0.55", L.Data(), d, L.Data(), d));
      long n1 = w->GetEntries(sel + "&&!" + bad + sl + Form("&&%s-(%g)>0.25&&%s-(%g)<0.55", L.Data(), d, L.Data(), d));
      long d0 = w->GetEntries(sel + sl + Form("&&abs(%s-(%g))<0.15", L.Data(), d)), d1 = w->GetEntries(sel + "&&!" + bad + sl + Form("&&abs(%s-(%g))<0.15", L.Data(), d));
      printf("  Brho %.2f-%.2f: triton window %ld -> %ld (%.0f %% kept), deuteron peak %ld -> %ld (%.0f %% kept)\n", bl[b], bl[b + 1], n0, n1,
             100. * n1 / std::max(n0, 1L), d0, d1, 100. * d1 / std::max(d0, 1L));
      if (b == 0) {
         auto *lg = new TLegend(0.12, 0.78, 0.7, 0.89);
         lg->AddEntry(h[0], "before cuts", "l"); lg->AddEntry(h[1], "after cuts", "l"); lg->Draw();
      }
   }
   c->SaveAs(png);
}
