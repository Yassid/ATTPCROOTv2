/// @file export_gui_C17.C
/// @brief Dump the campaign's accepted events to JSON for the standalone HTML kinematics viewer.
///
///   root -b -q 'export_gui_C17.C'
///
/// Writes gui/C17_gui_data.json: for each of the four campaign cells, parallel arrays of the
/// quantities the viewer plots. Parallel arrays rather than an array of objects because the key
/// names would otherwise be repeated ~80000 times and triple the file.
///
/// Per accepted event (chi2/ndf < 5):
///     lev   0 = g.s., 1 = 1/2+ 217, 2 = 5/2+ 332
///     th    reconstructed theta_lab [deg]
///     ke    recoil KE vertex-corrected to E0 = 135 MeV [MeV]
///     ex    excitation energy at E_beam(z_reco)  -- the vertex-corrected reconstruction
///     exr   excitation energy at the constant E0 -- the tree's exReco, i.e. NO correction
///     w     FRESCO dsigma/dOmega at this event's theta_cm [mb/sr], the angular weight
///
/// The elastic has no FRESCO calculation, so it is given the 1/2+ SHAPE here and the viewer sets
/// its normalisation from an explicit N_elastic/N_217 control. The 1/2+ : 5/2+ ratio of 2.38 is
/// real and comes out of the weights on its own.
///
/// E_beam(z) is measured from this campaign's own truth, exactly as in inel_summary_C17.C.

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TSystem.h"
#include "TTree.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace
{
double exOmega2(double x, double y, double z)
{
   return std::sqrt(std::max(0., x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z));
}
double exEx(double m1, double m2, double m3, double m4, double K_proj, double thetalabDeg, double K_eject)
{
   const double thetalab = thetalabDeg * TMath::DegToRad();
   const double Et1 = K_proj + m1, Et3 = K_eject + m3;
   const double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   const double uu = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   const double arg = (std::cos(thetalab) * exOmega2(s, m1 * m1, m2 * m2) * exOmega2(uu, m2 * m2, m3 * m3) -
                       (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - uu)) /
                         (2 * m2 * m2) +
                      s + uu - m2 * m2;
   if (arg <= 0)
      return -1e9;
   return std::sqrt(arg) - m4;
}
double exKE(double m1, double m2, double m3, double m4, double Ebeam, double thLabDeg, double ex, double keMax)
{
   const int nS = 400;
   double prevKE = -1, prevF = 0;
   for (int i = 1; i <= nS; ++i) {
      const double ke = keMax * i / nS;
      const double v = exEx(m1, m2, m3, m4, Ebeam, thLabDeg, ke);
      if (v < -1e8) {
         prevKE = -1;
         continue;
      }
      const double fv = v - ex;
      if (prevKE > 0 && prevF * fv <= 0) {
         double lo = prevKE, hi = ke, flo = prevF;
         for (int it = 0; it < 70; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double fm = exEx(m1, m2, m3, m4, Ebeam, thLabDeg, mid) - ex;
            if (flo * fm <= 0)
               hi = mid;
            else {
               lo = mid;
               flo = fm;
            }
         }
         return 0.5 * (lo + hi);
      }
      prevKE = ke;
      prevF = fv;
   }
   return -1;
}
bool exReadFresco(const std::string &path, std::vector<double> &th, std::vector<double> &xs)
{
   std::ifstream in(path);
   if (!in.good())
      return false;
   std::string line;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#' || line[0] == '@')
         continue;
      double a, b;
      if (sscanf(line.c_str(), "%lf %lf", &a, &b) == 2) {
         th.push_back(a);
         xs.push_back(b);
      }
   }
   return th.size() > 10;
}
double exInterp(const std::vector<double> &th, const std::vector<double> &xs, double x)
{
   if (th.empty() || x < th.front() || x > th.back())
      return 0;
   const size_t i = std::min((size_t)(x - th.front()), th.size() - 2);
   const double f = (x - th[i]) / (th[i + 1] - th[i]);
   return xs[i] * (1 - f) + xs[i + 1] * f;
}
} // namespace

void export_gui_C17(TString root = "/media/yassid/Seagate Hub/ATTPC/C17_inel", TString frescoDir = "./fresco/",
                    TString outFile = "./gui/C17_gui_data.json", Double_t chi2Cut = 5.0)
{
   const int nL = 3;
   const double ExGen[nL] = {0.0, 0.217, 0.332};
   const char *stTag[nL] = {"gs", "ex217", "ex332"};
   const char *cfgs[4] = {"pp_b285", "pp_b400", "dd_b285", "dd_b400"};

   std::vector<double> fth[2], fxs[2];
   if (!exReadFresco((frescoDir + "c17pp_217keV.out").Data(), fth[0], fxs[0]) ||
       !exReadFresco((frescoDir + "c17pp_332keV.out").Data(), fth[1], fxs[1])) {
      printf("\033[1;31mno FRESCO tables in %s\033[0m\n", frescoDir.Data());
      return;
   }
   double sig4pi[2] = {0, 0};
   for (int k = 0; k < 2; ++k)
      for (size_t i = 0; i < fth[k].size(); ++i)
         sig4pi[k] += fxs[k][i] * 2 * TMath::Pi() * std::sin(fth[k][i] * TMath::DegToRad()) * (TMath::Pi() / 180.0);

   const double uAmu = 931.49401;
   const double mBeam = 17.0225787 * uAmu;
   const double E0 = 135.0;

   TString dir = gSystem->DirName(outFile);
   gSystem->mkdir(dir, kTRUE);
   std::ofstream out(outFile.Data());
   out << "{\n";
   out << "\"meta\":{\"E0\":" << E0 << ",\"sig217\":" << sig4pi[0] << ",\"sig332\":" << sig4pi[1]
       << ",\"levels\":[0,0.217,0.332]},\n";
   out << "\"cfg\":{\n";

   for (int c = 0; c < 4; ++c) {
      TString cfgTag(cfgs[c]);
      TString chan(cfgTag);
      chan.Remove(2);
      TString btag(cfgTag);
      btag.Remove(0, 3);
      const double mL = (chan == "dd" ? 2.0141018 : 1.007825) * uAmu;
      const double keMax = (chan == "dd") ? 60.0 : 34.0;

      std::vector<int> vlev;
      std::vector<double> vth, vke, vex, vexr, vw;
      double ebz_a = 0, ebz_b = 0;
      bool haveEbz = false;
      double nGen[nL] = {0, 0, 0};

      for (int l = 0; l < nL; ++l) {
         TString dirq = TString("\"") + root + "/" + cfgTag + "\"";
         TString pat = dirq + "/exres_" + chan + "_" + stTag[l] + "_" + btag + "_s*.root";
         TString found = gSystem->GetFromPipe("ls -1 " + pat + " 2>/dev/null | head -1");
         found = found.Strip(TString::kBoth);
         if (found.IsNull()) {
            printf("\033[1;31m  %s %s MISSING\033[0m\n", cfgTag.Data(), stTag[l]);
            continue;
         }
         TString accLog = dirq + "/" + chan + "_" + stTag[l] + "_" + btag + "_s*_acc.log";
         TString ng = gSystem->GetFromPipe("grep -h 'generated reactions' " + accLog +
                                           " 2>/dev/null | head -1 | awk '{print $3}'");
         nGen[l] = atof(ng.Strip(TString::kBoth).Data());

         TFile *f = TFile::Open(found);
         TTree *t = f ? (TTree *)f->Get("res") : nullptr;
         if (!t)
            continue;
         double thT, thR, keT, keR, cmT, c2n, zT, zR, exR;
         t->SetBranchAddress("thTrue", &thT);
         t->SetBranchAddress("thReco", &thR);
         t->SetBranchAddress("keTrue", &keT);
         t->SetBranchAddress("keReco", &keR);
         t->SetBranchAddress("cmTrue", &cmT);
         t->SetBranchAddress("chi2ndf", &c2n);
         t->SetBranchAddress("zTrue", &zT);
         t->SetBranchAddress("zReco", &zR);
         t->SetBranchAddress("exReco", &exR);

         if (!haveEbz) {
            TGraph g;
            for (Long64_t i = 0; i < t->GetEntries(); ++i) {
               t->GetEntry(i);
               double lo = E0 - 40, hi = E0 + 40;
               if ((exEx(mBeam, mL, mL, mBeam, lo, thT, keT) - ExGen[l]) *
                      (exEx(mBeam, mL, mL, mBeam, hi, thT, keT) - ExGen[l]) >
                   0)
                  continue;
               for (int it = 0; it < 60; ++it) {
                  const double mid = 0.5 * (lo + hi);
                  if ((exEx(mBeam, mL, mL, mBeam, lo, thT, keT) - ExGen[l]) *
                         (exEx(mBeam, mL, mL, mBeam, mid, thT, keT) - ExGen[l]) <=
                      0)
                     hi = mid;
                  else
                     lo = mid;
               }
               g.SetPoint(g.GetN(), zT, 0.5 * (lo + hi));
            }
            if (g.GetN() > 100) {
               TF1 lin("lin", "pol1");
               g.Fit(&lin, "QN");
               ebz_a = lin.GetParameter(0);
               ebz_b = lin.GetParameter(1);
               haveEbz = true;
            }
         }

         for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            if (c2n >= chi2Cut || !haveEbz)
               continue;
            const double w = exInterp(fth[l == 0 ? 0 : l - 1], fxs[l == 0 ? 0 : l - 1], cmT);
            if (w <= 0)
               continue;
            const double exH = exEx(mBeam, mL, mL, mBeam, ebz_a + ebz_b * zR, thR, keR);
            if (exH < -1e8)
               continue;
            const double keC = exKE(mBeam, mL, mL, mBeam, E0, thR, exH, keMax);
            if (keC <= 0)
               continue;
            vlev.push_back(l);
            vth.push_back(thR);
            vke.push_back(keC);
            vex.push_back(exH);
            vexr.push_back(exR);
            vw.push_back(w);
         }
         f->Close();
      }

      // acceptance-folded detected fraction, so the viewer can put counts on a real footing
      double detFrac[nL] = {1, 1, 1};
      {
         const double cmLo = 10.0, cmHi = 178.0;
         const double dcosTot = std::cos(cmLo * TMath::DegToRad()) - std::cos(cmHi * TMath::DegToRad());
         for (int l = 0; l < nL; ++l) {
            if (nGen[l] <= 0)
               continue;
            long nacc = 0;
            for (size_t i = 0; i < vlev.size(); ++i)
               if (vlev[i] == l)
                  ++nacc;
            detFrac[l] = nacc / nGen[l] * 1.0;
            (void)dcosTot;
         }
      }

      auto arrI = [&](const char *name, const std::vector<int> &v) {
         out << "\"" << name << "\":[";
         for (size_t i = 0; i < v.size(); ++i)
            out << (i ? "," : "") << v[i];
         out << "]";
      };
      auto arrD = [&](const char *name, const std::vector<double> &v, int nd) {
         out << "\"" << name << "\":[";
         char buf[32];
         for (size_t i = 0; i < v.size(); ++i) {
            snprintf(buf, sizeof(buf), "%.*f", nd, v[i]);
            out << (i ? "," : "") << buf;
         }
         out << "]";
      };

      out << "\"" << cfgTag << "\":{";
      out << "\"ebz\":[" << ebz_a << "," << ebz_b << "],";
      out << "\"detFrac\":[" << detFrac[0] << "," << detFrac[1] << "," << detFrac[2] << "],";
      out << "\"keMax\":" << keMax << ",";
      arrI("lev", vlev);
      out << ",";
      arrD("th", vth, 2);
      out << ",";
      arrD("ke", vke, 3);
      out << ",";
      arrD("ex", vex, 4);
      out << ",";
      arrD("exr", vexr, 4);
      out << ",";
      arrD("w", vw, 4);
      out << "}" << (c < 3 ? ",\n" : "\n");
      printf("  %-8s %7zu events   E_beam(z) = %.3f %+.6f z   detFrac %.3f/%.3f/%.3f\n", cfgTag.Data(), vlev.size(),
             ebz_a, ebz_b, detFrac[0], detFrac[1], detFrac[2]);
   }
   out << "}}\n";
   out.close();
   printf("\n  wrote %s (%.1f MB)\n\n", outFile.Data(),
          gSystem->GetPathInfo(outFile, nullptr, (Long_t *)nullptr, nullptr, nullptr) == 0 ? 0.0 : 0.0);
}
