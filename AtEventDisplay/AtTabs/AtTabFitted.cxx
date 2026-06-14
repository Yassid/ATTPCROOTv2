#include "AtTabFitted.h"

#include "AtFittedTrack.h"
#include "AtTabInfo.h"
#include "AtTrackingEvent.h"
#include "AtViewerManager.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>
#include <TEveLine.h>
#include <TEveManager.h>
#include <TString.h>

#include <cmath>
#include <iostream>

ClassImp(AtTabFitted);

namespace {

Color_t TrackColor(int i)
{
   static const std::vector<Color_t> kColors = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                                                kOrange + 7, kCyan + 2, kYellow - 2, kViolet + 1};
   return kColors[i % kColors.size()];
}

} // namespace

AtTabFitted::AtTabFitted(const std::string &name, const std::string &branchName)
   : AtTabBase(name), fBranchName(branchName)
{
   if (AtViewerManager::Instance() == nullptr)
      throw "AtViewerManager must be initialized before creating AtTabFitted!";
   fEntry = &AtViewerManager::Instance()->GetCurrentEntry();
   fEntry->Attach(this);
   // Unique Eve container per branch so two tabs (e.g. UKF + genfit) do not collide.
   fEveTracking = std::make_unique<TEveEventManager>(branchName.c_str());
}

AtTabFitted::~AtTabFitted()
{
   if (fEntry != nullptr)
      fEntry->Detach(this);
}

void AtTabFitted::InitTab()
{
   std::cout << " ===== AtTabFitted::Init =====" << std::endl;

   gEve->AddEvent(fEveTracking.get());
   // Pin the vertex marker set so RemoveElements does not free it.
   fVertexSet->SetDestroyOnZeroRefCnt(false);
   fVertexSet->SetMarkerStyle(29);
   fVertexSet->SetMarkerSize(2.0);
   fVertexSet->SetMarkerColor(kRed + 2);
   fEveTracking->AddElement(fVertexSet.get());

   fInitialized = true;

   std::cout << " AtTabFitted::Init : Initialization complete!" << std::endl;
}

void AtTabFitted::Update(DataHandling::AtSubject *sub)
{
   if (sub != fEntry)
      return;
   if (!fInitialized)
      return;
   UpdateFittedElements();
   gEve->Redraw3D(false);
}

void AtTabFitted::UpdateFittedElements()
{
   // Read AtTrackingEvent directly from FairRootManager. We bypass the
   // AtBranch index/name plumbing because the branch-id round-trip
   // (SetBranchName → GetBranchId → GetBranchName) was resolving to the
   // wrong global branch in chained-source setups.
   if (fTrackingArray == nullptr) {
      fTrackingArray = dynamic_cast<TClonesArray *>(
         FairRootManager::Instance()->GetObject(fBranchName.c_str()));
      if (fTrackingArray == nullptr) {
         std::cout << "[AtTabFitted] branch '" << fBranchName << "' not in this file — "
                      "tab inactive\n";
         return;
      }
   }
   AtTrackingEvent *te = nullptr;
   if (fTrackingArray->GetEntries() > 0)
      te = dynamic_cast<AtTrackingEvent *>(fTrackingArray->At(0));
   if (te == nullptr) {
      fEveTracking->RemoveElements();
      fFittedTrackLines.clear();
      fEveTracking->AddElement(fVertexSet.get());
      fVertexSet->Reset();
      return;
   }

   const auto &fits = te->GetFittedTracks();

   // Clear old contents. RemoveElements deletes Eve-owned TEveLines from
   // the previous event; fVertexSet is pinned and detached — re-add it.
   // Stale TEveLine* in fFittedTrackLines are dangling here; we only clear.
   fEveTracking->RemoveElements();
   fFittedTrackLines.clear();
   fEveTracking->AddElement(fVertexSet.get());
   fVertexSet->Reset();
   fFittedTrackLines.reserve(fits.size());

   // Smoothed-trajectory mode: draw the actual fitted polyline (vertex + the points
   // the Kalman filter integrated, GetSmoothedPositions). Fitter-agnostic, so two
   // tabs can overlay UKF and genfit for the same event. The fitters store positions
   // in the lab frame (z_lab = ZPadPlane - z_digi); AtTabMain draws the hits in the
   // digi frame, so map back z_digi = ZPadPlane - z_lab when fZPadPlane is set.
   if (fDrawSmoothed) {
      auto zDisp = [&](double zLab) { return (fZPadPlane != 0.0) ? (fZPadPlane - zLab) : fZSign * zLab; };
      int nDrawn = 0;
      for (size_t i = 0; i < fits.size(); ++i) {
         const auto &ft = *fits[i];
         const auto &sm = ft.GetSmoothedPositions();
         if (sm.empty())
            continue;
         const auto &pV = ft.GetTrackPropertiesStruct().initialPositionXtr;
         const bool hasVertex = pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0;

         auto *line = new TEveLine(TString::Format("%s_Fit_%zu", fBranchName.c_str(), i));
         const Color_t col = (fFixedColor >= 0) ? fFixedColor : TrackColor((int)i);
         line->SetLineColor(col);
         line->SetMainColor(col);
         line->SetLineWidth(fLineWidthCfg);
         line->SetLineStyle(fLineStyle);

         if (hasVertex)
            line->SetNextPoint(pV.X() * 0.1, pV.Y() * 0.1, zDisp(pV.Z()) * 0.1);
         for (const auto &sp : sm)
            line->SetNextPoint(sp.X() * 0.1, sp.Y() * 0.1, zDisp(sp.Z()) * 0.1);

         fEveTracking->AddElement(line);
         fFittedTrackLines.push_back(line);
         if (fDrawVertex && hasVertex)
            fVertexSet->SetNextPoint(pV.X() * 0.1, pV.Y() * 0.1, zDisp(pV.Z()) * 0.1);
         ++nDrawn;
      }
      std::cout << "[AtTabFitted:" << fBranchName << "] drew " << nDrawn << " of " << fits.size()
                << " smoothed trajectories\n";
      return;
   }

   // Analytical helix from the fit's stored parameters:
   //   vertex = initialPositionXtr (mm)
   //   direction = Kinematics.{theta, phi}  (lab frame)
   //   |p| = sqrt(KE·(KE + 2m))             (KE in MeV, m from amu)
   //   sign(q) = ParticleInfo.charge
   //   B = +fBField_T·ẑ                     (settable; PUMA default 4 T)
   // Cyclotron radius: R[mm] = p[GeV/c] / (0.299792458·1e-3·|q|·B[T])
   // Lorentz force F = q v × B for B = +Bẑ:
   //   q>0 → particle curls clockwise viewed from +z  (φ_lab decreases)
   //   q<0 → counterclockwise                          (φ_lab increases)
   // Centripetal direction at vertex (q>0, B>0): rotate v_xy by −90°.
   // We sample the helix in fine arc-length steps so it draws as a smooth
   // curve regardless of how many smoothed clusters the UKF kept. Display
   // is in cm so we divide mm by 10 (matches AtTabMain).
   const double m_per_amu = 0.9314940954;        // GeV/c²
   const double c_const   = 0.299792458e-3;       // GeV/(c·T·mm)
   const double ds_mm     = std::max(1.0, fMaxArc_mm / std::max(1, fHelixSamples));

   int nDrawn = 0;
   for (size_t i = 0; i < fits.size(); ++i) {
      const auto &ft = *fits[i];
      const auto &props = ft.GetTrackPropertiesStruct();
      const auto &kin = ft.GetKinematics();
      const auto &pinfo = ft.GetParticleInfo();
      const auto &pV = props.initialPositionXtr;

      const bool hasVertex = pV.X() != 0 || pV.Y() != 0 || pV.Z() != 0;
      const bool hasKin = (kin.kineticEnergy > 0) && (pinfo.mass > 0)
                          && (pinfo.charge != 0)
                          && std::isfinite(kin.theta) && std::isfinite(kin.phi);
      if (!hasVertex || !hasKin)
         continue;

      const double m_GeV  = pinfo.mass * m_per_amu;
      const double KE_GeV = kin.kineticEnergy * 1e-3;
      const double p_GeV  = std::sqrt(KE_GeV * (KE_GeV + 2.0 * m_GeV));
      const double sinT   = std::sin(kin.theta);
      const double cosT   = std::cos(kin.theta);
      const double cosP   = std::cos(kin.phi);
      const double sinP   = std::sin(kin.phi);
      const double p_perp = p_GeV * sinT;

      if (p_perp <= 0)
         continue;

      const double R_mm  = p_GeV / (c_const * std::abs((double)pinfo.charge) * std::abs(fBField_T));
      const double q_sgn = (pinfo.charge > 0 ? +1.0 : -1.0);
      const double B_sgn = (fBField_T > 0 ? +1.0 : -1.0);
      const double rot   = -q_sgn * B_sgn; // dφ_around_C / d(s_xy/R)
                                           // q>0,B>0 → −1 (clockwise from +z)

      // Center of the circle in (x,y), at distance R perpendicular to v_xy.
      // F_unit at vertex (q>0,B>0) = (sin φ, −cos φ, 0); centripetal points
      // FROM particle TO center, so C = V + R · F_unit.
      const double cx = pV.X() + R_mm * (q_sgn * B_sgn) *  sinP;
      const double cy = pV.Y() - R_mm * (q_sgn * B_sgn) *  cosP;

      // Initial angle around the center: V − C points opposite to F_unit.
      const double u0 = std::atan2(pV.Y() - cy, pV.X() - cx);

      auto *line = new TEveLine(TString::Format("Fit_%zu", i));
      line->SetLineColor(TrackColor((int)i));
      line->SetLineWidth(2);
      line->SetMainColor(TrackColor((int)i));

      double s = 0.0;
      while (s <= fMaxArc_mm) {
         const double s_xy = s * sinT;
         const double u    = u0 + rot * (s_xy / R_mm);
         const double x    = cx + R_mm * std::cos(u);
         const double y    = cy + R_mm * std::sin(u);
         const double z    = pV.Z() + s * cosT;
         line->SetNextPoint(x * 0.1, y * 0.1, fZSign * z * 0.1);
         s += ds_mm;
      }

      fEveTracking->AddElement(line);
      fFittedTrackLines.push_back(line);

      if (fDrawVertex)
         fVertexSet->SetNextPoint(pV.X() * 0.1, pV.Y() * 0.1, fZSign * pV.Z() * 0.1);

      ++nDrawn;
   }
   std::cout << "[AtTabFitted] drew " << nDrawn << " of " << fits.size()
             << " fitted tracks\n";
}
