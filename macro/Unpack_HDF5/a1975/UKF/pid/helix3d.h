/// @file helix3d.h
/// @brief Bench-macro alias of the framework estimator: the algorithm now lives in
///        AtTools/AtPID/AtHelix3DCore.h (class AtTools::AtHelix3D wraps it). Including the core here
///        keeps every a1975 bench macro on the SAME code as the framework -- edit the core, not this.
///        The last stand-alone macro version is archived at
///        /mnt/f/a1975/pid3d_multi/helix3d_macro_final_before_port.h (2026-10-30).
///        Difference from that version: ties in the spanning tree are broken by hit POSITION instead
///        of input index (exact input-order independence), and the optional line-likeness test uses
///        a closed-form eigenvalue instead of ROOT::Matrix.
#pragma once
#include "../../../../../AtTools/AtPID/AtHelix3DCore.h"
namespace h3d = AtTools::Helix3D;
