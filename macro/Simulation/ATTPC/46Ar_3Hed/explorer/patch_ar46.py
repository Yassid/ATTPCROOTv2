#!/usr/bin/env python3
"""Patch a generated 46Ar(3He,d)47K explorer page for this reaction's beam.

    ./patch_ar46.py ~/ar46_3Hed_explorer.html

WHY A PATCH AND NOT AN EDIT TO THE TEMPLATE. explorer_template.html is shared by every channel
(12Be, 15C, a1975 (d,t), ...). All of those have a beam that loses a few MeV crossing the gas, so
one constant E_beam per page is a fair approximation for them. 46Ar is not like that: it enters at
598 MeV and loses 0.957 MeV/cm, i.e. 95.7 MeV over the metre of drift. Editing the shared template
to fix 46Ar would change every other channel's page; patching the generated file changes only this
one. Same reason a1975's add_keoff.py exists.

FOUR THINGS ARE WRONG WITH AN UNPATCHED PAGE, and three of them are silent:

  1. The E_beam slider runs 80-260 MeV. 46Ar needs ~598, so the control cannot even reach the
     right value -- this is the one failure that is obvious on sight.
  2. kine2b() is called with ONE constant beam energy for every event. Every excitation energy is
     then wrong by up to +-48 MeV of beam energy depending on where in the chamber the reaction
     happened, which at this beam energy is several MeV of E_x.
  3. Each kinematic locus is drawn at that same single energy, so it is a thin line through the
     middle of what is really a wide band, and the data looks scattered around a curve that is
     itself only true at one vertex position.
  4. The three 47K levels are merged in the cache with no way to look at one at a time.

Every replacement below asserts it applied EXACTLY once. There is no JS engine on this box to
lint the result, so a silently-missed replacement would ship a page that loads and lies; a
changed template must fail here loudly instead.
"""
import sys
import pathlib


def sub_once(s: str, old: str, new: str, what: str) -> str:
    n = s.count(old)
    if n != 1:
        raise SystemExit(f"PATCH FAILED [{what}]: anchor found {n} times, expected exactly 1.\n"
                         f"The template has changed -- re-read it before editing this script.")
    return s.replace(old, new)


def main(path: str) -> int:
    p = pathlib.Path(path)
    s = p.read_text(encoding="utf-8")

    if "__AR46_PATCHED__" in s:
        print(f"  {p} is already patched -- nothing to do")
        return 0

    # -- 1. the slider has to be able to reach 598 --------------------------------------------
    s = sub_once(s,
                 '<input type="range" id="ebeamR" min="80" max="260" step="0.25">',
                 '<input type="range" id="ebeamR" min="450" max="700" step="0.5">',
                 "ebeam slider range")

    # -- 2. a control to look at one 47K level at a time ---------------------------------------
    s = sub_once(s,
                 '      <button class="primary" id="zeroBtn">solve E<sub>beam</sub> for g.s. = 0</button>\n'
                 '    </div>',
                 '      <button class="primary" id="zeroBtn">solve E<sub>beam</sub> for g.s. = 0</button>\n'
                 '    </div>\n'
                 '\n'
                 '    <div class="grp">\n'
                 '      <div class="eyebrow">47K level (MC truth)</div>\n'
                 '      <div class="row"><label for="stateSel">show</label>\n'
                 '        <select id="stateSel">\n'
                 '          <option value="all">all three</option>\n'
                 '          <option value="0">g.s. 1/2+ (1s1/2)</option>\n'
                 '          <option value="0.36">0.36 3/2+ (0d3/2)</option>\n'
                 '          <option value="2.02">2.02 7/2- (0f7/2)</option>\n'
                 '        </select></div>\n'
                 '      <div class="row"><span class="mono num" style="font-size:11px;color:var(--ink-3)">'
                 'truth, not a measurement</span></div>\n'
                 '    </div>',
                 "level selector control")

    # -- 3. per-event beam energy, and the level cut -------------------------------------------
    s = sub_once(s,
                 "const refColor = i => css(i === 0 ? '--accent' : i === 1 ? '--ref2' : '--refdim');",
                 "const refColor = i => css(i === 0 ? '--accent' : i === 1 ? '--ref2' : '--refdim');\n"
                 "\n"
                 "/* __AR46_PATCHED__ 46Ar loses CFG.dEdz MeV per cm of drift, so the beam energy at the\n"
                 "   reaction is a per-EVENT quantity, not a page-wide constant. vz here is the vertex z\n"
                 "   ALREADY un-mirrored by kin_Ar46.C -- do not flip it again. Falls back to the constant\n"
                 "   when a cache carries no vz, so the page still runs on one. */\n"
                 "const AR46_DEDZ = CFG.dEdz || 0;\n"
                 "const ebAt = (s, d, i) =>\n"
                 "  (AR46_DEDZ && d && d.vz) ? s.ebeam - AR46_DEDZ * d.vz[i] : s.ebeam;\n"
                 "let STATE_SEL = 'all';",
                 "per-event beam energy helper")

    s = sub_once(s,
                 "      && (!HAS_RUN || !d.run || (d.run[i] >= s.runLo && d.run[i] <= s.runHi));\n}",
                 "      && (!HAS_RUN || !d.run || (d.run[i] >= s.runLo && d.run[i] <= s.runHi))\n"
                 "      // MC-truth level cut. `state` is the generated E_x, so this selects on truth and\n"
                 "      // is a diagnostic, never something the experiment could do.\n"
                 "      && (STATE_SEL === 'all' || !d.state || Math.abs(d.state[i] - (+STATE_SEL)) < 0.01);\n}",
                 "level cut in pass()")

    s = sub_once(s,
                 "    const [ex] = kine2b(s.ebeam, thCorr(s, d.th[i], d.ke[i])*Math.PI/180, d.ke[i]);",
                 "    const [ex] = kine2b(ebAt(s, d, i), thCorr(s, d.th[i], d.ke[i])*Math.PI/180, d.ke[i]);",
                 "exOnly beam energy")

    s = sub_once(s,
                 "    const [ex, tcm] = kine2b(s.ebeam, th*Math.PI/180, ke);",
                 "    const [ex, tcm] = kine2b(ebAt(s, d, i), th*Math.PI/180, ke);",
                 "compute beam energy")

    # -- 4. the loci are a BAND, because the vertex position varies -----------------------------
    old_locus = """      const pts = kinLocus(s.ebeam, ex);
      if (!pts.length) continue;
      g.strokeStyle = col; g.lineWidth = 2; g.setLineDash(dash); g.beginPath();
      let started = false;
      for (const [th, ke] of pts) {
        if (th < s.thMin || th > s.thMax || ke < s.keMin || ke > s.keMax) { started = false; continue; }
        const x = PAD.l + (th - s.thMin)/(s.thMax - s.thMin)*pw,
              y = PAD.t + ph - (ke - s.keMin)/(s.keMax - s.keMin)*ph;
        started ? g.lineTo(x, y) : g.moveTo(x, y); started = true;
      }
      g.stroke(); g.setLineDash([]);"""
    new_locus = """      // ONE LOCUS IS A BAND, NOT A LINE. The beam loses AR46_DEDZ MeV/cm, so a reaction at the
      // entrance and one at the pad plane sit on different curves. The solid curve is the middle
      // of the current vertex-z window and the two thin ones are its edges: narrowing the vz
      // slider collapses the band onto the line, which is the fastest way to see how much of the
      // apparent spread is only where the reaction happened rather than detector resolution.
      const zLo = (s.vzLo !== undefined) ? s.vzLo : 0, zHi = (s.vzHi !== undefined) ? s.vzHi : 0;
      const arms = AR46_DEDZ > 0 ? [[zLo, 0.30, 1], [zHi, 0.30, 1], [0.5*(zLo + zHi), 1.0, 2]]
                                 : [[0, 1.0, 2]];
      for (const [zz, alpha, lw] of arms) {
        const pts = kinLocus(s.ebeam - AR46_DEDZ*zz, ex);
        if (!pts.length) continue;
        g.save(); g.globalAlpha = alpha;
        g.strokeStyle = col; g.lineWidth = lw; g.setLineDash(dash); g.beginPath();
        let started = false;
        for (const [th, ke] of pts) {
          if (th < s.thMin || th > s.thMax || ke < s.keMin || ke > s.keMax) { started = false; continue; }
          const x = PAD.l + (th - s.thMin)/(s.thMax - s.thMin)*pw,
                y = PAD.t + ph - (ke - s.keMin)/(s.keMax - s.keMin)*ph;
          started ? g.lineTo(x, y) : g.moveTo(x, y); started = true;
        }
        g.stroke(); g.setLineDash([]); g.restore();
      }"""
    s = sub_once(s, old_locus, new_locus, "kinematic loci as a band")

    # -- 5. wire the selector, inside the main script so scope is not a question ----------------
    s = sub_once(s,
                 "  setEbeam(CFG.ebeam0);\n",
                 "  setEbeam(CFG.ebeam0);\n"
                 "  { const el = document.getElementById('stateSel');\n"
                 "    if (el && !el.dataset.wired) { el.dataset.wired = '1';\n"
                 "      el.addEventListener('change', e => { STATE_SEL = e.target.value; render(); }); } }\n",
                 "level selector wiring")

    # -- 6. axis and cut ranges, framed on THIS reaction ---------------------------------------
    # The template's defaults are built for 12Be/15C recoil PROTONS: the KE-vs-theta map spans
    # theta 0-95 deg and KE 0-40 MeV. 46Ar(3He,d) deuterons come out at 58.9-133.6 deg with up to
    # 57 MeV, so on the stock ranges most of the data sits off the map entirely and the backward
    # half -- where the transfer peaks and where the resolution is best -- is simply not drawn.
    # Framed on the data, with the kinematic envelope for reference: theta_cm 15 deg puts the
    # deuteron at 131.4 deg / 5.1 MeV and theta_cm 80 deg at 59.8 deg / 55.3 MeV.
    for what, old, new in [
        ("map theta range",
         '<input type="number" id="thMin" step="5" value="0"><input type="number" id="thMax" step="5" value="95">',
         '<input type="number" id="thMin" step="5" value="55"><input type="number" id="thMax" step="5" value="140">'),
        ("map KE range",
         '<input type="number" id="keMin" step="1" value="0"><input type="number" id="keMax" step="1" value="40">',
         '<input type="number" id="keMin" step="1" value="0"><input type="number" id="keMax" step="1" value="60">'),
        # The KE CUT keeps everything kinematically possible (max 55.3 MeV) and drops only the
        # pathological tail -- 0.3-1.3% of fits come back with KE up to 19832 MeV, and ~98% of
        # those pass chi2/ndf < 5, so the quality cut cannot remove them. 70 is deliberately
        # permissive: a default cut that removed real data would be worse than showing junk.
        ("KE cut",
         '<input type="number" id="keLo" step="0.5" value="0"><input type="number" id="keHi" step="0.5" value="1000">',
         '<input type="number" id="keLo" step="0.5" value="0"><input type="number" id="keHi" step="0.5" value="70">'),
        # 47K has its three levels inside 0-2.02 MeV; -5..25 wasted three quarters of the axis.
        ("Ex range",
         '<input type="number" id="exLo" step="1" value="-5"><input type="number" id="exHi" step="1" value="25">',
         '<input type="number" id="exLo" step="1" value="-5"><input type="number" id="exHi" step="1" value="10">'),
        ("Ex bins",
         '<input type="number" id="exBins" step="10" value="500">',
         '<input type="number" id="exBins" step="10" value="300">'),
    ]:
        s = sub_once(s, old, new, what)

    # -- 7. a KE vs theta_cm panel ------------------------------------------------------------
    # The template ships KE vs theta_LAB and E_x vs theta_cm, but not KE vs theta_cm -- which is
    # the frame the DWBA is quoted in and the one the proposal's 15-80 deg window refers to.
    #
    # THE TWO theta_cm CONVENTIONS ARE SUPPLEMENTS AND MUST NOT BE MIXED. kine2b() returns
    # tcm = pi - acos(...), i.e. the angle of the RESIDUAL (47K) -- the DWBA convention, and what
    # the E_x vs theta_cm panel already plots. kinLocus() instead sweeps the EJECTILE's CM angle.
    # Drawing one against the other mirrors the locus about 90 deg, which looks entirely plausible
    # on a symmetric-looking distribution. kinLocusCm below returns (180 - tc) so the curve and
    # the points are in the same convention as everything else on the page.
    s = sub_once(s,
                 """    if (ke > 0 && th >= 0) pts.push([th, ke]);
  }
  return pts;
}""",
                 """    if (ke > 0 && th >= 0) pts.push([th, ke]);
  }
  return pts;
}

/** ejectile KE vs theta_cm locus, in the RESIDUAL (DWBA) convention that kine2b returns --
    see the note in patch_ar46.py. 180 - tc, not tc. */
function kinLocusCm(Eb, Ex) {
  const m4 = M4 + Ex, E1 = Eb + M1;
  const P = Math.sqrt(E1*E1 - M1*M1), W = E1 + M2, s = W*W - P*P;
  const pts = [];
  if (s <= (M3 + m4)*(M3 + m4)) return pts;
  const rs = Math.sqrt(s), E3cm = (s + M3*M3 - m4*m4)/(2*rs);
  const p3cm = Math.sqrt(Math.max(0, E3cm*E3cm - M3*M3));
  const beta = P/W, gamma = W/rs;
  for (let tc = 0; tc <= 180; tc += 0.5) {
    const c = Math.cos(tc*Math.PI/180);
    const ke = gamma*(E3cm + beta*p3cm*c) - M3;
    if (ke > 0) pts.push([180 - tc, ke]);
  }
  return pts;
}""",
                 "kinLocusCm")

    # the card, placed next to KE vs theta_lab rather than at the end
    s = sub_once(s,
                 """      <div class="card">
        <h2>E<sub>x</sub> vs &theta;<sub>cm</sub></h2>""",
                 """      <div class="card">
        <h2>ejectile KE vs &theta;<sub>cm</sub></h2>
        <div class="legend"><span>&theta;<sub>cm</sub> of the <b>residual</b> (<sup>47</sup>K) &mdash; the DWBA convention; proposal window 15&ndash;80&deg;</span></div>
        <div class="canvas-wrap"><canvas id="cKC"></canvas><div class="tip" id="tipKC"></div></div>
      </div>

      <div class="card">
        <h2>E<sub>x</sub> vs &theta;<sub>cm</sub></h2>""",
                 "KE vs theta_cm card")

    s = sub_once(s,
                 "  const hEt = new Float64Array(s.cmBins*s.etExBins);",
                 "  const hEt = new Float64Array(s.cmBins*s.etExBins);\n"
                 "  const hKC = new Float64Array(s.cmBins*s.keBins);   // KE vs theta_cm",
                 "hKC allocation")

    s = sub_once(s,
                 "    if (cx >= 0 && cx < s.cmBins && be >= 0 && be < s.etExBins) hEt[be*s.cmBins + cx]++;",
                 "    if (cx >= 0 && cx < s.cmBins && be >= 0 && be < s.etExBins) hEt[be*s.cmBins + cx]++;\n"
                 "    // same theta_cm bin, same KE bin as the theta_lab map, so the two panels are\n"
                 "    // the same events under two frames rather than two different selections\n"
                 "    if (cx >= 0 && cx < s.cmBins && by >= 0 && by < s.keBins) hKC[by*s.cmBins + cx]++;",
                 "hKC fill")

    s = sub_once(s,
                 "  return {n, hEx, hSl, hKT, hKTc, hEt, hVz, vzProf, exW};",
                 "  return {n, hEx, hSl, hKT, hKTc, hEt, hKC, hVz, vzProf, exW};",
                 "hKC return")

    s = sub_once(s,
                 "  views.cEt = drawHeat($('cEt'), r.hEt,",
                 """  views.cKC = drawHeat($('cKC'), r.hKC, s.cmBins, s.keBins, 0, 180, s.keMin, s.keMax,
                       '\\\\theta_{cm} of ^{47}K  (deg)', 'ejectile KE  (MeV)', s, c => {
    if (!s.kinLines) return;
    const {g, pw, ph} = c;
    const zLo = (s.vzLo !== undefined) ? s.vzLo : 0, zHi = (s.vzHi !== undefined) ? s.vzHi : 0;
    const arms = AR46_DEDZ > 0 ? [[zLo, 0.30, 1], [zHi, 0.30, 1], [0.5*(zLo + zHi), 1.0, 2]]
                               : [[0, 1.0, 2]];
    REF_EX.forEach((r2, i) => {
      const col = refColor(i), dash = i < 2 ? [] : [5,4];
      for (const [zz, alpha, lw] of arms) {
        const pts = kinLocusCm(s.ebeam - AR46_DEDZ*zz, r2.ex);
        if (!pts.length) continue;
        g.save(); g.globalAlpha = alpha;
        g.strokeStyle = col; g.lineWidth = lw; g.setLineDash(dash); g.beginPath();
        let started = false;
        for (const [tc, ke] of pts) {
          if (ke < s.keMin || ke > s.keMax) { started = false; continue; }
          const x = PAD.l + tc/180*pw, y = PAD.t + ph - (ke - s.keMin)/(s.keMax - s.keMin)*ph;
          started ? g.lineTo(x, y) : g.moveTo(x, y); started = true;
        }
        g.stroke(); g.setLineDash([]); g.restore();
      }
    });
    // the proposal's own angular window, so the panel is read against it
    g.strokeStyle = css('--accent'); g.lineWidth = 1.5; g.setLineDash([4,4]);
    for (const v of [15, 80]) {
      const x = PAD.l + v/180*pw;
      g.beginPath(); g.moveTo(x, PAD.t); g.lineTo(x, PAD.t + ph); g.stroke();
    }
    g.setLineDash([]);
  }, 'ejectile KE vs \\\\theta_{cm}');

  views.cEt = drawHeat($('cEt'), r.hEt,""",
                 "KE vs theta_cm draw")

    s = sub_once(s,
                 "  hookTip('cKT','tipKT',heatTip); hookTip('cEt','tipEt',heatTip);",
                 "  hookTip('cKT','tipKT',heatTip); hookTip('cEt','tipEt',heatTip);\n"
                 "  hookTip('cKC','tipKC',heatTip);",
                 "KE vs theta_cm tip")

    s = sub_once(s,
                 "const names = {cEx:'Ex', cKT:'KE_vs_thetalab',",
                 "const names = {cEx:'Ex', cKT:'KE_vs_thetalab', cKC:'KE_vs_thetacm',",
                 "KE vs theta_cm png name")

    p.write_text(s, encoding="utf-8")
    print(f"  patched {p}")
    print("    E_beam slider 450-700 MeV; per-event E_beam(z) = ebeam - dEdz*vz;")
    print("    loci drawn as a band across the vertex-z window; 47K level selector added.")
    print("    ranges framed on the data: map theta 55-140 deg, KE 0-60 MeV, Ex -5..10, KE cut 0-70.")
    print("    added a KE vs theta_cm panel (residual/DWBA convention) with the 15-80 deg window marked.")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
