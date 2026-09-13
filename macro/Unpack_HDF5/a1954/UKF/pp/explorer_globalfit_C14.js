/* ============================================================================================
   GLOBAL FIT PANEL -- a1954 14C(p,p') add-on, injected into the a2091 explorer page by
   open_explorer_catima_C14.sh (just before init()). Reference in ROOT: globalfit_check_C14.C.

   One E_x spectrum (the cuts of the rail + its own theta_cm window) fitted with ONE model:
       14C levels      gaussians, centroid = literature + shift, sigma = s0 + ds*(E - 6.094)
       14N (p,n) blend gaussian at its measured centroid/width (a different reaction, not 14C)
       1n continuum    3-body phase space 13C + p' + n, reconstructed with the TWO-BODY E_x
       background      linear (two non-negative ramps, so it cannot go below zero)
   Only the amplitudes float (plus shift and/or s0 when ticked). POISSON likelihood, never chi2:
   chi2 discards empty bins and biases small amplitudes low (fit_pd_ps.C / fit_angles_ps_C14.C).

   The live phase space is regenerated whenever E_beam or the theta_cm window change and is cut on
   the same two-body theta_cm as the data, so it follows the controls; "phasespace_1n_C14.root"
   uses the shape baked from that macro instead (159.75 MeV, TRUE theta_cm 20-140, sigma 0.15).
   ========================================================================================== */
const GF_HPS = /*@@GF_HPS@@*/null;       // {lo, hi, v:[...]} from plots/phasespace_1n_C14.root

const GF_LEVELS = [                      // order = colour slot; colour follows the level, not its rank
  {E: 6.091, jp: '1<sup>&minus;</sup>', on: true},
  {E: 6.728, jp: '3<sup>&minus;</sup>', on: true},
  {E: 7.012, jp: '2<sup>+</sup>', on: true},
  {E: 7.341, jp: '2<sup>&minus;</sup>', on: true},
  {E: 8.317, jp: '2<sup>+</sup>', on: true},
  {E: 6.589, jp: '0<sup>+</sup>', on: false},
  {E: 6.903, jp: '0<sup>&minus;</sup>', on: false},
  {E: 9.801, jp: '', on: false}];
const GF_SN = 8.176;
const GF_M13C = 13.00335484*U, GF_MN = 1.00866492*U;     // atomic 13C: electrons balance 14C + H
const GF_SQ2PI = Math.sqrt(2*Math.PI);

/* ---- card, styles, controls --------------------------------------------------------------- */
(() => {
  const st = document.createElement('style');
  st.textContent = `
:root{--gf1:#2a78d6;--gf2:#eb6834;--gf3:#1baf7a;--gf4:#eda100;--gf5:#e87ba4;--gf6:#008300;--gf7:#4a3aa7;--gf8:#e34948}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--gf1:#3987e5;--gf2:#d95926;--gf3:#199e70;--gf4:#c98500;--gf5:#d55181;--gf6:#008300;--gf7:#9085e9;--gf8:#e66767}}
:root[data-theme="dark"]{--gf1:#3987e5;--gf2:#d95926;--gf3:#199e70;--gf4:#c98500;--gf5:#d55181;--gf6:#008300;--gf7:#9085e9;--gf8:#e66767}
.gf-ctl{display:flex;flex-wrap:wrap;align-items:center;gap:6px 18px;margin:6px 0 6px;font-size:12px;color:var(--ink-2)}
.gf-ctl .f{display:inline-flex;align-items:center;gap:5px;white-space:nowrap}
.gf-ctl input[type=number]{width:78px;padding:3px 5px}
.gf-ctl select{font:inherit;font-size:12px;padding:3px 4px;color:var(--ink);background:var(--panel-2);border:1px solid var(--line);border-radius:5px}
.gf-ctl .check{font-size:12px}
.gf-ctl .sep{width:1px;align-self:stretch;background:var(--line-2)}
.gf-ctl .lv i{display:inline-block;width:10px;height:10px;border-radius:2px}
.gf-wrap{aspect-ratio:16/8}
@media (max-width:700px){.gf-wrap{aspect-ratio:4/3}}
.gf-sum{display:flex;flex-wrap:wrap;gap:4px 20px;margin:8px 0 4px;font-size:12px;color:var(--ink-2);
  font-family:ui-monospace,SFMono-Regular,"JetBrains Mono","DejaVu Sans Mono",Menlo,monospace;font-variant-numeric:tabular-nums}
.gf-sum b{color:var(--ink);font-weight:600}
.gf-out{overflow-x:auto}
.gf-out table{border-collapse:collapse;width:100%;font-size:12px}
.gf-out th{font-weight:600;text-align:right;color:var(--ink-3);font-size:10px;letter-spacing:.08em;text-transform:uppercase;
  padding:5px 10px;border-bottom:1px solid var(--line);white-space:nowrap}
.gf-out td{text-align:right;padding:3px 10px;border-bottom:1px solid var(--line-2);white-space:nowrap;
  font-family:ui-monospace,SFMono-Regular,"JetBrains Mono","DejaVu Sans Mono",Menlo,monospace;font-variant-numeric:tabular-nums}
.gf-out th:first-child,.gf-out td:first-child{text-align:left;font-family:inherit}
.gf-out td:last-child,.gf-out th:last-child{text-align:left}
.gf-out i.sw{display:inline-block;width:13px;height:3px;border-radius:2px;margin-right:7px;vertical-align:middle}
.gf-flag{color:var(--warn)}
.gf-note{font-size:11.5px;color:var(--ink-3);margin:8px 0 2px;line-height:1.5;max-width:110ch}`;
  document.head.appendChild(st);

  const lv = GF_LEVELS.map((L, i) =>
    `<label class="check lv"><input type="checkbox" id="gfLv${i}"${L.on ? ' checked' : ''}>`
    + `<i style="background:var(--gf${i + 1})"></i>${L.E.toFixed(3)}${L.jp ? ' ' + L.jp : ''}</label>`).join('');
  const card = document.createElement('div');
  card.className = 'card wide'; card.id = 'cardGf';
  card.innerHTML = `
<h2 id="ttlGf">global fit &middot; <sup>14</sup>C levels + <sup>14</sup>N(p,n) blend + 1n phase space + background</h2>
<div class="legend" id="legGf"></div>
<div class="gf-ctl">
  <span class="f">&theta;<sub>cm</sub> <input type="number" id="gfCmLo" step="5" value="20">&ndash;<input type="number" id="gfCmHi" step="5" value="140"> &deg;</span>
  <span class="f">E<sub>x</sub> <input type="number" id="gfLo" step="0.1" value="5.0">&ndash;<input type="number" id="gfHi" step="0.1" value="10.4"> MeV</span>
  <span class="f">bin <input type="number" id="gfBin" step="10" value="50"> keV</span>
  <span class="sep"></span>
  <span class="f">shift <input type="number" id="gfShift" step="0.005" value="-0.011"> MeV <label class="check"><input type="checkbox" id="gfFloatShift"> float</label></span>
  <span class="f">&sigma; = <input type="number" id="gfSig0" step="0.002" value="0.132"> + <input type="number" id="gfDSig" step="0.001" value="0.0123">&middot;(E&minus;6.094) <label class="check"><input type="checkbox" id="gfFloatSig"> float &sigma;<sub>0</sub></label></span>
  <span class="sep"></span>
  <span class="f"><sup>14</sup>N &mu; <input type="number" id="gfMuN" step="0.01" value="9.178"> &sigma; <input type="number" id="gfSgN" step="0.01" value="0.296"></span>
  <span class="f">1n shape <select id="gfPs"><option value="live">live phase space</option><option value="file">phasespace_1n_C14.root</option></select></span>
</div>
<div class="gf-ctl">
  ${lv}
  <span class="sep"></span>
  <label class="check"><input type="checkbox" id="gfUseN14" checked> <sup>14</sup>N blend</label>
  <label class="check"><input type="checkbox" id="gfUsePs" checked> 1n continuum</label>
  <label class="check"><input type="checkbox" id="gfUseBg" checked> background</label>
  <label class="check"><input type="checkbox" id="gfLogY"> log counts</label>
</div>
<div class="canvas-wrap gf-wrap"><canvas id="cGf"></canvas><div class="tip" id="tipGf"></div></div>
<div class="gf-sum" id="gfSum"></div>
<div class="gf-out" id="gfOut"></div>
<p class="gf-note">Positions are fixed at literature + shift and widths follow &sigma;(E); only the amplitudes float
(and the shift / &sigma;<sub>0</sub> when ticked). Poisson likelihood; &chi;&sup2;<sub>&lambda;</sub> is the Baker&ndash;Cousins
likelihood ratio, Pearson &chi;&sup2; is shown alongside as the ROOT macros print it. Areas are RAW counts in the
&theta;<sub>cm</sub> window &mdash; no acceptance &mdash; with Hessian errors at the fitted shift/&sigma;<sub>0</sub>, floored at &radic;N.
Above S<sub>n</sub> the 1n continuum, the <sup>14</sup>N blend and the background are mutually degenerate, so their split is
not a measurement. The live phase space is 13C + p&prime; + n at the current E<sub>beam</sub>, reconstructed with the
two-body E<sub>x</sub>, cut on the same two-body &theta;<sub>cm</sub> as the data and smeared with &sigma;(E).</p>`;
  const vz = document.getElementById('cardVz');
  vz.parentNode.insertBefore(card, vz.nextSibling);

  IDS.push('gfCmLo', 'gfCmHi', 'gfLo', 'gfHi', 'gfBin', 'gfShift', 'gfSig0', 'gfDSig', 'gfMuN', 'gfSgN', 'gfPs');
  FLAGS.push('gfFloatShift', 'gfFloatSig', 'gfUseN14', 'gfUsePs', 'gfUseBg', 'gfLogY',
             ...GF_LEVELS.map((_, i) => 'gfLv' + i));
})();

/* ---- state -------------------------------------------------------------------------------- */
function gfState() {
  const f = id => parseFloat($(id).value), c = id => $(id).checked;
  const g = {cmLo: f('gfCmLo'), cmHi: f('gfCmHi'), lo: f('gfLo'), hi: f('gfHi'), binKeV: f('gfBin'),
             shift: f('gfShift'), sig0: f('gfSig0'), dSig: f('gfDSig'), muN: f('gfMuN'), sgN: f('gfSgN'),
             psSrc: $('gfPs').value, floatShift: c('gfFloatShift'), floatSig: c('gfFloatSig'),
             useN14: c('gfUseN14'), usePs: c('gfUsePs'), useBg: c('gfUseBg'), logY: c('gfLogY'),
             lv: GF_LEVELS.map((_, i) => c('gfLv' + i))};
  return gfSane(g);
}
function gfSane(g) {
  if (!isFinite(g.cmLo)) g.cmLo = 0;
  if (!(g.cmHi > g.cmLo)) { g.cmLo = 0; g.cmHi = 180; }
  if (!(isFinite(g.lo) && g.hi > g.lo)) { g.lo = 5.0; g.hi = 10.4; }
  if (!isFinite(g.binKeV)) g.binKeV = 50;
  g.binKeV = Math.max(5, Math.min(500, g.binKeV));
  g.nb = Math.max(5, Math.min(2000, Math.round((g.hi - g.lo)/(g.binKeV/1000))));
  g.w = (g.hi - g.lo)/g.nb;
  if (!isFinite(g.shift)) g.shift = 0;
  if (!(g.sig0 > 0.005)) g.sig0 = 0.132;
  if (!isFinite(g.dSig)) g.dSig = 0;
  if (!isFinite(g.muN)) g.muN = 9.178;
  if (!(g.sgN > 0.005)) g.sgN = 0.296;
  return g;
}
const gfSigma = (g, sig0, E) => Math.max(0.01, sig0 + g.dSig*(E - 6.094));

/* ---- the spectrum ------------------------------------------------------------------------- */
function gfHist(s, g) {
  const d = DATA(), n = new Float64Array(g.nb);
  for (let i = 0; i < d.ke.length; ++i) {
    if (!pass(d, i, s)) continue;
    const [ex, tcm] = kine2b(s.ebeam, thCorr(s, d.th[i], d.ke[i])*Math.PI/180, d.ke[i]);
    if (!isFinite(ex) || tcm < g.cmLo || tcm > g.cmHi) continue;
    const b = Math.floor((ex - g.lo)/g.w);
    if (b >= 0 && b < g.nb) n[b]++;
  }
  return n;
}

/* ---- 1n phase space: GENBOD for three bodies (what TGenPhaseSpace does), deterministic seed --- */
const GF_PSF = {lo: 3, hi: 17, w: 0.01};
const gfPsCache = {key: null, fine: null};
function gfRng(seed) {
  let a = seed >>> 0;
  return () => { a = (a + 0x6D2B79F5) >>> 0; let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0)/4294967296; };
}
const gfPdk = (a, b, c) => { const x = (a - b - c)*(a + b + c)*(a - b + c)*(a + b - c); return x > 0 ? Math.sqrt(x)/(2*a) : 0; };
function gfPsFine(eb, cmLo, cmHi, trueCm) {
  const key = [eb, cmLo, cmHi, !!trueCm].join('|');
  if (gfPsCache.key === key) return gfPsCache.fine;
  const nF = Math.round((GF_PSF.hi - GF_PSF.lo)/GF_PSF.w), fine = new Float64Array(nF);
  const mC = GF_M13C, mp = M3, mn = GF_MN;
  const E1 = eb + M1, P = Math.sqrt(E1*E1 - M1*M1), W = E1 + M2, Mw = Math.sqrt(W*W - P*P);
  const bz = P/W, gz = W/Mw;
  if (Mw > mC + mp + mn) {
    const rnd = gfRng(20260826), NEV = 200000;
    const iso = () => { const c = 2*rnd() - 1, ph = 2*Math.PI*rnd(), sn = Math.sqrt(1 - c*c);
                        return [sn*Math.cos(ph), sn*Math.sin(ph), c]; };
    for (let k = 0; k < NEV; ++k) {
      const m01 = mC + mp + rnd()*(Mw - mC - mp - mn);           // invariant mass of 13C + p'
      const p2 = gfPdk(Mw, m01, mn), q = gfPdk(m01, mC, mp), wt = p2*q;
      if (!(wt > 0)) continue;
      const u = iso(), E01 = Math.sqrt(p2*p2 + m01*m01);
      const bx = u[0]*p2/E01, by = u[1]*p2/E01, bzz = u[2]*p2/E01;
      const v = iso(), Ep = Math.sqrt(q*q + mp*mp);
      let px = v[0]*q, py = v[1]*q, pz = v[2]*q;
      const b2 = bx*bx + by*by + bzz*bzz, g1 = 1/Math.sqrt(1 - b2), bp = bx*px + by*py + bzz*pz;
      const ff = b2 > 0 ? (g1 - 1)*bp/b2 + g1*Ep : 0, Ec = g1*(Ep + bp);
      px += ff*bx; py += ff*by; pz += ff*bzz;                     // proton in the overall CM
      const pT = Math.hypot(px, py);
      const El = gz*(Ec + bz*pz), pzl = gz*(pz + bz*Ec);           // and in the lab
      const [ex, tcm] = kine2b(eb, Math.atan2(pT, pzl), El - mp);
      if (!isFinite(ex)) continue;
      const tc = trueCm ? Math.atan2(pT, pz)*180/Math.PI : tcm;
      if (tc < cmLo || tc > cmHi) continue;
      const b = Math.floor((ex - GF_PSF.lo)/GF_PSF.w);
      if (b >= 0 && b < nF) fine[b] += wt;
    }
  }
  gfPsCache.key = key; gfPsCache.fine = fine;
  return fine;
}
/** density of the smeared fine histogram at x (arbitrary normalisation) */
function gfPsSmeared(fine, x, sigAt) {
  let s = 0;
  for (let j = 0; j < fine.length; ++j) {
    const fj = fine[j];
    if (!fj) continue;
    const xj = GF_PSF.lo + (j + 0.5)*GF_PSF.w, sg = sigAt(xj), z = (x - xj)/sg;
    if (z > 6 || z < -6) continue;
    s += fj*Math.exp(-0.5*z*z)/sg;
  }
  return s;
}
/** TH1::Interpolate on the baked hPS */
function gfHpsAt(x) {
  const H = GF_HPS, nb = H.v.length, bw = (H.hi - H.lo)/nb, c = i => H.lo + (i + 0.5)*bw;
  if (x <= c(0)) return H.v[0];
  if (x >= c(nb - 1)) return H.v[nb - 1];
  let i = Math.floor((x - H.lo)/bw);
  if (x < c(i)) i -= 1;
  const t = (x - c(i))/bw;
  return H.v[i] + t*(H.v[i + 1] - H.v[i]);
}

/* ---- amplitude fit: Poisson ML with non-negative amplitudes --------------------------------
   EM (monotone, keeps every amplitude >= 0) to get close, then projected Newton on the free set
   to converge. The model is linear in the amplitudes, so this IS the Minuit "L" minimum with the
   SetParLimits(0, ...) bounds the macros use.                                                  */
function gfModelVec(F, a, m) {
  m.fill(0);
  for (let k = 0; k < F.length; ++k) { const ak = a[k], f = F[k]; if (ak) for (let b = 0; b < m.length; ++b) m[b] += ak*f[b]; }
  for (let b = 0; b < m.length; ++b) if (m[b] < 1e-12) m[b] = 1e-12;
  return m;
}
function gfNll(n, m) { let L = 0; for (let b = 0; b < n.length; ++b) L += m[b] - (n[b] > 0 ? n[b]*Math.log(m[b]) : 0); return L; }
function gfFitAmps(n, F, warm) {
  const K = F.length, nb = n.length;
  const S = F.map(f => { let s = 0; for (const v of f) s += v; return s; });
  let N = 0; for (const v of n) N += v;
  let a = S.map((sk, k) => { const cold = Math.max(1e-6, N/K/Math.max(sk, 1e-12));
                             return warm && warm.length === K ? Math.max(warm[k], 0.02*cold) : cold; });
  const m = new Float64Array(nb);
  let Lprev = Infinity;
  for (let it = 0; it < 400; ++it) {
    gfModelVec(F, a, m);
    for (let k = 0; k < K; ++k) {
      if (!(S[k] > 0)) { a[k] = 0; continue; }
      let num = 0; const f = F[k];
      for (let b = 0; b < nb; ++b) if (n[b]) num += f[b]*n[b]/m[b];
      a[k] *= num/S[k];
    }
    if (it % 20 === 19) { const L = gfNll(n, gfModelVec(F, a, m)); if (Lprev - L < 1e-7) break; Lprev = L; }
  }
  const mt = new Float64Array(nb);
  for (let it = 0; it < 100; ++it) {
    gfModelVec(F, a, m);
    const L0 = gfNll(n, m), grad = new Array(K).fill(0);
    for (let k = 0; k < K; ++k) { const f = F[k]; for (let b = 0; b < nb; ++b) grad[k] += f[b]*(1 - n[b]/m[b]); }
    const free = [];
    for (let k = 0; k < K; ++k) if (S[k] > 0 && (a[k] > 0 || grad[k] < 0)) free.push(k);
    if (!free.length) break;
    const H = free.map(j => free.map(k => { let s = 0; const fj = F[j], fk = F[k];
      for (let b = 0; b < nb; ++b) if (n[b]) s += fj[b]*fk[b]*n[b]/(m[b]*m[b]); return s; }));
    free.forEach((_, i) => { H[i][i] = H[i][i]*(1 + 1e-10) + 1e-12; });
    const step = solve(H, free.map(j => -grad[j]));
    if (!step) break;
    let t = 1, aN = null, Ln = L0;
    for (let ls = 0; ls < 30; ++ls) {
      const trial = a.slice();
      free.forEach((j, i) => { trial[j] = Math.max(0, a[j] + t*step[i]); });
      const L = gfNll(n, gfModelVec(F, trial, mt));
      if (L <= L0 + 1e-12) { aN = trial; Ln = L; break; }
      t *= 0.5;
    }
    if (!aN) break;
    a = aN;
    if (L0 - Ln < 1e-10) break;
  }
  gfModelVec(F, a, m);
  // covariance over the components off their bound (observed information, as HESSE)
  const act = [];
  for (let k = 0; k < K; ++k) if (a[k]*S[k] > 1e-4) act.push(k);
  let cov = null;
  if (act.length) {
    const H = act.map(j => act.map(k => { let s = 0; const fj = F[j], fk = F[k];
      for (let b = 0; b < nb; ++b) if (n[b]) s += fj[b]*fk[b]*n[b]/(m[b]*m[b]); return s; }));
    const inv = act.map(() => new Array(act.length).fill(0));
    let ok = true;
    for (let c = 0; c < act.length && ok; ++c) {
      const col = solve(H.map(r => r.slice()), act.map((_, i) => i === c ? 1 : 0));
      if (!col) { ok = false; break; }
      col.forEach((v, r) => { inv[r][c] = v; });
    }
    if (ok) { cov = {act, inv}; }
  }
  return {a, m, nll: gfNll(n, m), cov, S};
}

/* ---- one complete fit: build shapes, optionally profile shift / s0, collect results -------- */
function gfRun(s, g, opt) {
  opt = opt || {};
  const t0 = performance.now();
  const n = gfHist(s, g), nb = g.nb;
  const x = Float64Array.from({length: nb}, (_, b) => g.lo + (b + 0.5)*g.w);
  let N = 0; for (const v of n) N += v;
  const psLive = g.usePs && g.psSrc === 'live';
  const psFile = g.usePs && g.psSrc === 'file' && GF_HPS;
  const fine = psLive ? gfPsFine(s.ebeam, g.cmLo, g.cmHi, !!opt.trueCm) : null;

  const build = (shift, sig0) => {
    const C = [];
    GF_LEVELS.forEach((L, i) => {
      if (!g.lv[i]) return;
      const mu = L.E + shift, sg = gfSigma(g, sig0, L.E);
      const f = Float64Array.from(x, xx => Math.exp(-0.5*((xx - mu)/sg)**2));
      C.push({kind: 'level', i, E: L.E, mu, sg, f, at: xx => Math.exp(-0.5*((xx - mu)/sg)**2),
              label: L.E.toFixed(3) + (L.jp ? ' ' + L.jp : ''), color: `--gf${i + 1}`, per: sg*GF_SQ2PI/g.w});
    });
    if (g.useN14) {
      const mu = g.muN, sg = g.sgN;
      C.push({kind: 'n14', mu, sg, f: Float64Array.from(x, xx => Math.exp(-0.5*((xx - mu)/sg)**2)),
              at: xx => Math.exp(-0.5*((xx - mu)/sg)**2), label: '<sup>14</sup>N(p,n) blend', color: '--ink-2',
              per: sg*GF_SQ2PI/g.w});
    }
    if (psLive || psFile) {
      const sigAt = opt.psSig ? (() => opt.psSig) : (xx => gfSigma(g, sig0, xx));
      const raw = psLive ? (xx => gfPsSmeared(fine, xx, sigAt)) : gfHpsAt;
      const f = Float64Array.from(x, raw);
      let sum = 0; for (const v of f) sum += v;
      if (sum > 0) {
        for (let b = 0; b < nb; ++b) f[b] /= sum;
        C.push({kind: 'ps', f, at: xx => raw(xx)/sum, label: '1n phase space', color: '--refdim', per: 1});
      }
    }
    if (g.useBg) {
      const r0 = xx => (g.hi - xx)/(g.hi - g.lo), r1 = xx => (xx - g.lo)/(g.hi - g.lo);
      C.push({kind: 'bg0', f: Float64Array.from(x, r0), at: r0, label: 'background', color: '--ink-3', per: 0});
      C.push({kind: 'bg1', f: Float64Array.from(x, r1), at: r1, label: 'background', color: '--ink-3', per: 0});
    }
    return C.filter(c => { let sm = 0; for (const v of c.f) sm += v; return sm > 1e-9; });
  };

  let shift = g.shift, sig0 = g.sig0, warm = null, nEval = 0;
  const prof = (sh, sg) => { const r = gfFitAmps(n, build(sh, sg).map(c => c.f), warm); warm = r.a; ++nEval; return r.nll; };
  const golden = (fn, a, b, tol) => {
    const gr = (Math.sqrt(5) - 1)/2;
    let c = b - gr*(b - a), d = a + gr*(b - a), fc = fn(c), fd = fn(d);
    while (b - a > tol) {
      if (fc < fd) { b = d; d = c; fd = fc; c = b - gr*(b - a); fc = fn(c); }
      else { a = c; c = d; fc = fd; d = a + gr*(b - a); fd = fn(d); }
    }
    return 0.5*(a + b);
  };
  const scan = (fn, lo, hi, step) => {
    let best = lo, bv = Infinity;
    for (let v = lo; v <= hi + 1e-9; v += step) { const L = fn(v); if (L < bv) { bv = L; best = v; } }
    return best;
  };
  if (N > 0 && (g.floatShift || g.floatSig)) {
    for (let round = 0; round < 2; ++round) {
      if (g.floatShift) {
        const c0 = round ? shift : scan(v => prof(v, sig0), g.shift - 0.3, g.shift + 0.3, 0.02);
        shift = golden(v => prof(v, sig0), c0 - 0.025, c0 + 0.025, 4e-4);
      }
      if (g.floatSig) {
        const c0 = round ? sig0 : scan(v => prof(shift, v), 0.04, 0.40, 0.02);
        sig0 = golden(v => prof(shift, v), Math.max(0.015, c0 - 0.025), c0 + 0.025, 3e-4);
      }
      if (!(g.floatShift && g.floatSig)) break;
    }
  }
  const C = build(shift, sig0);
  const r = gfFitAmps(n, C.map(c => c.f), warm);
  // profile-likelihood curvature for the nonlinear parameters (delta -lnL = 1/2)
  const curv = (fn, v, d) => { const L0 = r.nll, Lp = fn(v + d), Lm = fn(v - d), c2 = Lp + Lm - 2*L0;
                               return c2 > 0 ? d/Math.sqrt(c2) : NaN; };
  const shiftErr = g.floatShift && N > 0 ? curv(v => prof(v, sig0), shift, 0.01) : NaN;
  const sig0Err = g.floatSig && N > 0 ? curv(v => prof(shift, v), sig0, 0.005) : NaN;

  // results per component; background's two ramps are reported as one line
  const idxAct = new Map(); if (r.cov) r.cov.act.forEach((k, i) => idxAct.set(k, i));
  const varOf = ks => { if (!r.cov) return NaN; let v = 0;
    for (const [j, wj] of ks) for (const [k, wk] of ks) {
      if (!idxAct.has(j) || !idxAct.has(k)) continue;
      v += wj*wk*r.cov.inv[idxAct.get(j)][idxAct.get(k)]; }
    return v; };
  C.forEach((c, k) => { c.a = r.a[k]; });
  const rows = [];
  C.forEach((c, k) => {
    if (c.kind === 'bg1') return;
    let counts, err;
    if (c.kind === 'bg0') {
      const k1 = C.findIndex(o => o.kind === 'bg1');
      counts = r.a[k]*r.S[k] + (k1 >= 0 ? r.a[k1]*r.S[k1] : 0);
      err = Math.sqrt(Math.max(0, varOf([[k, r.S[k]], ...(k1 >= 0 ? [[k1, r.S[k1]]] : [])])));
    } else {
      const per = c.kind === 'ps' ? r.S[k] : c.per;
      counts = r.a[k]*per;
      err = Math.sqrt(Math.max(0, varOf([[k, per]])));
    }
    const railed = counts < 0.5;
    const floored = counts > 0 && !(err >= Math.sqrt(counts));
    if (floored) err = Math.sqrt(counts);
    rows.push({c, counts, err, railed, floored});
  });

  let pear = 0, bc = 0, used = 0;
  for (let b = 0; b < nb; ++b) {
    const m = r.m[b], o = n[b];
    if (m > 0) { pear += (o - m)*(o - m)/m; ++used; }
    bc += 2*(m - o + (o > 0 ? o*Math.log(o/m) : 0));
  }
  const npar = C.length + (g.floatShift ? 1 : 0) + (g.floatSig ? 1 : 0);
  const ndf = used - npar;
  let modelSum = 0; for (const v of r.m) modelSum += v;
  return {g, n, x, N, C, rows, m: r.m, nll: r.nll, bc, pear, ndf, npar, shift, sig0, shiftErr, sig0Err,
          modelSum, covOk: !!r.cov, nEval, ms: performance.now() - t0,
          psNote: psLive ? `live, E_beam ${s.ebeam.toFixed(2)} MeV` : psFile ? 'phasespace_1n_C14.root' : 'off'};
}

/* ---- drawing ------------------------------------------------------------------------------ */
function gfDraw(cv, R) {
  const c = setup(cv); const {g: ctx, w, h, pw, ph} = c, G = R.g;
  const phMain = Math.round(ph*0.73), gap = 22, yP0 = PAD.t + phMain + gap, phP = ph - phMain - gap;
  const X = v => PAD.l + (v - G.lo)/(G.hi - G.lo)*pw;
  let max = 0;
  for (let b = 0; b < G.nb; ++b) max = Math.max(max, R.n[b] + Math.sqrt(R.n[b]), R.m[b]);
  const logY = G.logY && max > 0;
  const ylo = logY ? 0.5 : 0, yhi = logY ? max*1.8 : (max*1.12 || 1);
  const Y = v => {
    if (!logY) return PAD.t + phMain - (v - ylo)/(yhi - ylo)*phMain;
    return PAD.t + phMain - (Math.log10(Math.max(v, 0.3)) - Math.log10(ylo))/(Math.log10(yhi) - Math.log10(ylo))*phMain;
  };
  const cm = Object.assign({}, c, {ph: phMain});
  const title = `${nucl(CFG.tag || '')}   ${G.cmLo}\u00B0 < \\theta_{cm} < ${G.cmHi}\u00B0`;
  axes(cm, G.lo, G.hi, ylo, yhi, '', `counts / ${Math.round(G.w*1000)} keV`, logY, title);

  ctx.save();
  ctx.beginPath(); ctx.rect(PAD.l, PAD.t, pw, phMain); ctx.clip();
  // S_n
  if (GF_SN > G.lo && GF_SN < G.hi) {
    ctx.strokeStyle = css('--accent'); ctx.lineWidth = 1; ctx.setLineDash([3, 4]);
    ctx.beginPath(); ctx.moveTo(X(GF_SN), PAD.t); ctx.lineTo(X(GF_SN), PAD.t + phMain); ctx.stroke(); ctx.setLineDash([]);
  }
  const NPX = Math.max(200, Math.round(pw));
  const curve = (fn, color, width, dash, fill) => {
    const pts = [];
    for (let p = 0; p <= NPX; ++p) { const xx = G.lo + p/NPX*(G.hi - G.lo); pts.push([X(xx), Y(fn(xx))]); }
    if (fill) {
      ctx.beginPath(); ctx.moveTo(pts[0][0], Y(ylo));
      for (const [px, py] of pts) ctx.lineTo(px, py);
      ctx.lineTo(pts[pts.length - 1][0], Y(ylo)); ctx.closePath();
      ctx.fillStyle = css(color); ctx.globalAlpha = 0.16; ctx.fill(); ctx.globalAlpha = 1;
    }
    ctx.strokeStyle = css(color); ctx.lineWidth = width; ctx.setLineDash(dash || []);
    ctx.beginPath(); pts.forEach(([px, py], i) => i ? ctx.lineTo(px, py) : ctx.moveTo(px, py)); ctx.stroke();
    ctx.setLineDash([]);
  };
  const ps = R.C.find(o => o.kind === 'ps'), b0 = R.C.find(o => o.kind === 'bg0'), b1 = R.C.find(o => o.kind === 'bg1');
  const bgAt = xx => (b0 ? b0.a*b0.at(xx) : 0) + (b1 ? b1.a*b1.at(xx) : 0);
  if (b0 || b1) curve(bgAt, '--ink-3', 1.5, [2, 3]);
  if (ps) curve(xx => ps.a*ps.at(xx), '--refdim', 1.5, null, true);
  for (const o of R.C) if (o.kind === 'n14') curve(xx => o.a*o.at(xx), '--ink-2', 1.5, [6, 3]);
  for (const o of R.C) if (o.kind === 'level' && o.a > 0) curve(xx => o.a*o.at(xx), o.color, 1.6, [5, 3]);
  curve(xx => R.C.reduce((sm, o) => sm + o.a*o.at(xx), 0), '--fit', 2.25);
  // data: points with sqrt(N) bars
  ctx.strokeStyle = css('--ink'); ctx.fillStyle = css('--ink'); ctx.lineWidth = 1;
  for (let b = 0; b < G.nb; ++b) {
    const v = R.n[b]; if (!(v > 0)) continue;
    const px = X(R.x[b]), e = Math.sqrt(v);
    ctx.beginPath(); ctx.moveTo(px, Y(Math.max(v - e, logY ? 0.3 : 0))); ctx.lineTo(px, Y(v + e)); ctx.stroke();
    ctx.beginPath(); ctx.arc(px, Y(v), 2.2, 0, 2*Math.PI); ctx.fill();
  }
  ctx.restore();
  if (GF_SN > G.lo && GF_SN < G.hi)
    drawRich(ctx, 'S_n', X(GF_SN) + 4, PAD.t + 14, {size: 11, color: css('--accent-ink')});
  // fit quality, top right, so an exported figure carries it
  const q = [`N = ${R.N}`, `\\chi^2_{\\lambda}/ndf = ${(R.bc/R.ndf).toFixed(2)}`, `Pearson \\chi^2/ndf = ${(R.pear/R.ndf).toFixed(2)}`];
  q.forEach((t, i) => drawRich(ctx, t, PAD.l + pw - 8, PAD.t + 16 + i*15, {size: 11.5, align: 'right', color: css('--ink-2')}));

  // pull strip
  ctx.fillStyle = css('--plot'); ctx.fillRect(PAD.l - 60, PAD.t + phMain + 1, pw + 70, gap + phP + 2);
  let pm = 3;
  const pulls = Array.from(R.n, (v, b) => R.m[b] > 0 ? (v - R.m[b])/Math.sqrt(R.m[b]) : 0);
  for (const p of pulls) pm = Math.max(pm, Math.min(8, Math.ceil(Math.abs(p))));
  const YP = v => yP0 + phP/2 - v/pm*(phP/2);
  ctx.font = '10px ui-monospace,SFMono-Regular,"DejaVu Sans Mono",monospace';
  ctx.strokeStyle = css('--grid'); ctx.lineWidth = 1; ctx.fillStyle = css('--ink-3');
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (const t of [-pm, -2, 0, 2, pm]) {
    ctx.setLineDash(t === 0 ? [] : [2, 3]);
    ctx.beginPath(); ctx.moveTo(PAD.l, YP(t)); ctx.lineTo(PAD.l + pw, YP(t)); ctx.stroke();
    ctx.fillText(String(t), PAD.l - 6, YP(t));
  }
  ctx.setLineDash([]);
  ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (const t of ticks(G.lo, G.hi, 7)) {
    const px = X(t);
    ctx.beginPath(); ctx.moveTo(px, yP0); ctx.lineTo(px, yP0 + phP); ctx.stroke();
    ctx.fillText(String(t), px, yP0 + phP + 6);
  }
  const bwpx = Math.max(1, pw/G.nb - 2);
  ctx.fillStyle = css('--data');
  pulls.forEach((p, b) => {
    const pc = Math.max(-pm, Math.min(pm, p)), px = X(R.x[b]);
    ctx.fillRect(px - bwpx/2, Math.min(YP(0), YP(pc)), bwpx, Math.abs(YP(pc) - YP(0)));
  });
  ctx.strokeStyle = css('--line'); ctx.lineWidth = 1.2; ctx.strokeRect(PAD.l, yP0, pw, phP);
  ctx.save(); ctx.translate(16, yP0 + phP/2); ctx.rotate(-Math.PI/2);
  drawRich(ctx, 'pull', 0, 0, {size: 12, align: 'center', color: css('--ink-2')});
  ctx.restore();
  drawRich(ctx, 'E_x  (MeV)', PAD.l + pw/2, h - 12, {size: 13, align: 'center', color: css('--ink-2')});

  // CSV rows for "save data"
  const head = ['bin_centre_MeV', 'counts', 'model', ...R.C.map(o => o.kind === 'level' ? 'lvl_' + o.E.toFixed(3) : o.kind)];
  const csv = [`# global fit  theta_cm ${G.cmLo}-${G.cmHi}  Ex ${G.lo}-${G.hi}  bin ${G.w.toFixed(4)} MeV`,
               `# shift ${R.shift.toFixed(4)}  sigma0 ${R.sig0.toFixed(4)}  dsigma ${G.dSig}  continuum ${R.psNote}`,
               `# chi2_lambda ${R.bc.toFixed(3)}  pearson ${R.pear.toFixed(3)}  ndf ${R.ndf}`,
               ...R.rows.map(rw => `# area ${rw.c.kind === 'level' ? rw.c.E.toFixed(3) : rw.c.kind} ${rw.counts.toFixed(2)} +- ${rw.err.toFixed(2)}`),
               head.join(',')];
  for (let b = 0; b < G.nb; ++b)
    csv.push([R.x[b].toFixed(5), R.n[b], R.m[b].toFixed(4), ...R.C.map(o => (o.a*o.f[b]).toFixed(4))].join(','));
  return {kind: 'gfit', csv, R, phMain, yP0, phP, ...c};
}

const gfFmt = (v, d) => isFinite(v) ? v.toFixed(d) : '&mdash;';
function gfReport(R) {
  const G = R.g;
  const sh = G.floatShift ? `${gfFmt(R.shift, 3)} &plusmn; ${gfFmt(R.shiftErr, 3)}` : `${gfFmt(R.shift, 3)} (fixed)`;
  const sg = G.floatSig ? `${gfFmt(R.sig0, 3)} &plusmn; ${gfFmt(R.sig0Err, 3)}` : `${gfFmt(R.sig0, 3)} (fixed)`;
  $('gfSum').innerHTML =
      `<span>N <b>${R.N}</b> in ${G.lo}&ndash;${G.hi} MeV</span>`
    + `<span>&chi;&sup2;<sub>&lambda;</sub>/ndf <b>${gfFmt(R.bc/R.ndf, 2)}</b></span>`
    + `<span>Pearson &chi;&sup2;/ndf <b>${gfFmt(R.pear/R.ndf, 2)}</b></span>`
    + `<span>ndf ${R.ndf} (${R.npar} par)</span>`
    + `<span>shift <b>${sh}</b> MeV</span><span>&sigma;<sub>0</sub> <b>${sg}</b> MeV</span>`
    + `<span>1n shape: ${R.psNote}</span>`
    + (R.covOk ? '' : `<span class="gf-flag">covariance singular &mdash; errors are &radic;N only</span>`)
    + `<span style="color:var(--ink-3)">${R.ms.toFixed(0)} ms</span>`;
  const rowsHtml = R.rows.map(rw => {
    const o = rw.c, isLv = o.kind === 'level';
    const sw = o.kind === 'ps' ? `background:var(--refdim)` : `background:var(${o.color})`;
    const flag = rw.railed ? '<span class="gf-flag">at zero (bound)</span>'
               : rw.floored ? 'error floored at &radic;N' : '';
    return `<tr><td><i class="sw" style="${sw}"></i>${o.kind === 'bg0' ? 'background (linear)' : o.label}</td>`
      + `<td>${isLv ? o.E.toFixed(3) : o.kind === 'n14' ? '&mdash;' : ''}</td>`
      + `<td>${isLv || o.kind === 'n14' ? o.mu.toFixed(3) : ''}</td>`
      + `<td>${isLv || o.kind === 'n14' ? o.sg.toFixed(3) : ''}</td>`
      + `<td>${rw.counts.toFixed(1)}</td><td>${rw.err.toFixed(1)}</td><td>${flag}</td></tr>`;
  }).join('');
  $('gfOut').innerHTML = `<table><thead><tr><th>component</th><th>E lit (MeV)</th><th>centroid</th><th>&sigma;</th>`
    + `<th>counts</th><th>&plusmn;</th><th></th></tr></thead><tbody>${rowsHtml}`
    + `<tr><td>model total / data</td><td></td><td></td><td></td><td>${R.modelSum.toFixed(1)}</td><td>${R.N}</td><td></td></tr>`
    + `</tbody></table>`;
  const lg = [`<span><i class="swatch" style="background:var(--ink);height:6px;width:6px;border-radius:3px"></i>data &plusmn;&radic;N</span>`,
              `<span><i class="swatch" style="background:var(--fit)"></i>total</span>`];
  for (const o of R.C) if (o.kind === 'level') lg.push(`<span><i class="swatch" style="background:var(${o.color})"></i>${o.label}</span>`);
  if (R.C.some(o => o.kind === 'n14')) lg.push(`<span><i class="swatch" style="background:var(--ink-2)"></i><sup>14</sup>N blend</span>`);
  if (R.C.some(o => o.kind === 'ps')) lg.push(`<span><i class="swatch" style="background:var(--refdim);height:8px"></i>1n phase space</span>`);
  if (R.C.some(o => o.kind === 'bg0')) lg.push(`<span><i class="swatch" style="background:var(--ink-3)"></i>background</span>`);
  lg.push(`<span><i class="swatch" style="background:var(--accent)"></i>S<sub>n</sub> = ${GF_SN}</span>`);
  $('legGf').innerHTML = lg.join('');
}

const gfCache = {key: null, R: null};
function globalFit(s) {
  const g = gfState();
  const fitKey = JSON.stringify([ACTIVE, s.ebeam, s.chi2, s.thLo, s.thHi, s.keLo, s.keHi, s.vzLo, s.vzHi,
                                 s.kcOn, s.kcDenom, s.kcPivot, Object.assign({}, g, {logY: false})]);
  if (gfCache.key !== fitKey) { gfCache.R = gfRun(s, g); gfCache.key = fitKey; }
  gfCache.R.g.logY = g.logY;
  views.cGf = gfDraw($('cGf'), gfCache.R);
  gfReport(gfCache.R);
}

hookTip('cGf', 'tipGf', (v, px, py) => {
  const R = v.R, G = R.g, b = Math.floor((px - PAD.l)/v.pw*G.nb);
  if (b < 0 || b >= G.nb) return null;
  const lines = [`E_x  ${(G.lo + b*G.w).toFixed(3)} .. ${(G.lo + (b + 1)*G.w).toFixed(3)} MeV`,
                 `counts ${R.n[b]}   model ${R.m[b].toFixed(1)}`,
                 `pull   ${((R.n[b] - R.m[b])/Math.sqrt(R.m[b])).toFixed(2)}`];
  const parts = R.C.map(o => [o, o.a*o.f[b]]).filter(([, y]) => y > 0.05).sort((p, q) => q[1] - p[1]);
  let bg = 0;
  for (const [o, y] of parts) {
    if (o.kind === 'bg0' || o.kind === 'bg1') { bg += y; continue; }
    lines.push(`  ${(o.kind === 'level' ? o.E.toFixed(3) : o.kind === 'n14' ? '14N blend' : '1n PS').padEnd(9)} ${y.toFixed(1)}`);
  }
  if (bg > 0.05) lines.push(`  ${'bg'.padEnd(9)} ${bg.toFixed(1)}`);
  return lines.join('\n');
});

/* ---- headless self-test: open the page with #gftest ---------------------------------------- */
function gfSelfTest() {
  const mN = /[#&]N=([0-9.]+)/.exec(location.hash);      // e.g. #gftest&N=6428.5714 for cat5_tc
  if (mN) { $('kcDenom').value = mN[1]; syncKcSliders(); }
  const s = readState(), out = {ebeam: s.ebeam, kcOn: s.kcOn, kcDenom: s.kcDenom};
  const base = gfSane(Object.assign(gfState(), {cmLo: 20, cmHi: 140, lo: 5.0, hi: 10.4, binKeV: 50, shift: -0.011,
    sig0: 0.132, dSig: 0.0123, muN: 9.178, sgN: 0.296, floatShift: false, floatSig: false, useN14: true, usePs: true,
    useBg: true, lv: [true, true, true, true, true, false, false, false]}));
  const summ = R => ({N: R.N, nll: +R.nll.toFixed(4), bc: +R.bc.toFixed(4), pear: +R.pear.toFixed(3), ndf: R.ndf,
    shift: +R.shift.toFixed(4), shiftErr: +(+R.shiftErr).toFixed(4), sig0: +R.sig0.toFixed(4), sig0Err: +(+R.sig0Err).toFixed(4),
    ms: Math.round(R.ms), nEval: R.nEval,
    rows: R.rows.map(rw => [rw.c.kind === 'level' ? rw.c.E : rw.c.kind, +rw.counts.toFixed(2), +rw.err.toFixed(2)])});
  const Rf = gfRun(s, Object.assign({}, base, {psSrc: 'file'}));
  out.bins = Array.from(Rf.n);
  out.fileFit = summ(Rf);
  // phase-space closure against the ROOT histogram: same generation window and smearing
  if (GF_HPS) {
    const fine = gfPsFine(159.75, 20, 140, true);
    const nbH = GF_HPS.v.length, bw = (GF_HPS.hi - GF_HPS.lo)/nbH;
    const js = Array.from({length: nbH}, (_, i) => gfPsSmeared(fine, GF_HPS.lo + (i + 0.5)*bw, () => 0.15));
    const norm = a => { const t = a.reduce((p, q) => p + q, 0); return a.map(v => v/t); };
    const A = norm(js), B = norm(GF_HPS.v.slice());
    const mean = a => a.reduce((p, v, i) => p + v*(GF_HPS.lo + (i + 0.5)*bw), 0);
    const pk = Math.max(...B);
    out.psClosure = {meanJS: +mean(A).toFixed(4), meanROOT: +mean(B).toFixed(4),
      maxAbsDiffOverPeak: +(Math.max(...A.map((v, i) => Math.abs(v - B[i])))/pk).toFixed(4),
      rmsJS: +Math.sqrt(A.reduce((p, v, i) => p + v*(GF_HPS.lo + (i + 0.5)*bw - mean(A))**2, 0)).toFixed(4),
      rmsROOT: +Math.sqrt(B.reduce((p, v, i) => p + v*(GF_HPS.lo + (i + 0.5)*bw - mean(B))**2, 0)).toFixed(4),
      firstAbove1pct: [A.findIndex(v => v > 0.01*Math.max(...A)), B.findIndex(v => v > 0.01*pk)].map(i => +(GF_HPS.lo + (i + 0.5)*bw).toFixed(3))};
    gfPsCache.key = null;
  }
  out.liveFixed = summ(gfRun(s, Object.assign({}, base, {psSrc: 'live'})));
  out.liveFloat = summ(gfRun(s, Object.assign({}, base, {psSrc: 'live', floatShift: true, floatSig: true})));
  const pre = document.createElement('pre'); pre.id = 'gfTest';
  pre.textContent = 'GFTEST ' + JSON.stringify(out) + ' GFEND';
  document.body.appendChild(pre);
}
if (/gftest/.test(location.hash)) setTimeout(gfSelfTest, 0);
if (/gfdark/.test(location.hash)) document.documentElement.setAttribute('data-theme', 'dark');
if (/gflight/.test(location.hash)) document.documentElement.setAttribute('data-theme', 'light');
