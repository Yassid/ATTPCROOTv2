# Figures of the trigger report, from the tables in trig/ (ana1_full, ana2_full, dsigma_fid, fid_checks) and report_data.py. Writes PNGs to trig/report/.
import numpy as np, json, re, os, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
D = '/mnt/f/a1975/sim_rebuild/trig/'; O = D + 'report/'; os.makedirs(O, exist_ok=True)
BL, OR, AQ, YE, INK, GR, LN = '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#0b0b0b', '#52514e', '#d9d8d4'
plt.rcParams.update({'font.size': 12, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.edgecolor': GR, 'xtick.color': GR, 'ytick.color': GR, 'axes.labelcolor': INK, 'axes.grid': True, 'grid.color': LN, 'grid.linewidth': .6, 'axes.axisbelow': True, 'legend.frameon': False, 'savefig.dpi': 130, 'savefig.facecolor': 'white', 'axes.titlesize': 12, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left'})
R = json.load(open(D + 'report_data.json')); T = json.load(open(D + 'report_traces.json'))
# 1 traces
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), sharey=True)
for a, k, t in [(ax[0], 'front', 'Front: run 17 evt 794, vz 262, r_res 31 mm'), (ax[1], 'far', 'Far: run 19 evt 24429, vz 753, r_res 11 mm')]:
    d = T[k]; a.plot(d['x'], d['ic'], color=GR, lw=1, label='IC'); a.plot(d['x'], d['mesh'], color=BL, lw=2, label='mesh'); a.axhline(320, color=OR, lw=1, ls='--'); a.set_title(t); a.set_xlabel('FRIB bin (80 ns)'); a.set_ylim(-100, 1400)
ax[0].set_ylabel('ADC, baseline subtracted'); ax[0].legend(loc='upper right'); ax[1].text(1210, 345, 'mesh edge, 320 ADC', color=INK, ha='left', fontsize=10)
plt.tight_layout(); plt.savefig(O + 'fig1_traces.png'); plt.close()
# 2 mesh height
h = T['hist']; e = np.array(h['e']); n = np.array(h['n'], float); fig, a = plt.subplots(figsize=(11, 3.4)); a.bar(e[:-1], np.where(n > 0, n, .5), width=18, align='edge', color=BL); a.set_yscale('log'); a.set_xlabel('mesh peak height [ADC, 5-bin mean]'); a.set_ylabel('triggered events / 20 ADC'); a.set_xlim(0, 1000); a.grid(axis='x', visible=False)
a.axvline(320, color=OR, lw=1, ls='--'); a.text(165, 1500, '280-300: 128\n300-320: 1798\n320-340: 15386', color=INK, fontsize=10); a.text(90, 9e4, 'no mesh pulse: 22 % of triggers', color=INK, fontsize=10)
plt.tight_layout(); plt.savefig(O + 'fig2_mesh_height.png'); plt.close()
# 3 early fire vs r_res, 4 timing vs vz
B = [b for b in R['bins']]; x = np.array([.5 * (b['a'] + b['b']) for b in B]); f = np.array([b['lt900']['f'] for b in B]); er = np.array([b['lt900']['e'] for b in B]); xe = np.array([.5 * (b['b'] - b['a']) for b in B])
fig, a = plt.subplots(figsize=(11, 3.9)); a.errorbar(x, f, er, xe, fmt='o', color=BL, ms=7, lw=1.5, capsize=0, label='vz < 900 mm'); b0 = B[0]['ge900']; a.errorbar([2.5], [b0['f']], [b0['e']], [2.5], fmt='s', mfc='white', color=OR, ms=7, lw=1.5, label='vz >= 900 mm, early by drift alone')
a.axvspan(0, 20, color='#eeeeea', zorder=0); a.text(10, .98, 'inside the beam hole, r < 20 mm', ha='center', va='top', color=GR, fontsize=10); a.axvline(28, color=INK, lw=1, ls='--'); a.text(28.4, .62, 'fiducial cut, 28 mm', color=INK, fontsize=10); a.set_xlim(0, 45); a.set_ylim(0, 1.03); a.set_xlabel('predicted residual radius at the pad plane, r_res [mm]'); a.set_ylabel('fraction with mesh within\n100 bins of the IC'); a.legend(loc='lower right')
plt.tight_layout(); plt.savefig(O + 'fig3_early_vs_rres.png'); plt.close()
s = open(D + 'ana1_full.txt').read(); rw = re.findall(r'vz\s+(\d+)-(\d+): n\s+(\d+) residual p16 p50 p84 \[\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\]', s); rw = np.array(rw, float); z = .5 * (rw[:, 0] + rw[:, 1])
fig, a = plt.subplots(figsize=(11, 3.7)); a.errorbar(z, rw[:, 4], [rw[:, 4] - rw[:, 3], rw[:, 5] - rw[:, 4]], fmt='o', color=BL, ms=7, lw=1.5, label='median, p16-p84'); zz = np.array([100, 1000]); a.plot(zz, -(1000 - zz) / .88, color=OR, ls='--', lw=1.5, label='mesh fires at the reaction time'); a.axhline(0, color=INK, lw=1); a.text(115, 25, 'mesh fires with the proton', color=INK, fontsize=10); a.set_xlabel('vertex position, vz [mm]'); a.set_ylabel('mesh crossing minus time\nfrom GET track start [bins]'); a.set_xlim(100, 1000); a.legend(loc='lower right')
plt.tight_layout(); plt.savefig(O + 'fig4_timing_vs_vz.png'); plt.close()
# 5 bound early fraction, vz x theta_cm
s2 = open(D + 'ana2_full.txt').read(); blk = s2.split('bound Ex')[1].split('Ex 3-5')[0]; rows = re.findall(r'^\s+(\d+)-(\d+)\s+(.*)$', blk, re.M); N = []; F = []; Y = []
for a0, a1, r in rows: v = np.array(r.split(), float).reshape(-1, 3); N.append(v[:, 0]); F.append(v[:, 1]); Y.append(a0 + '-' + a1)
N, F = np.array(N), np.array(F); C = ['10-20', '20-25', '25-30', '30-35', '35-40']; fig, a = plt.subplots(figsize=(7.4, 4.3)); M = np.ma.masked_where(N < 5, F); cm = plt.get_cmap('Blues').copy(); cm.set_bad('#eeeeea'); a.imshow(M, cmap=cm, vmin=0, vmax=1.15, aspect='auto', origin='lower'); a.grid(False); a.spines[:].set_visible(False); a.tick_params(length=0)
for i in range(N.shape[0]):
    for j in range(N.shape[1]):
        a.text(j, i, f'{F[i,j]:.2f}\nn={int(N[i,j])}' if N[i, j] >= 5 else f'n={int(N[i,j])}', ha='center', va='center', fontsize=10, color='white' if (N[i, j] >= 5 and F[i, j] > .7) else INK)
a.set_xticks(range(5), C); a.set_yticks(range(len(Y)), Y); a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel('vz [mm]'); plt.tight_layout(); plt.savefig(O + 'fig5_bound_cells.png'); plt.close()
# 6 dsigma vs cut, 7 angular distribution
tx = open(D + 'dsigma_fid.txt').read(); P = {}
for l in tx.splitlines():
    m = re.match(r'r_res>=\s*(\d+) dsdo(.*?)\|', l)
    if m: P[int(m[1])] = [None if c == '-' else tuple(map(float, re.match(r'([\d.]+)\(([\d.]+)\)', c).groups())) for c in m[2].split()]
Yv = {int(m[0]): list(map(int, m[1].split())) for m in re.findall(r'r_res>=\s*(\d+) dsdo.*\n\s+Y\s+([\d ]+?)\s+eps', tx)}
fig, a = plt.subplots(figsize=(11, 4)); cuts = [0, 16, 20, 22, 24, 26, 28, 30]
for b, lab, c, off in [(2, '20-25', YE, -.45), (3, '25-30', OR, -.15), (4, '30-35', AQ, .15), (5, '35-40', BL, .45)]:
    p = [(r, *P[r][b]) for r in cuts if P[r][b] and Yv[r][b] >= 10]; p = np.array(p); a.errorbar(p[:, 0] + off, p[:, 1], p[:, 2], fmt='o-', color=c, ms=6, lw=1.5, capsize=0, label=lab + r'$^\circ$')
a.set_xlabel('cut: r_res >= [mm]'); a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); a.set_ylim(0, 12.5); a.legend(title=r'$\theta_{cm}$', loc='upper left', ncol=4); plt.tight_layout(); plt.savefig(O + 'fig6_dsigma_vs_cut.png'); plt.close()
TB = np.arange(10, 50, 5.); mid = .5 * (TB[:-1] + TB[1:]); L2 = np.array([120.9, 91.5, 62.6, 38.5, 22., 12.8, 8.6]); S = .251
fig, a = plt.subplots(figsize=(11, 4.6)); a.step(TB, np.r_[S * L2, S * L2[-1]], where='post', color=INK, lw=1.5, label=f'l=2 sum x {S}, fitted to 25-40$^\\circ$ only')
A0 = np.array([3.68, 5.85, 6.82, 3.49]); E0 = np.array([.62, .64, .76, .49]); a.errorbar(mid[2:6] - .5, A0, E0, fmt='s', mfc='white', color=GR, ms=7, lw=1.5, label='adopted 2026-09-24, vz 200-400')
for r, c, mk, d in [(26, BL, 'o', 0), (28, OR, 'D', .5)]:
    idx = [i for i in range(3, 6)]; a.errorbar(mid[idx] + d, [P[r][i][0] for i in idx], [P[r][i][1] for i in idx], fmt=mk, color=c, ms=7, lw=1.5, label=f'r_res >= {r} mm')
a.errorbar([mid[6] + .5], [P[28][6][0]], [P[28][6][1]], fmt='D', mfc='white', color=OR, ms=7, lw=1.5, label='40-45$^\\circ$ at r_res >= 28, not in the fit'); a.set_yscale('log'); a.set_xlim(10, 45); a.set_ylim(.2, 60); a.set_xlabel(r'$\theta_{cm}$ [deg]'); a.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); a.legend(loc='upper right', fontsize=10); plt.tight_layout(); plt.savefig(O + 'fig7_angdist.png'); plt.close()
print({r: Yv[r] for r in Yv})
