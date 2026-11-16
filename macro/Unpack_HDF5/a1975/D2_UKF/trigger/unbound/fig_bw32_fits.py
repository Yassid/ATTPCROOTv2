# Per-bin BW fits at offset -0.32 (unb_bw_force32.py keep): all-free vs 2.980 tied at the YSOX ratio, one figure per vz window.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw_force.py').read().split("\nif __name__ == '__main__':")[0])
keep = np.load(U + 'unb_bw_force32_keep.npy', allow_pickle=True).item()
cols = ['tab:green', 'tab:olive', 'tab:orange', 'tab:purple', 'tab:cyan', 'tab:brown']
VAR = ['all 6 free, offset fixed -0.32', '2.980 tied at YSOX ratio, offset -0.32']
for vz0 in (100, 500):
   fig = plt.figure(figsize=(20, 11.5)); gs = fig.add_gridspec(4, NB, height_ratios=[3, 1, 3, 1], hspace=.12, wspace=.18)
   for r, lab in enumerate(VAR):
      L, g, Y, H, A, qs, mode, fo, c = keep[(vz0, lab)]; S = shapes(g[0], g[1], 1.0)
      for i in range(NB):
         q = qs[i]; ps = q[0] * np.sqrt(np.maximum(cen - SN, 0)) * (1 + g[2] * (cen - 3))
         comps = [(ST[s][3], Y[i][s] * S[s]) for s in range(NS)]          # Y already split for the tied pair
         tot = A[i] * (ps + sum(cc for _, cc in comps))
         ax = fig.add_subplot(gs[2 * r, i]); rx = fig.add_subplot(gs[2 * r + 1, i], sharex=ax)
         ax.errorbar(cen, H[i], np.sqrt(np.maximum(H[i], 1)), fmt='o', ms=3, color='k', lw=1)
         ax.plot(cen, tot, 'r-', lw=1.6, label=f'total, -2lnL {pnll(tot, H[i]):.1f}/{len(cen)}')
         ax.plot(cen, A[i] * ps, 'b--', lw=1, label='1n phase space x tilt')
         for s, (nm, cc) in enumerate(comps):
            if Y[i][s] > 0.5: ax.plot(cen, A[i] * cc, color=cols[s], lw=1.3, label=f'{nm}: Y {Y[i][s]:.0f}')
         ax2 = ax.twinx(); ax2.plot(cen, A[i], color='gray', ls=':', lw=1); ax2.set_ylim(0, 1.6); ax2.set_yticks([0, .5, 1]); ax2.tick_params(labelsize=7, colors='gray')
         ax.set_title(('ALL FREE' if r == 0 else '2.980 TIED (YSOX)') + f', $\\theta_{{cm}}$ {TB[i]:.0f}-{TB[i+1]:.0f} (N {int(H[i].sum())})', fontsize=9)
         ax.legend(fontsize=6.5, loc='upper left'); ax.grid(alpha=.3); ax.set_ylim(0, max(H[i].max() + 3 * np.sqrt(H[i].max()), 1)); plt.setp(ax.get_xticklabels(), visible=False)
         pull = (H[i] - tot) / np.sqrt(np.maximum(tot, 1)); rx.bar(cen, pull, width=.09, color='gray'); rx.axhline(0, color='k', lw=.6); rx.set_ylim(-3.5, 3.5); rx.grid(alpha=.3)
         rx.set_ylabel('pull', fontsize=8); rx.set_xlabel(r'E$_x$ [MeV] (calibrated frame; bound blend at 0)' if r == 1 else '')
         if i == 0: ax.set_ylabel('counts / 100 keV')
   La, Lb = keep[(vz0, VAR[0])][0], keep[(vz0, VAR[1])][0]; gA, gB = keep[(vz0, VAR[0])][1], keep[(vz0, VAR[1])][1]
   fig.suptitle(f'vz >= {vz0} mm, BW offset FIXED -0.32 MeV (bound-blend centroid +0.32 from YSOX/WBP/PLB C2S). Top: all 6 states free, -2lnL {La:.1f}, '
                f'sigma {gA[1]:.3f}, tilt {gA[2]:+.2f}. Bottom: Y(2.980) = YSOX ratio x Y(2.763), -2lnL {Lb:.1f} (+{Lb - La:.1f}, 4 fewer pars), sigma {gB[1]:.3f}. '
                'Dotted = acceptance.', fontsize=10, y=0.995)
   fig.savefig(U + f'unb_bw32_fits_vz{vz0}_2026-11-17.png', dpi=100, bbox_inches='tight'); print('saved', vz0)
