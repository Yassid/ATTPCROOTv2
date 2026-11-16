# Per-bin BW fits from unb_bw.py (reads unb_bw_keep.npy): data, total model, each state (x acceptance), phase space, pulls.
import os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw.py').read().split("\nif __name__ == '__main__':")[0])
tag = os.environ.get('TAG', 'bw'); keep = np.load(U + f'unb_{tag}_keep.npy', allow_pickle=True).item()
cols = ['tab:green', 'tab:olive', 'tab:orange', 'tab:purple', 'tab:cyan', 'tab:brown']
vzs = sorted(keep); fig = plt.figure(figsize=(20, 5.6 * len(vzs))); gs = fig.add_gridspec(2 * len(vzs), NB, height_ratios=[3, 1] * len(vzs), hspace=.1, wspace=.18)
for r, vz0 in enumerate(vzs):
   k = keep[vz0]; H, A, x, ng = k['H'], k['A'], k['x'], k['ng']; S = shapes(x[0], x[1], x[3] if ng == 4 else 1.0)
   for i in range(NB):
      o = ng + (NS + 1) * i; ps = x[o] * np.sqrt(np.maximum(cen - SN, 0)) * (1 + x[2] * (cen - 3)); Y = x[o + 1:o + 1 + NS]
      tot = A[i] * (ps + Y @ S)
      ax = fig.add_subplot(gs[2 * r, i]); rx = fig.add_subplot(gs[2 * r + 1, i], sharex=ax)
      ax.errorbar(cen, H[i], np.sqrt(np.maximum(H[i], 1)), fmt='o', ms=3, color='k', lw=1)
      ax.plot(cen, tot, 'r-', lw=1.6, label=f'total, -2lnL {pnll(tot, H[i]):.1f}/{len(cen)}')
      ax.plot(cen, A[i] * ps, 'b--', lw=1, label='1n phase space x tilt')
      for s in range(NS):
         if Y[s] > 0.5: ax.plot(cen, A[i] * Y[s] * S[s], color=cols[s], lw=1.2, label=f'{ST[s][3]}: Y {Y[s]:.0f}')
      ax2 = ax.twinx(); ax2.plot(cen, A[i], color='gray', ls=':', lw=1); ax2.set_ylim(0, 1.6); ax2.set_yticks([0, .5, 1]); ax2.tick_params(labelsize=7, colors='gray')
      ax.set_title(f'vz $\\geq$ {vz0} mm, $\\theta_{{cm}}$ {TB[i]:.0f}-{TB[i+1]:.0f} (N {int(H[i].sum())}); dotted = acceptance', fontsize=9)
      ax.legend(fontsize=6.5, loc='upper left'); ax.grid(alpha=.3); ax.set_ylim(0, max(H[i].max() + 3 * np.sqrt(H[i].max()), 1)); plt.setp(ax.get_xticklabels(), visible=False)
      pull = (H[i] - tot) / np.sqrt(np.maximum(tot, 1)); rx.bar(cen, pull, width=.09, color='gray'); rx.axhline(0, color='k', lw=.6); rx.set_ylim(-3.5, 3.5); rx.grid(alpha=.3)
      rx.set_ylabel('pull', fontsize=8); rx.set_xlabel(r'E$_x$ [MeV] (calibrated frame)' if r == len(vzs) - 1 else '')
      if i == 0: ax.set_ylabel('counts / 100 keV')
   SUM = globals().get('SUM', []) + [f'vz >= {vz0}: offset {x[0]:+.3f} MeV, sigma_res {x[1]:.3f}, tilt {x[2]:+.3f}']
fig.suptitle('BW E0 = table + offset, Gamma fixed (table), penetrabilities.   ' + '   |   '.join(SUM), fontsize=10, y=0.995)
fig.savefig(U + f'unb_{tag}_fits_2026-11-17.png', dpi=100, bbox_inches='tight'); print('saved')
