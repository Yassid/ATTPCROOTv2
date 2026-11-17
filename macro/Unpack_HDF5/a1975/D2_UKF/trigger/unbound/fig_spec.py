# Inclusive Ex spectrum of the working unbound selection (theta_lab 89-THI, KE >= KEMIN, chi2, vz windows), with the BW fit (offset -0.32, all free)
# summed over the theta_cm bins: total, each state and the phase space, all x acceptance. Data shown from -1.5 (bound peak) to 6 MeV.
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw_force.py').read().split("\nif __name__ == '__main__':")[0])
keep = np.load(U + f'unb_bw_force32{SUF}_keep.npy', allow_pickle=True).item(); LAB = 'all 6 free, offset fixed -0.32'
cols = ['tab:green', 'tab:olive', 'tab:orange', 'tab:purple', 'tab:cyan', 'tab:brown']
fig, ax = plt.subplots(1, 2, figsize=(17, 6))
for j, vz0 in enumerate((100, 500)):
   L, g, Y, H, A, qs, mode, fo, c = keep[(vz0, LAB)]; S = shapes(g[0], g[1], 1.0); a = ax[j]
   ww = dsel & (vz >= vz0) & (vz < 900) & (dc >= TB[0]) & (dc < TB[-1])
   e2 = np.arange(-1.5, 6.0001, .1); h2 = np.histogram(dx[ww], e2)[0]
   a.errorbar(.5 * (e2[1:] + e2[:-1]), h2, np.sqrt(np.maximum(h2, 1)), fmt='o', ms=3, color='k', lw=1, label=f'data, theta_cm {TB[0]:.0f}-{TB[-1]:.0f} (N {ww.sum()})')
   tot = np.zeros_like(cen); comp = np.zeros((NS, len(cen))); psum = np.zeros_like(cen)
   for i in range(NB):
      ps = qs[i][0] * np.sqrt(np.maximum(cen - SN, 0)) * (1 + g[2] * (cen - 3)); psum += A[i] * ps
      for s in range(NS): comp[s] += A[i] * Y[i][s] * S[s]
   tot = psum + comp.sum(0)
   a.plot(cen, tot, 'r-', lw=2, label=f'BW fit (sum of {NB} bins), -2lnL {L:.0f}/{H.size - (3 - 1) - NB * (NS + 1)}')
   a.plot(cen, psum, 'b--', lw=1.2, label='1n phase space x tilt')
   for s in range(NS): a.plot(cen, comp[s], color=cols[s], lw=1.4, label=f'{ST[s][3]}: {comp[s].sum():.0f} counts in window')
   for s in range(NS): a.axvline(ST[s][0] + g[0], color=cols[s], ls=':', lw=.8)
   a.axvspan(-1.5, EXLO, color='gray', alpha=.08); a.text(-1.3, h2.max() * .95, 'bound peak\n(not fitted here)', fontsize=8, color='gray')
   a.set_xlabel(r'E$_x$ [MeV] (calibrated: bound blend at 0, so state energies sit 0.32 MeV low)'); a.set_ylabel('counts / 100 keV')
   a.set_title(f'vz $\\geq$ {vz0} mm, $\\theta_{{lab}}$ 89-{THI:.0f}, KE $\\geq$ {KEMIN} MeV; offset {g[0]:+.2f}, $\\sigma$ {g[1]:.3f}', fontsize=10)
   a.set_xlim(-1.5, 6); a.set_ylim(0, None); a.grid(alpha=.3); a.legend(fontsize=7.5, loc='upper right')
fig.suptitle('16C(d,p)17C, Ex spectrum of the unbound selection with the Breit-Wigner fit (C17_dp_fits states, Gamma fixed); NOT ADOPTED', fontsize=11)
fig.tight_layout(); fig.savefig(U + f'unb_spec{SUF}_2026-11-17.png', dpi=110); print('saved')
