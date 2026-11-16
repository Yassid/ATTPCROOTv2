# Same check with the offset FIXED at -0.32 MeV = minus the bound-blend centroid (+0.317..0.321 for YSOX, WBP and PLB 811 C2S x ADWA over theta_cm
# 10-45: the 5/2+ at 0.331 carries 93-96 % of the bound peak). unb_bw_force.py run() reused.
exec(open('/mnt/f/a1975/sim_rebuild/unb_hdbtau/unb_bw_force.py').read().split("\nif __name__ == '__main__':")[0])
V = [('all 6 free, offset fixed -0.32', 'free', -0.32, None), ('2.763 dropped, offset fixed -0.32', 'drop', -0.32, None),
     ('2.980 tied at YSOX ratio, offset -0.32', 'tie', -0.32, .2043 / .0684), ('2.980 tied at WBP ratio, offset -0.32', 'tie', -0.32, .2029 / .0372)]
keep = {}
for vz0 in (100, 500):
   print(f'\n=== vz >= {vz0}   (free-offset baseline: -2lnL 102.6 at vz>=100, 99.9 at vz>=500)      summed yields: 2.763  2.980  3.661  4.231  4.841')
   for lab, mode, fo, c in V:
      L, ndf, g, Y, H, A, qs = run(vz0, mode, fo, c); Ys = Y.sum(0)
      print(f'  {lab:42s} -2lnL {L:6.1f}/{ndf}  sigma {g[1]:.3f} tilt {g[2]:+.3f} | ' + '  '.join(f'{v:6.0f}' for v in Ys[:5]))
      keep[(vz0, lab)] = (L, g, Y, H, A, qs, mode, fo, c)
np.save(U + f'unb_bw_force32{SUF}_keep.npy', keep, allow_pickle=True)
