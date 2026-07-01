#!/usr/bin/env python
"""Compare 17C Ex spectra: default vs loose vs adaptive TC clustering.
Adaptive = default per event, except events flagged as recovered multi-turn spirals
(loose_max >= RATIO*def_max AND loose_turns >= TURNS) use the loose result.
"""
import numpy as np, uproot, polars as pl
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

RATIO, TURNS = 1.3, 1.5
D = "spyral_compare/"

def load(cache):
    a = uproot.open(D + cache)['pk'].arrays(library='np')
    return pl.DataFrame({k: a[k] for k in a})

def_df = load("cache_def.root")
loose_df = load("cache_loose.root")
flag = pl.read_csv(D + "spiral_flag.csv")
flagged = set(flag.filter((pl.col('loose_max') >= RATIO*pl.col('def_max')) &
                          (pl.col('loose_turns') >= TURNS) & (pl.col('def_max') >= 20))['event'].to_list())
print(f"flagged spiral events: {len(flagged)}")

# adaptive: default rows for non-flagged events + loose rows for flagged events
adapt = pl.concat([def_df.filter(~pl.col('event').is_in(list(flagged))),
                   loose_df.filter(pl.col('event').is_in(list(flagged)))])

def spec(df, name):
    ex = df['ex'].to_numpy()
    n = len(ex)
    gs = ((ex > -1.5) & (ex < 1.5)).sum()      # near ground state
    return ex, n, gs

bins = np.linspace(-5, 15, 100)
fig, ax = plt.subplots(1, 2, figsize=(14, 5))
rows = []
for df, name, col in [(def_df, 'default', 'tab:blue'), (loose_df, 'loose', 'tab:orange'),
                      (adapt, 'adaptive', 'tab:green')]:
    ex, n, gs = spec(df, name)
    ax[0].hist(ex, bins=bins, histtype='step', lw=2, color=col, label=f'{name} (n={n}, gs={gs})')
    rows.append((name, n, gs))
ax[0].set_xlabel('Ex(17C) [MeV]'); ax[0].set_ylabel('proton candidates'); ax[0].legend()
ax[0].set_title('17C Ex spectrum: default vs loose vs adaptive')
ax[0].axvline(0, color='k', ls=':')

# difference adaptive - default
exd = np.histogram(def_df['ex'].to_numpy(), bins=bins)[0]
exa = np.histogram(adapt['ex'].to_numpy(), bins=bins)[0]
ctr = 0.5*(bins[1:]+bins[:-1])
ax[1].bar(ctr, exa-exd, width=(bins[1]-bins[0]), color='tab:green', alpha=0.7)
ax[1].set_xlabel('Ex(17C) [MeV]'); ax[1].set_ylabel('adaptive - default')
ax[1].set_title('Net change from adaptive recovery'); ax[1].axhline(0, color='k')
plt.tight_layout(); plt.savefig(D + "plots/ex_compare.png", dpi=110)

print("\n%-10s %8s %8s" % ("variant", "Ncand", "near-gs"))
for nm, n, gs in rows:
    print("%-10s %8d %8d" % (nm, n, gs))
print(f"\nadaptive adds {rows[2][1]-rows[0][1]:+d} candidates vs default ({rows[2][2]-rows[0][2]:+d} near g.s.)")
print(f"loose    adds {rows[1][1]-rows[0][1]:+d} candidates vs default ({rows[1][2]-rows[0][2]:+d} near g.s.)")
print("saved", D + "plots/ex_compare.png")
