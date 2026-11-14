import numpy as np, uproot, io, contextlib, runpy, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()): G=runpy.run_path('/mnt/f/a1975/sim_rebuild/trig/dsigma_fid.py')
run,L2,TH,TB=G['run'],G['L2'],G['TH'],G['TB']; mid=.5*(TB[:-1]+TB[1:])
X=uproot.open('/mnt/f/a1975/sim_rebuild/r280/p_gs_census.root')['c'].arrays(['thl','vz','ke','ex','tcm','sv'],library='np'); M1=16.0147013*931.49401
def rres(ke,th,vz):
    t=np.radians(th); p=np.sqrt(ke**2+2*ke*938.272); E=183.871-0.011071*vz; pb=np.sqrt(E**2+2*E*M1); return (1000-vz)*np.tan(np.arctan2(p*np.sin(t),pb-p*np.cos(t)))
def fit1(x,er,m,u): S=np.sum((x*m/er**2)[u])/np.sum((m*m/er**2)[u]); return S,np.sqrt(1/np.sum((m*m/er**2)[u])),np.sum((((x-S*m)/er)**2)[u])
R=28; x,er,Y,B,e,Le=run(R)
b=(X['ex']>-1.5)&(X['ex']<1.0)&(X['vz']>200)&(X['vz']<900)&(X['thl']>=89)&(X['thl']<=161)&(X['ke']>=1.2)&(rres(X['ke'],X['thl'],X['vz'])>=R)
sv=np.array([np.mean(X['sv'][b&(X['tcm']>=a)&(X['tcm']<c)]>0) if (b&(X['tcm']>=a)&(X['tcm']<c)).sum()>20 else 1.0 for a,c in zip(TB[:-1],TB[1:])])
xc,ec=x/sv,er/sv; u=(Y>=3)&(TB[:-1]>=25)&(TB[1:]<=45); S=fit1(xc,ec,L2,u)
PLBc=np.array([0.03,0.64,0.62]); plb=TH@PLBc   # states: gs 3/2+, 1/2+, 5/2+ as in dsigma_fid2
fig,ax=plt.subplots(figsize=(8,6))
ax.errorbar(mid[u]-.3,x[u],er[u],fmt='o',mfc='white',color='0.5',capsize=2,label='R>=28, previous (40-45 without Spyral term)')
ax.errorbar(mid[u]+.3,xc[u],ec[u],fmt='D',color='k',capsize=2,label='R>=28, Spyral-validity term in acceptance')
ax.errorbar(mid[u]+.6,xc[u]/0.84,ec[u]/0.84,fmt='s',color='tab:purple',capsize=2,label='same / 0.84 (environment loss, re-weighted)')
tt=np.r_[TB[:-1],TB[-1]]
ax.step(tt,np.r_[S[0]*L2,S[0]*L2[-1]],where='post',color='tab:red',label=f'l=2 sum x {S[0]:.3f}({S[1]:.3f}), chi2 {S[2]:.1f}/{u.sum()-1}')
ax.step(tt,np.r_[plb,plb[-1]],where='post',color='tab:blue',ls='--',label='PLB 811 C2S (.03,.64,.62) on our FRESCO curves')
ax.set_yscale('log'); ax.set_xlim(20,47); ax.set_ylim(.2,40); ax.grid(alpha=.3,which='both'); ax.legend(fontsize=8)
ax.set_xlabel(r'$\theta_{cm}$ [deg]'); ax.set_ylabel(r'd$\sigma$/d$\Omega$ [mb/sr]'); ax.set_title(r'$^{16}$C(d,p)$^{17}$C bound, trigger-safe fiducial r$_{res}\geq$28 mm (NOT adopted)',fontsize=10)
plt.tight_layout(); plt.savefig('/mnt/f/a1975/sim_rebuild/trig/fiducial_result_spyterm.png',dpi=100)
print('data/PLB:', np.round((xc/plb)[u],2), ' with /0.84:', np.round((xc/plb/0.84)[u],2))
