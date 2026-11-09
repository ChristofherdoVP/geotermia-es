from config import *
import numpy as np, harmonica as hm, warnings; warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
plt.rcParams.update({'font.size':9,'axes.titlesize':10})
from bouligand import bouligand, spectrum, t, H
P=DATA
g=hm.load_oasis_montaj_grid(P+'M_SemIGRF.grd').coarsen(easting=10,northing=10,boundary='trim').mean()
x=g.easting.values; y=g.northing.values; dxa=500.
def spec_aero(cx,cy,L=150e3):
    ix=(x>=cx-L/2)&(x<cx+L/2); iy=(y>=cy-L/2)&(y<cy+L/2); w=g.values[np.ix_(iy,ix)]
    w=np.where(np.isnan(w),np.nanmean(w),w); w=w-w.mean()
    yy,xx=np.mgrid[:w.shape[0],:w.shape[1]]; A=np.c_[xx.ravel(),yy.ravel(),np.ones(xx.size)]; w=w-(A@np.linalg.lstsq(A,w.ravel(),rcond=None)[0]).reshape(w.shape)
    Pw=np.abs(np.fft.fft2(w*np.outer(np.hanning(w.shape[0]),np.hanning(w.shape[1]))))**2
    k1=np.fft.fftfreq(w.shape[1],dxa/1e3)*2*np.pi; k0=np.fft.fftfreq(w.shape[0],dxa/1e3)*2*np.pi
    K=np.hypot(*np.meshgrid(k1,k0)); dk=2*np.pi/(L/1e3); b=np.arange(dk/2,4,dk); i=np.digitize(K.ravel(),b)
    kr=np.array([K.ravel()[i==j].mean() for j in range(1,len(b))]); lp=np.array([np.log(Pw.ravel()[i==j].mean()) for j in range(1,len(b))])
    return kr,lp
fig,ax=plt.subplots(1,2,figsize=(10,4))
res=[]
for (cx,cy),c in zip([(339e3,7909e3),(354e3,7889e3),(330e3,7840e3)],['C0','C1','C2']):
    kr,lp=spec_aero(cx,cy); m=(kr>=0.3)&(kr<=3.0)
    r=least_squares(lambda p: bouligand(kr[m],p[0],p[1],60.0,p[2])-lp[m],[3,1.0,lp[m][0]],bounds=([0.5,0.01,-np.inf],[6,10,np.inf]))
    J=r.jac; cov=np.linalg.pinv(J.T@J)*(r.fun**2).sum()/(m.sum()-3); e=np.sqrt(np.diag(cov)); res.append((cx,cy,r.x[0],e[0],r.x[1],e[1]))
    ax[0].plot(kr,lp-lp[m][0],'.',color=c,ms=4,label=f'window at ({cx/1e3:.0f}, {cy/1e3:.0f}) km: β = {r.x[0]:.2f} ± {e[0]:.2f}')
    ax[0].plot(kr[m],bouligand(kr[m],*r.x[:2],60.0,r.x[2])-lp[m][0],'-',color=c,lw=1)
ax[0].axvspan(0.3,3.0,color='0.93',zorder=0); ax[0].set_xlabel('k (rad km$^{-1}$)'); ax[0].set_ylabel('ln P (normalized)')
ax[0].set_title('(a) Aeromagnetic data, 150-km windows'); ax[0].legend(fontsize=7,frameon=False); ax[0].set_xlim(0,3.5)
xc,yc=t.transform(-40.63,-19.54); kr,lp,sd=spectrum(xc,yc,300e3); m=kr<=0.4
ax[1].plot(kr[m],lp[m]-lp[m][0],'ko',ms=4,label='EMAG2v3, 300-km window on Colatina')
beta=2.75
for dz,c in [(15,'C3'),(30,'C1'),(60,'C2'),(500,'C0')]:
    r=least_squares(lambda p: bouligand(kr[m],beta,p[0],dz,p[1])-lp[m],[6,lp[m][0]],bounds=([0.1,-np.inf],[30,np.inf]))
    rms=np.sqrt((r.fun**2).mean())
    lab=(f'base {r.x[0]+dz-H:.0f} km bsl' if dz<500 else 'no base (Δz → ∞)')+f'; rms = {rms:.2f}'
    ax[1].plot(kr[m],bouligand(kr[m],beta,r.x[0],dz,r.x[1])-lp[m][0],'-',color=c,lw=1.2,label=lab)
    print(dz,'rms',round(rms,3),'zt',round(r.x[0]-H,2))
ax[1].set_xlabel('k (rad km$^{-1}$)'); ax[1].set_title('(b) Long wavelengths (β = 2.75; top and level fitted)'); ax[1].legend(fontsize=7,frameon=False)
for a in ax: a.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('../figures/Fig4_spectra.png',dpi=300)
for r_ in res: print('aeromag window (%.0f,%.0f): beta %.2f±%.2f, zt %.2f±%.2f km below flight'%(r_[0]/1e3,r_[1]/1e3,*r_[2:]))
