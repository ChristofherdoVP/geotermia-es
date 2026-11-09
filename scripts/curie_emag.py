import numpy as np, xarray as xr, warnings; warnings.filterwarnings('ignore')
from pyproj import Transformer; from scipy.interpolate import griddata
da=xr.open_dataarray('emag2v3_recorte.nc')
LON,LAT=np.meshgrid(da.lon,da.lat); t=Transformer.from_crs(4326,32724,always_xy=True)
X,Y=t.transform(LON,LAT); m=~np.isnan(da.values)
xc,yc=t.transform(-40.63,-19.54); dx=2000.
def window(L):
  xs=np.arange(xc-L/2,xc+L/2,dx); ys=np.arange(yc-L/2,yc+L/2,dx); GX,GY=np.meshgrid(xs,ys)
  sel=m&(abs(X-xc)<L/2+3e4)&(abs(Y-yc)<L/2+3e4)
  return griddata((X[sel],Y[sel]),da.values[sel],(GX,GY),method='cubic')
def spec(w,L):
  w=np.where(np.isnan(w),np.nanmean(w),w); w=w-w.mean()
  Pw=np.abs(np.fft.fft2(w*np.outer(np.hanning(w.shape[0]),np.hanning(w.shape[1]))))**2
  k1=np.fft.fftfreq(w.shape[1],dx/1e3)*2*np.pi; k0=np.fft.fftfreq(w.shape[0],dx/1e3)*2*np.pi
  K=np.hypot(*np.meshgrid(k1,k0)); dk=2*np.pi/(L/1e3); b=np.arange(dk/2,K.max(),dk); i=np.digitize(K.ravel(),b)
  kr=np.array([K.ravel()[i==j].mean() for j in range(1,len(b))]); pr=np.array([Pw.ravel()[i==j].mean() for j in range(1,len(b))])
  return kr,pr
out={}
for L in [200e3,300e3,400e3]:
  w=window(L); kr,pr=spec(w,L); out[L]=(kr,pr)
  print(f'--- janela {int(L/1e3)} km  NaN {100*np.isnan(w).mean():.1f}%  Δk={kr[1]-kr[0]:.3f} rad/km')
  zbs=[]
  for zmax in [0.10,0.15,0.20,0.25]:
    b=kr<=zmax; a=(kr>=0.25)&(kr<=0.6)
    if b.sum()<3: continue
    z0=-np.polyfit(kr[b],np.log(np.sqrt(pr[b])/kr[b]),1)[0]; zt=-np.polyfit(kr[a],np.log(np.sqrt(pr[a])),1)[0]-4  # desconta 4 km de altitude
    zb=2*z0-zt; zbs.append(zb)
    print(f'  Z0 ajuste k<={zmax:.2f} ({b.sum()} pts): Z0={z0:5.1f}  Zt={zt:4.1f}  Zb={zb:5.1f} km (abaixo da superfície; Z0 inclui 4 km de altitude: Zb_corr={zb-4:5.1f})')
  print('  faixa Zb_corr: %.1f a %.1f km'%(min(zbs)-4,max(zbs)-4))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(11,4.5))
for L,(kr,pr) in out.items():
  ax[0].plot(kr,np.log(np.sqrt(pr)),'.-',label=f'{int(L/1e3)} km'); ax[1].plot(kr,np.log(np.sqrt(pr)/kr),'.-',label=f'{int(L/1e3)} km')
ax[0].set_xlim(0,1); ax[1].set_xlim(0,0.4); ax[0].set_xlabel('k (rad/km)'); ax[1].set_xlabel('k (rad/km)')
ax[0].set_ylabel('ln √P'); ax[1].set_ylabel('ln(√P/k)'); ax[0].set_title('Topo das fontes'); ax[1].set_title('Centroide')
for a in ax: a.legend(frameon=False); a.spines[['top','right']].set_visible(False)
plt.suptitle('EMAG2v3, janelas centradas em Colatina'); plt.tight_layout(); plt.savefig('curie_emag2_espectros.png',dpi=140)
