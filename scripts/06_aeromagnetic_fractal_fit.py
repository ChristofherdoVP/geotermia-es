from config import *
from bouligand import bouligand, least_squares, np, H, spectrum, t
import harmonica as hm, warnings; warnings.filterwarnings('ignore')
# 1) β e topo a partir do aeromag de alta resolução (grade 500 m), janela 150 km, k alto
g=hm.load_oasis_montaj_grid(DATA+'M_SemIGRF.grd').coarsen(easting=10,northing=10,boundary='trim').mean()
x=g.easting.values; y=g.northing.values; dxa=500.
def spec_aero(cx,cy,L=150e3):
    ix=(x>=cx-L/2)&(x<cx+L/2); iy=(y>=cy-L/2)&(y<cy+L/2); w=g.values[np.ix_(iy,ix)]
    w=np.where(np.isnan(w),np.nanmean(w),w); w=w-w.mean()
    yy,xx=np.mgrid[:w.shape[0],:w.shape[1]]; A=np.c_[xx.ravel(),yy.ravel(),np.ones(xx.size)]; w=w-(A@np.linalg.lstsq(A,w.ravel(),rcond=None)[0]).reshape(w.shape)
    P=np.abs(np.fft.fft2(w*np.outer(np.hanning(w.shape[0]),np.hanning(w.shape[1]))))**2
    k1=np.fft.fftfreq(w.shape[1],dxa/1e3)*2*np.pi; k0=np.fft.fftfreq(w.shape[0],dxa/1e3)*2*np.pi
    K=np.hypot(*np.meshgrid(k1,k0)); dk=2*np.pi/(L/1e3); b=np.arange(dk/2,4,dk); i=np.digitize(K.ravel(),b)
    kr=np.array([K.ravel()[i==j].mean() for j in range(1,len(b))]); lp=np.array([np.log(P.ravel()[i==j].mean()) for j in range(1,len(b))])
    return kr,lp
for (cx,cy) in [(339e3,7909e3),(354e3,7889e3),(330e3,7840e3)]:
  try: kr,lp=spec_aero(cx,cy)
  except Exception as e: continue
  m=(kr>=0.3)&(kr<=3.0)
  res=lambda p: bouligand(kr[m],p[0],p[1],60.0,p[2])-lp[m]   # base fixa funda: insensível em k alto
  r=least_squares(res,[3,1.0,lp[m][0]],bounds=([0.5,0.01,-np.inf],[6,10,np.inf]))
  J=r.jac; cov=np.linalg.pinv(J.T@J)*(r.fun**2).sum()/(m.sum()-3); e=np.sqrt(np.diag(cov))
  print('aeromag janela centro (%d,%d) km: β=%.2f±%.2f  topo=%.2f±%.2f km (abaixo do plano de voo)  rms=%.2f'%(cx/1e3,cy/1e3,r.x[0],e[0],r.x[1],e[1],np.sqrt((r.fun**2).mean())))
