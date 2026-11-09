import numpy as np, xarray as xr, warnings; warnings.filterwarnings('ignore')
from scipy.special import gamma, kv
from scipy.optimize import least_squares
from pyproj import Transformer; from scipy.interpolate import griddata
H=4.0  # altura de observação do EMAG2v3 (km acima do geoide)
def bouligand(k,beta,zt,dz,C):
    x=k*dz
    A=np.sqrt(np.pi)/gamma(1+0.5*beta)*(0.5*np.cosh(x)*gamma(0.5*(1+beta))-kv(-0.5*(1+beta),x)*np.power(0.5*x,0.5*(1+beta)))
    return C-2*k*zt-(beta-1)*np.log(k)-x+np.log(A)
# grade UTM 2 km
da=xr.open_dataarray('emag2v3_recorte.nc'); LON,LAT=np.meshgrid(da.lon,da.lat)
t=Transformer.from_crs(4326,32724,always_xy=True); X,Y=t.transform(LON,LAT); ok=~np.isnan(da.values)
dx=2000.; gx=np.arange(-50e3,850e3,dx); gy=np.arange(7350e3,8450e3,dx)
GX,GY=np.meshgrid(gx,gy); G=griddata((X[ok],Y[ok]),da.values[ok],(GX,GY),method='cubic')
print('grade',G.shape,'NaN%',round(100*np.isnan(G).mean(),2))
def spectrum(cx,cy,L=300e3):
    ix=(gx>=cx-L/2)&(gx<cx+L/2); iy=(gy>=cy-L/2)&(gy<cy+L/2); w=G[np.ix_(iy,ix)]
    if np.isnan(w).mean()>0.02: return None
    w=np.where(np.isnan(w),np.nanmean(w),w); w=w-w.mean()
    # remove plano
    yy,xx=np.mgrid[:w.shape[0],:w.shape[1]]; A=np.c_[xx.ravel(),yy.ravel(),np.ones(xx.size)]; w=w-(A@np.linalg.lstsq(A,w.ravel(),rcond=None)[0]).reshape(w.shape)
    Pw=np.abs(np.fft.fft2(w*np.outer(np.hanning(w.shape[0]),np.hanning(w.shape[1]))))**2
    k1=np.fft.fftfreq(w.shape[1],dx/1e3)*2*np.pi; k0=np.fft.fftfreq(w.shape[0],dx/1e3)*2*np.pi
    K=np.hypot(*np.meshgrid(k1,k0)); dk=2*np.pi/(L/1e3); b=np.arange(dk/2,0.6+dk,dk); i=np.digitize(K.ravel(),b)
    kr=np.array([K.ravel()[i==j].mean() for j in range(1,len(b))]); pr=np.array([Pw.ravel()[i==j].mean() for j in range(1,len(b))])
    sd=np.array([np.log(Pw.ravel()[i==j]).std()/np.sqrt((i==j).sum()) for j in range(1,len(b))])
    return kr,np.log(pr),sd
def fit(kr,lp,sd,beta_prior=(3.0,1.0),beta_fixed=None):
    def res(p):
        b,zt,dz,C=p
        if beta_fixed is not None: b=beta_fixed
        r=(bouligand(kr,b,zt,dz,C)-lp)/np.maximum(sd,0.05)
        return np.r_[r,(b-beta_prior[0])/beta_prior[1]] if beta_fixed is None else r
    p0=[3.0,H+1,25.0,lp[0]]
    s=least_squares(res,p0,bounds=([0.5,0.5,1,-np.inf],[5,H+15,100,np.inf]))
    J=s.jac; dof=max(len(s.fun)-4,1); s2=(s.fun**2).sum()/dof
    try: cov=np.linalg.pinv(J.T@J)*s2; err=np.sqrt(np.diag(cov))
    except: err=np.full(4,np.nan)
    b,zt,dz,C=s.x; 
    if beta_fixed is not None: b=beta_fixed
    zb=zt+dz-H; ezb=np.sqrt(err[1]**2+err[2]**2+2*cov[1,2]) if np.all(np.isfinite(err)) else np.nan
    return b,zt-H,zb,ezb,err[0]
if __name__=='__main__':
    xc,yc=t.transform(-40.63,-19.54)
    kr,lp,sd=spectrum(xc,yc)
    print('Colatina, janela 300 km, %d pontos espectrais'%len(kr))
    b,zt,zb,ezb,eb=fit(kr,lp,sd); print('  ajuste conjunto (prior β=3±1): β=%.2f±%.2f  topo=%.1f km  Zb=%.1f±%.1f km'%(b,eb,zt,zb,ezb))
    for bf in [2,2.5,3,3.5,4]:
        b,zt,zb,ezb,_=fit(kr,lp,sd,beta_fixed=bf); print('  β fixo %.1f: topo=%.1f  Zb=%.1f±%.1f km'%(bf,zt,zb,ezb))
    for bp in [(2,1),(3,0.5),(4,1),(3,3)]:
        b,zt,zb,ezb,eb=fit(kr,lp,sd,beta_prior=bp); print('  prior β=%s: β=%.2f  Zb=%.1f±%.1f km'%(bp,b,zb,ezb))
