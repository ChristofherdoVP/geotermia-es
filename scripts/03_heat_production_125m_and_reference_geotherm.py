from config import *
import numpy as np, rasterio, geopandas as gpd, warnings; warnings.filterwarnings('ignore')
from rasterio.warp import reproject, Resampling
from scipy.ndimage import gaussian_filter
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
# 1) A0 em 125 m (média por área a partir da grade de 50 m)
with rasterio.open('calor_radiogenico_uWm3_ES_v2.tif') as s:
    tr=rasterio.transform.from_origin(s.bounds.left,s.bounds.top,125,125)
    W=int(np.ceil((s.bounds.right-s.bounds.left)/125)); Hh=int(np.ceil((s.bounds.top-s.bounds.bottom)/125))
    A=np.full((Hh,W),np.nan,'float32')
    reproject(rasterio.band(s,1),A,src_transform=s.transform,src_crs=s.crs,dst_transform=tr,dst_crs=s.crs,resampling=Resampling.average,src_nodata=np.nan,dst_nodata=np.nan)
prof=dict(driver='GTiff',height=Hh,width=W,count=1,dtype='float32',crs='EPSG:32724',transform=tr,compress='deflate',nodata=np.nan)
with rasterio.open('A0_125m.tif','w',**prof) as d: d.write(A,1)
m=~np.isnan(A); print('A0 125 m: mediana %.2f  média %.2f µW/m³'%(np.median(A[m]),A[m].mean()))
# 2) suavização lateral (σ = 5 km) — o calor de corpos pequenos não controla a geoterma regional
def smooth(a,sig_km):
    s=sig_km*1000/125; w=gaussian_filter(np.where(m,1.,0.),s); v=gaussian_filter(np.where(m,a,0.),s); return np.where(m,v/np.maximum(w,1e-6),np.nan)
As=smooth(A,5.0)
# 3) modelo condutivo 1D com A(z)=A0 exp(-z/hr)
def model(A0,qr,k,hr,T0=23.0):
    q0=qr+A0*hr                          # mW/m² (µW/m³ * km = mW/m²)
    T=lambda z: T0+qr*z/k+(A0*hr**2/k)*(1-np.exp(-z/hr))   # °C, z em km, q em mW/m², k em W/m/K
    return q0,T
def z_of_T(A0,qr,k,hr,Tc=580.,T0=23.):
    z=np.full(A0.shape,0.0); lo=np.zeros_like(A0); hi=np.full_like(A0,120.)
    for _ in range(40):
        mid=(lo+hi)/2; T=T0+qr*mid/k+(A0*hr**2/k)*(1-np.exp(-mid/hr)); hi=np.where(T>Tc,mid,hi); lo=np.where(T>Tc,lo,mid)
    return (lo+hi)/2
base=dict(qr=35.,k=2.8,hr=10.)
_z=z_of_T
q0,T=model(As,**base); z580=np.where(m,z_of_T(As,**base),np.nan); T5=T(5.); TM=T(36.)
def st(x): v=x[m]; return '%.1f [%.1f–%.1f]'%(np.median(v),np.percentile(v,5),np.percentile(v,95))
print('CASO BASE qr=35 mW/m², k=2.8 W/m/K, hr=10 km, T0=23 °C')
print('  q0 (mW/m²):',st(q0)); print('  gradiente sup. (°C/km):',st(q0/base['k'])); print('  T a 5 km (°C):',st(T5)); print('  T no Moho 36 km (°C):',st(TM)); print('  prof. isoterma 580 °C (km):',st(z580))
print('SENSIBILIDADE (mediana da prof. 580 °C | q0 mediano):')
for qr in [30,35,40]:
  for k in [2.4,2.8,3.2]:
    for hr in [8,10,12]:
      if (qr!=35)+(k!=2.8)+(hr!=10)>1: continue
      zz=z_of_T(As,qr,k,hr); qq=qr+As*hr; print('  qr=%d k=%.1f hr=%d: z580=%.1f km  q0=%.1f'%(qr,k,hr,np.median(zz[m]),np.median(qq[m])))
for nm,a in [('q0_mWm2',q0),('T5km_C',T5),('z580_km',z580),('TMoho36_C',TM)]:
    with rasterio.open(f'geoterma_{nm}_125m.tif','w',**prof) as d: d.write(a.astype('float32'),1)
# figura
es=gpd.read_file(DATA+'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp').to_crs(32724)
cor=gpd.read_file('lineamento_colatina_corredor.gpkg')
x0,y1=tr.c,tr.f; ext=[x0/1e3,(x0+W*125)/1e3,(y1-Hh*125)/1e3,y1/1e3]
fig,ax=plt.subplots(1,3,figsize=(15,7.5))
for a,(z,tt,cm,lab) in zip(ax,[(q0,'Fluxo de calor superficial','inferno','mW/m²'),(T5,'Temperatura a 5 km','magma','°C'),(z580,'Profundidade da isoterma de 580 °C','viridis_r','km')]):
    v=z[m]; im=a.imshow(z[::4,::4],extent=ext,cmap=cm,vmin=np.percentile(v,2),vmax=np.percentile(v,98),interpolation='nearest')
    es.geometry.boundary.scale(1e-3,1e-3,origin=(0,0)).plot(ax=a,color='k',lw=.6)
    cor.geometry.boundary.scale(1e-3,1e-3,origin=(0,0)).plot(ax=a,color='w',lw=1,ls='--')
    a.set_title(tt,fontsize=11); a.set_aspect('equal'); a.set_xlabel('E (km, UTM 24S)'); fig.colorbar(im,ax=a,shrink=.6,orientation='horizontal',pad=.08,label=lab)
    a.spines[['top','right']].set_visible(False)
ax[0].set_ylabel('N (km)'); [a.set_ylabel('') for a in ax[1:]]
fig.suptitle('Modelo condutivo 1D (qr = 35 mW/m², k = 2,8 W/m·K, hr = 10 km, T0 = 23 °C; A0 suavizado σ = 5 km)  —  tracejado: corredor do Lineamento Colatina',fontsize=10)
plt.tight_layout(); plt.savefig('geoterma_v1.png',dpi=140)
