from config import *
import numpy as np, rasterio, pandas as pd, geopandas as gpd, warnings; warnings.filterwarnings('ignore')
from scipy.ndimage import gaussian_filter
from pyproj import Transformer
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rng=np.random.default_rng(42)
with rasterio.open('A0_125m.tif') as d: A=d.read(1); tr=d.transform
m=~np.isnan(A); s=40; w=gaussian_filter(m*1.,s); As=np.where(m,gaussian_filter(np.where(m,A,0),s)/np.maximum(w,1e-6),np.nan)
f=4; Ad=As[::f,::f]; md=~np.isnan(Ad); a=Ad[md]            # 500 m para o Monte Carlo
t=Transformer.from_crs(4326,32724,always_xy=True)
xc,yc=t.transform(-41.10,-20.85); ic,jc=~tr*(xc,yc); A_cach=As[int(jc),int(ic)]
print('A0 suavizado em Cachoeiro: %.2f µW/m³'%A_cach)
N=4000
P=pd.DataFrame(dict(qr=rng.uniform(10,40,N),hr=rng.uniform(2,12,N),k=rng.uniform(2.2,3.6,N),rho=rng.normal(1,0.04,N),T0=rng.uniform(18,25,N),zm=rng.uniform(33,40,N)))
P['q_cach']=P.qr+P.rho*A_cach*P.hr
sig=np.hypot(4,6)   # incerteza da medida (4) + representatividade/qualidade (6)
P['w']=np.exp(-0.5*((P.q_cach-26)/sig)**2); P['w']/=P.w.sum()
ess=1/np.sum(P.w**2); print('amostras efetivas (posterior): %.0f de %d'%(ess,N))
def wq(x,w,qs): i=np.argsort(x); c=np.cumsum(w[i]); return [x[i][np.searchsorted(c,q)] for q in qs]
for c in ['qr','hr','k']:
    print('  %s: prior mediana %.1f | posterior mediana %.1f [%.1f–%.1f]'%(c,np.median(P[c]),*wq(P[c].values,P.w.values,[.5,.05,.95])))
def z580(qr,k,hr,a,T0,Tc=580):
    lo=np.zeros_like(a); hi=np.full_like(a,150.)
    for _ in range(30):
        mid=(lo+hi)/2; T=T0+qr*mid/k+(a*hr**2/k)*(1-np.exp(-mid/hr)); hi=np.where(T>Tc,mid,hi); lo=np.where(T>Tc,lo,mid)
    return (lo+hi)/2
idx=rng.choice(N,300,p=P.w.values)          # reamostragem da posterior
idp=rng.choice(N,300)                        # prior
def run(ids):
    Q=[];Z=[];TM=[]
    for i in ids:
        p=P.iloc[i]; aa=a*p.rho
        Q.append(p.qr+aa*p.hr); Z.append(z580(p.qr,p.k,p.hr,aa,p.T0)); TM.append(p.T0+p.qr*p.zm/p.k+(aa*p.hr**2/p.k)*(1-np.exp(-p.zm/p.hr)))
    return [np.array(x) for x in (Q,Z,TM)]
res={}
for nm,ids in [('prior',idp),('posterior',idx)]:
    Q,Z,TM=run(ids); res[nm]=(Q,Z,TM)
    print(nm.upper(),': q0 %.0f±%.0f mW/m² | grad %.1f °C/km | T Moho %.0f±%.0f °C | z580 mediana %.0f km (P(z580<36 km)=%.2f)'%(
        np.median(Q),np.median(Q.std(0)),np.median(Q/np.array([P.k[i] for i in (idx if nm=='posterior' else idp)])[:,None]),np.median(TM),np.median(TM.std(0)),np.median(Z),np.mean(Z<np.array([P.zm[i] for i in (idx if nm=='posterior' else idp)])[:,None])))
Q,Z,TM=res['posterior']
def tomap(v): o=np.full(Ad.shape,np.nan); o[md]=v; return o
maps={'q0_media':tomap(Q.mean(0)),'q0_dp':tomap(Q.std(0)),'TMoho_media':tomap(TM.mean(0)),'z580_P50':tomap(np.median(Z,0))}
tr4=rasterio.Affine(tr.a*f,0,tr.c,0,tr.e*f,tr.f)
for k_,v in maps.items():
    with rasterio.open(f'mc_{k_}_500m.tif','w',driver='GTiff',height=v.shape[0],width=v.shape[1],count=1,dtype='float32',crs='EPSG:32724',transform=tr4,nodata=np.nan,compress='deflate') as d: d.write(v.astype('float32'),1)
P.to_csv('mc_parametros.csv',index=False)
# figura
g=pd.read_csv('ghfdb_regiao.csv',low_memory=False)
g=g[g.environment.str.contains('onshore',na=False)&~g.name.astype(str).str.match(r'^Br\d+')]
gx,gy=t.transform(g.lon.values,g.lat.values)
es=gpd.read_file(DATA+'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp').to_crs(32724)
cor=gpd.read_file('lineamento_colatina_corredor.gpkg')
ext=[tr.c/1e3,(tr.c+Ad.shape[1]*tr.a*f)/1e3,(tr.f+Ad.shape[0]*tr.e*f)/1e3,tr.f/1e3]
fig,ax=plt.subplots(1,3,figsize=(16,7.5))
spec=[(maps['q0_media'],'Fluxo de calor (média posterior)','inferno','mW/m²',(20,45)),
      (maps['q0_dp'],'Incerteza do fluxo (1σ)','Greys','mW/m²',None),
      (maps['TMoho_media'],'Temperatura no Moho (33–40 km)','magma','°C',None)]
for i,(a_,(z,tt,cm,lab,lim)) in enumerate(zip(ax,spec)):
    v=z[~np.isnan(z)]; vmin,vmax=lim if lim else (np.percentile(v,2),np.percentile(v,98))
    im=a_.imshow(z,extent=ext,cmap=cm,vmin=vmin,vmax=vmax,interpolation='nearest')
    es.geometry.boundary.scale(1e-3,1e-3,origin=(0,0)).plot(ax=a_,color='k',lw=.6)
    cor.geometry.boundary.scale(1e-3,1e-3,origin=(0,0)).plot(ax=a_,color='c',lw=1,ls='--')
    if i==0:
        sc=a_.scatter(np.array(gx)/1e3,np.array(gy)/1e3,c=g.q,cmap='inferno',vmin=vmin,vmax=vmax,s=70,edgecolor='w',linewidth=1.2,zorder=5)
        for x_,y_,q_,n_ in zip(np.array(gx)/1e3,np.array(gy)/1e3,g.q,g.name):
            if 150<x_<470 and 7680<y_<8050: a_.annotate(f'{n_}\n{q_:.0f}',(x_,y_),xytext=(5,-12),textcoords='offset points',fontsize=7,color='k')
    a_.set_xlim(150,470); a_.set_ylim(7600,8050); a_.set_aspect('equal'); a_.set_title(tt,fontsize=11); a_.set_xlabel('E (km, UTM 24S)')
    fig.colorbar(im,ax=a_,shrink=.6,orientation='horizontal',pad=.08,label=lab); a_.spines[['top','right']].set_visible(False)
ax[0].set_ylabel('N (km)'); [a.set_ylabel('') for a in ax[1:]]
fig.suptitle('Monte Carlo condicionado ao fluxo medido em Cachoeiro de Itapemirim (26 ± 4 mW/m²; GHFDB 2024)  •  círculos: medidas do GHFDB  •  tracejado: corredor do Colatina',fontsize=10)
plt.tight_layout(); plt.savefig('montecarlo_termico.png',dpi=140)
