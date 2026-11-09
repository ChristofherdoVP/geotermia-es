from config import *
import numpy as np, pandas as pd, geopandas as gpd, rasterio, warnings; warnings.filterwarnings('ignore')
from pyproj import Transformer
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.titlesize':10})
P=pd.read_csv('mc_parametros.csv'); rng=np.random.default_rng(1)
P_=DATA
es=gpd.read_file(P_+'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp').to_crs(32724)
cor=gpd.read_file('lineamento_colatina_corredor.gpkg'); t=Transformer.from_crs(4326,32724,always_xy=True)
km=lambda g: g.geometry.scale(1e-3,1e-3,origin=(0,0))
g=pd.read_csv('ghfdb_regiao.csv',low_memory=False); g=g[g.environment.str.contains('onshore',na=False)&~g.name.astype(str).str.match(r'^Br\d+')]
gx,gy=t.transform(g.lon.values,g.lat.values); gx=np.array(gx)/1e3; gy=np.array(gy)/1e3
fig=plt.figure(figsize=(12,6.2)); gs=fig.add_gridspec(1,3,width_ratios=[1,1,0.9])
for i,(f,tt,cm,lab,lim) in enumerate([('mc_q0_media_500m.tif','(a) Surface heat flow (posterior mean)','inferno','mW m$^{-2}$',(20,45)),('mc_TMoho_media_500m.tif','(b) Moho temperature (posterior mean)','magma','°C',None)]):
    a=fig.add_subplot(gs[i])
    with rasterio.open(f) as d: z=d.read(1); tr=d.transform
    ext=[tr.c/1e3,(tr.c+z.shape[1]*tr.a)/1e3,(tr.f+z.shape[0]*tr.e)/1e3,tr.f/1e3]; v=z[~np.isnan(z)]
    vmin,vmax=lim if lim else (np.percentile(v,2),np.percentile(v,98))
    im=a.imshow(z,extent=ext,cmap=cm,vmin=vmin,vmax=vmax,interpolation='nearest')
    km(es).boundary.plot(ax=a,color='k',lw=.6); km(cor).boundary.plot(ax=a,color='c',lw=1.1,ls='--')
    if i==0:
        a.scatter(gx,gy,c=g.q,cmap=cm,vmin=vmin,vmax=vmax,s=60,edgecolor='w',linewidth=1.2,zorder=5)
        for x_,y_,q_ in zip(gx,gy,g.q):
            if 195<x_<450 and 7660<y_<8030: a.annotate(f'{q_:.0f}',(x_,y_),xytext=(6,-3),textcoords='offset points',fontsize=8)
    a.set_xlim(195,450); a.set_ylim(7630,8030); a.set_aspect('equal'); a.set_xlabel('Easting (km)'); a.set_title(tt)
    a.set_ylabel('Northing (km)' if i==0 else '')
    fig.colorbar(im,ax=a,shrink=.7,orientation='horizontal',pad=.09,label=lab); a.spines[['top','right']].set_visible(False)
# (c) geotermas para A0 mediano do estado
with rasterio.open('A0_125m.tif') as d: A=d.read(1)
a0=np.nanmedian(A); z=np.linspace(0,100,201)
def T(p,a0): return p.T0+p.qr*z/p.k+(p.rho*a0*p.hr**2/p.k)*(1-np.exp(-z/p.hr))
a=fig.add_subplot(gs[2])
for nm,w,c in [('prior',None,'0.6'),('posterior',P.w.values,'C3')]:
    ids=rng.choice(len(P),400,p=w) if w is not None else rng.choice(len(P),400)
    TT=np.array([T(P.iloc[i],a0) for i in ids]); lo,med,hi=np.percentile(TT,[5,50,95],axis=0)
    a.fill_betweenx(z,lo,hi,color=c,alpha=.25,lw=0); a.plot(med,z,color=c,lw=1.5,label=f'{nm} (median, 5–95%)')
a.axhspan(33,40,color='C0',alpha=.15); a.text(40,36.5,'Moho (33–40 km)',fontsize=8,va='center',color='C0')
a.axvline(580,color='k',ls=':'); a.text(590,6,'580 °C (Curie,\nmagnetite)',fontsize=8,va='top')
a.set_ylim(100,0); a.set_xlim(0,1100); a.set_xlabel('Temperature (°C)'); a.set_ylabel('Depth (km)')
a.set_title(f'(c) Geotherms for median A$_0$ = {a0:.1f} µW m$^{{-3}}$'); a.legend(fontsize=8,frameon=False,loc='lower left'); a.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('../figures/Fig5_thermal.png',dpi=300)
# estatísticas para o texto
w=P.w.values
def wq(x,qs): i=np.argsort(x); c=np.cumsum(w[i]); return [x[i][np.searchsorted(c,q)] for q in qs]
for c in ['qr','hr','k','T0','zm']: print(c,'post 5/50/95:',np.round(wq(P[c].values,[.05,.5,.95]),2),' prior 5/50/95:',np.round(np.percentile(P[c],[5,50,95]),2))
ids=rng.choice(len(P),2000,p=w); 
Tm=np.array([P.T0[i]+P.qr[i]*P.zm[i]/P.k[i]+(P.rho[i]*a0*P.hr[i]**2/P.k[i])*(1-np.exp(-P.zm[i]/P.hr[i])) for i in ids])
q0=P.qr.values[ids]+P.rho.values[ids]*a0*P.hr.values[ids]
print('median-A0 column: q0 5/50/95',np.round(np.percentile(q0,[5,50,95]),1),' T Moho 5/50/95',np.round(np.percentile(Tm,[5,50,95]),0))
for f in ['mc_q0_media_500m.tif','mc_TMoho_media_500m.tif','mc_z580_P50_500m.tif']:
    with rasterio.open(f) as d: v=d.read(1); v=v[~np.isnan(v)]; print(f,'state 5/50/95',np.round(np.percentile(v,[5,50,95]),1))
