from config import *
import os; os.makedirs('../figures', exist_ok=True)
import numpy as np, pandas as pd, geopandas as gpd, rasterio, harmonica as hm, xarray as xr, warnings; warnings.filterwarnings('ignore')
from rasterio import features
from pyproj import Transformer
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'font.family':'DejaVu Sans'})
P=DATA
es=gpd.read_file(P+'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp').to_crs(32724)
cor=gpd.read_file('lineamento_colatina_corredor.gpkg'); t=Transformer.from_crs(4326,32724,always_xy=True)
km=lambda g: g.geometry.scale(1e-3,1e-3,origin=(0,0))
def deco(a,lab=True):
    km(es).boundary.plot(ax=a,color='k',lw=.6); km(cor).boundary.plot(ax=a,color='c',lw=1.1,ls='--')
    a.set_xlim(195,450); a.set_ylim(7630,8030); a.set_aspect('equal'); a.set_xlabel('Easting (km)')
    if lab: a.set_ylabel('Northing (km)')
    a.spines[['top','right']].set_visible(False)
g=pd.read_csv('ghfdb_regiao.csv',low_memory=False)
g=g[g.environment.str.contains('onshore',na=False)&~g.name.astype(str).str.match(r'^Br\d+')]
gx,gy=t.transform(g.lon.values,g.lat.values); gx=np.array(gx)/1e3; gy=np.array(gy)/1e3
# ---------- Fig 1: geology
lit=gpd.read_file(P+'Espirito Santo/Litologia/EspiritoSanto_lito.shp').to_crs(32724)
def grupo(s):
    if s.startswith(('NP3_gamma_1','NP23_gamma_1')): return 'G1 pre-collisional'
    if s.startswith('NP3_gamma_2'): return 'G2 syn-collisional (S-type)'
    if 'gamma_5' in s: return 'G5 post-collisional'
    if s.startswith('PP'): return 'Paleoproterozoic basement'
    if s.startswith(('NP','MP')): return 'Neoproterozoic paragneiss complexes'
    if s.startswith(('ENb','Q','N','E')): return 'Cenozoic cover'
    return 'Other'
lit['grupo']=lit.SIGLA_UNID.astype(str).apply(grupo)
cols={'G1 pre-collisional':'#8fb3d9','G2 syn-collisional (S-type)':'#e8796b','G5 post-collisional':'#a35bb5','Neoproterozoic paragneiss complexes':'#c9b27c','Paleoproterozoic basement':'#7d8c6a','Cenozoic cover':'#f2ecc7','Other':'#cccccc'}
fig,a=plt.subplots(figsize=(8.6,8))
for k,c in cols.items(): km(lit[lit.grupo==k]).plot(ax=a,color=c,lw=0)
deco(a); sc=a.scatter(gx,gy,c='k',s=28,marker='^',zorder=5)
for x_,y_,q_,n_ in zip(gx,gy,g.q,g.name):
    if 195<x_<450 and 7660<y_<8030: a.annotate(f'{n_} ({q_:.0f})',(x_,y_),xytext=(4,-10),textcoords='offset points',fontsize=7)
xc,yc=t.transform(-40.63,-19.54); a.plot(xc/1e3,yc/1e3,'k*',ms=9); a.annotate('Colatina',(xc/1e3,yc/1e3),xytext=(5,4),textcoords='offset points',fontsize=8)
h=[Patch(fc=c,label=k) for k,c in cols.items()]+[plt.Line2D([],[],color='c',ls='--',label='Colatina lineament corridor'),plt.Line2D([],[],color='k',marker='^',ls='',label='Heat-flow site (GHFDB; mW m$^{-2}$)')]
a.legend(handles=h,loc='center left',bbox_to_anchor=(1.0,0.5),fontsize=7.5,frameon=False)
plt.tight_layout(); plt.savefig('../figures/Fig1_geology.png',dpi=300); plt.close()
# ---------- Fig 2: eTh and A0
with rasterio.open('A0_125m.tif') as d: A=d.read(1); tr=d.transform
ext=[tr.c/1e3,(tr.c+A.shape[1]*125)/1e3,(tr.f-A.shape[0]*125)/1e3,tr.f/1e3]
Th=hm.load_oasis_montaj_grid(P+'GM_TH_PPM.grd')
_x=Th.easting.values; _y=Th.northing.values; _d=_x[1]-_x[0]
_tr=rasterio.transform.from_origin(_x[0]-_d/2,_y[-1]+_d/2,_d,_d)
_in=features.geometry_mask(es.geometry,(len(_y),len(_x)),_tr,invert=True)[::-1]
Th=Th.where(_in)
fig,ax=plt.subplots(1,2,figsize=(10,7))
v=Th.values; im=ax[0].pcolormesh(Th.easting[::5]/1e3,Th.northing[::5]/1e3,v[::5,::5],cmap='viridis',vmin=0,vmax=np.nanpercentile(v,98),shading='auto',rasterized=True)
deco(ax[0]); fig.colorbar(im,ax=ax[0],shrink=.6,orientation='horizontal',pad=.08,label='eTh (ppm)'); ax[0].set_title('(a) Equivalent thorium')
im=ax[1].imshow(A[::2,::2],extent=ext,cmap='inferno',vmin=0,vmax=np.nanpercentile(A,99),interpolation='nearest')
deco(ax[1],False); ax[1].set_ylabel(''); fig.colorbar(im,ax=ax[1],shrink=.6,orientation='horizontal',pad=.08,label='A$_0$ (µW m$^{-3}$)'); ax[1].set_title('(b) Radiogenic heat production')
plt.tight_layout(); plt.savefig('../figures/Fig2_eTh_A0.png',dpi=300); plt.close()
# ---------- Fig 3: A0 by lithology
order=['G1 pre-collisional','G2 syn-collisional (S-type)','G5 post-collisional','Neoproterozoic paragneiss complexes','Paleoproterozoic basement','Cenozoic cover']
code={k:i+1 for i,k in enumerate(order)}
L=features.rasterize(((geo,code.get(gr,0)) for geo,gr in zip(lit.geometry,lit.grupo)),out_shape=A.shape,transform=tr,fill=0)
C=features.geometry_mask(cor.geometry,A.shape,tr,invert=True)
short={'G1 pre-collisional':'G1\npre-collisional','G2 syn-collisional (S-type)':'G2\nsyn-collisional','G5 post-collisional':'G5\npost-collisional','Neoproterozoic paragneiss complexes':'Paragneiss\ncomplexes','Paleoproterozoic basement':'Paleoproterozoic\nbasement','Cenozoic cover':'Cenozoic\ncover'}
fig,a=plt.subplots(figsize=(8,4)); pos=0; tk=[];tl=[]
for k in order:
    m=(L==code[k])&~np.isnan(A)
    for j,(mm,c) in enumerate([(m&C,'#e08a2c'),(m&~C,'#bdbdbd')]):
        if mm.sum()>200:
            b=a.boxplot([A[mm]],positions=[pos+j*.42],widths=.36,showfliers=False,patch_artist=True); b['boxes'][0].set_facecolor(c); b['medians'][0].set_color('k')
    tk.append(pos+.21); tl.append(short[k]); pos+=1.3
a.axhspan(1.6,4.0,color='0.9',zorder=0); a.text((tk[3]+tk[4])/2,8.6,'shaded band: fresh-rock range of eastern\nAraçuaí metasediments (Schannor et al., 2021)',fontsize=7,va='top',ha='center')
a.set_xticks(tk); a.set_xticklabels(tl,fontsize=8); a.set_ylabel('A$_0$ (µW m$^{-3}$)')
a.legend([Patch(fc='#e08a2c'),Patch(fc='#bdbdbd')],['inside Colatina corridor','outside corridor'],frameon=False,loc='upper left',fontsize=8)
a.spines[['top','right']].set_visible(False); plt.tight_layout(); plt.savefig('../figures/Fig3_A0_lithology.png',dpi=300); plt.close()
print('figs 1-3 ok')
