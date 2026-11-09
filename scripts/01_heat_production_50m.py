from config import *
import harmonica as hm, numpy as np, geopandas as gpd, warnings, rasterio
from rasterio import features
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
P=DATA
es=gpd.read_file(P+'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp').to_crs(32724)
K,Th,U=[hm.load_oasis_montaj_grid(P+f).values.astype('float32') for f in ['GM_K_Percetual.grd','GM_TH_PPM.grd','GM_U_PPM.grd']]
g=hm.load_oasis_montaj_grid(P+'GM_U_PPM.grd'); x=g.easting.values; y=g.northing.values; dx=x[1]-x[0]
tr=rasterio.transform.from_origin(x[0]-dx/2, y[-1]+dx/2, dx, dx)   # linhas de cima para baixo
inside=features.geometry_mask(es.geometry,out_shape=(len(y),len(x)),transform=tr,invert=True)[::-1]
rho=2700.; c={'U':9.52,'Th':2.56,'K':3.48}
A=1e-5*rho*(c['U']*U+c['Th']*Th+c['K']*K)
m=inside&~np.isnan(A); v=A[m]
aES=es.area.sum()/1e6
print('célula',dx,'m | área ES coberta',round(m.sum()*dx*dx/1e6),'km2 de',round(aES),'=',round(100*m.sum()*dx*dx/1e6/aES,1),'%')
print('A<0 %',round(100*(v<0).mean(),2))
print('percentis (5,25,50,75,95,99):',np.round(np.percentile(v,[5,25,50,75,95,99]),2),'média',round(v.mean(),2),'max',round(v.max(),1))
for n,a in [('U',U),('Th',Th),('K',K)]: print(' contribuição',n,round(100*(1e-5*rho*c[n]*a[m]).mean()/v.mean(),1),'%')
Am=np.where(m,A,np.nan)
with rasterio.open('calor_radiogenico_uWm3_ES_v2.tif','w',driver='GTiff',height=len(y),width=len(x),count=1,dtype='float32',crs='EPSG:32724',transform=tr,compress='deflate',nodata=np.nan) as d: d.write(Am[::-1],1)
s=4; Ad=Am[::s,::s]
fig,ax=plt.subplots(figsize=(8,10))
im=ax.pcolormesh(x[::s]/1e3,y[::s]/1e3,Ad,cmap='inferno',vmin=0,vmax=float(np.percentile(v,99)),shading='auto',rasterized=True)
es.geometry.boundary.scale(1e-3,1e-3,origin=(0,0)).plot(ax=ax,color='k',lw=.6)
ax.set_aspect('equal'); ax.set_xlabel('E (km, UTM 24S)'); ax.set_ylabel('N (km)')
ax.set_title('Produção de calor radiogênico (superfície), ES\nRybach (1988), ρ = 2700 kg/m³, Projeto 1093 (grids completos)')
cb=fig.colorbar(im,ax=ax,extend='max',shrink=.7); cb.set_label('A (µW/m³)')
ax.spines[['top','right']].set_visible(False); plt.tight_layout(); plt.savefig('calor_radiogenico_v2.png',dpi=150)
