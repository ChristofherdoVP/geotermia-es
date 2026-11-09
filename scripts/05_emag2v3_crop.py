from config import *
import rasterio, numpy as np, xarray as xr
f=EMAG2_TIF
with rasterio.open(f) as d:
  w=rasterio.windows.from_bounds(360-47,-26,360-35,-13,d.transform); a=d.read(1,window=w).astype(float); tr=d.window_transform(w)
a[(a<-1e30)|(a>=99999)]=np.nan
ny,nx=a.shape; lon=tr.c+tr.a*(np.arange(nx)+.5)-360; lat=tr.f+tr.e*(np.arange(ny)+.5)
da=xr.DataArray(a,coords={'lat':lat,'lon':lon},dims=('lat','lon'),name='emag2v3_nT_4km')
print(da.shape,'lon',lon.min().round(2),lon.max().round(2),'lat',lat.min().round(2),lat.max().round(2))
print('NaN%',round(100*np.isnan(a).mean(),1),'pct 1/50/99',np.round(np.nanpercentile(a,[1,50,99]),1))
# valor perto de Colatina e comparação de extremos
print('Colatina',float(da.sel(lat=-19.54,lon=-40.63,method='nearest')))
da.to_netcdf('emag2v3_recorte.nc')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, geopandas as gpd
es=gpd.read_file(DATA+'Espirito Santo/Limites_Administrativos/Limite_Estadual_es.shp')
fig,ax=plt.subplots(figsize=(8,8))
v=np.nanpercentile(np.abs(a),98)
im=ax.pcolormesh(lon,lat,a,cmap='RdBu_r',vmin=-v,vmax=v,shading='auto')
es.boundary.plot(ax=ax,color='k',lw=.8); ax.plot(-40.63,-19.54,'k*',ms=10)
for L,c in [(200,'C2'),(300,'C1'),(400,'C3')]:
  h=L/2/111.; hx=h/np.cos(np.radians(19.54))
  ax.add_patch(plt.Rectangle((-40.63-hx,-19.54-h),2*hx,2*h,fill=False,ec=c,lw=1.5,label=f'janela {L} km'))
ax.set_aspect(1/np.cos(np.radians(19.5))); ax.legend(loc='lower right',frameon=False)
fig.colorbar(im,ax=ax,shrink=.75,label='anomalia (nT), 4 km'); ax.set_title('EMAG2v3 (Meyer et al. 2017), recorte ES/MG; ★ Colatina')
plt.tight_layout(); plt.savefig('emag2v3_recorte.png',dpi=130)
