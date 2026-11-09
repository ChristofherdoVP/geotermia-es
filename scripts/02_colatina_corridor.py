from config import *
import geopandas as gpd, numpy as np, shapely
from shapely.geometry import LineString
P=DATA+'Espirito Santo/'
g=gpd.read_file(P+'Estrutural/EspiritoSanto_estrutura.shp').to_crs(32724)
g=g[g.geometry.notna()]
# eixo Ecoporanga -> Vitória (dissertação: 250 km, 30 km de largura, N10-15W)
ax=gpd.GeoSeries([LineString([(-40.83,-18.37),(-40.34,-20.32)])],crs=4326).to_crs(32724).iloc[0]
c=np.array(ax.coords); az=np.degrees(np.arctan2(c[0,0]-c[1,0],c[0,1]-c[1,1])); print('eixo: comprimento %.0f km, azimute N%.0fW'%(ax.length/1e3,-az))
corr=ax.buffer(15000,cap_style='flat')
def azim(geom):
  c=np.array(geom.coords if geom.geom_type=='LineString' else list(geom.geoms)[0].coords); d=c[-1]-c[0]
  a=np.degrees(np.arctan2(d[0],d[1]))%180; return a
g['az']=g.geometry.apply(azim); g['nnw']=g.az.between(150,175)  # N5W a N30W
inside=g[g.intersects(corr)]
print('feições no corredor:',len(inside),'| diques:',(inside.ESTRUTURA=='Dique').sum(),'| NNW (N5-30W):',inside.nnw.sum())
print(inside[inside.nnw].ESTRUTURA.value_counts().head(8).to_dict())
gpd.GeoDataFrame({'nome':['Lineamento Colatina (corredor 30 km, Ecoporanga–Vitória)']},geometry=[corr],crs=32724).to_file('lineamento_colatina_corredor.gpkg')
gpd.GeoDataFrame({'nome':['eixo']},geometry=[ax],crs=32724).to_file('lineamento_colatina_eixo.gpkg')
sel=inside[inside.nnw|(inside.ESTRUTURA=='Dique')]; sel.to_file('lineamento_colatina_estruturas.gpkg')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
es=gpd.read_file(P+'Limites_Administrativos/Limite_Estadual_es.shp').to_crs(32724)
f,a=plt.subplots(figsize=(6.5,8)); es.boundary.plot(ax=a,color='k',lw=.7)
g.plot(ax=a,color='0.8',lw=.4); gpd.GeoSeries([corr]).plot(ax=a,fc='C1',alpha=.2,ec='C1')
sel.plot(ax=a,color='C3',lw=.8); g[g.ESTRUTURA=='Dique'].plot(ax=a,color='C0',lw=1.2)
a.set_title('Lineamento Colatina: corredor da dissertação (laranja),\nestruturas NNW no corredor (vermelho), diques (azul)',fontsize=10)
a.set_axis_off(); plt.tight_layout(); plt.savefig('lineamento_colatina.png',dpi=150)
