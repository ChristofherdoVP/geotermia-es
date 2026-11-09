from config import *
import geopandas as gpd, numpy as np, pandas as pd, rasterio, warnings; warnings.filterwarnings('ignore')
from rasterio import features
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
g=gpd.read_file(DATA+'Espirito Santo/Litologia/EspiritoSanto_lito.shp').to_crs(32724)
def grupo(s):
    if s.startswith('NP3_gamma_1') or s.startswith('NP23_gamma_1'): return 'G1 pré-colisional'
    if s.startswith('NP3_gamma_2'): return 'G2 sin-colisional'
    if 'gamma_5' in s: return 'G5 pós-colisional'
    if s.startswith('PP'): return 'Embasamento paleoproterozoico'
    if s.startswith(('NP','MP')): return 'Complexos paraderivados NP'
    if s.startswith(('ENb','Q','N','E')): return 'Coberturas cenozoicas'
    return 'Outros'
g['grupo']=g.SIGLA_UNID.astype(str).apply(grupo)
print(g.assign(a=g.area/1e6).groupby('grupo').a.sum().round(0).to_string())
order=['G1 pré-colisional','G2 sin-colisional','G5 pós-colisional','Complexos paraderivados NP','Embasamento paleoproterozoico','Coberturas cenozoicas','Outros']
code={k:i+1 for i,k in enumerate(order)}
with rasterio.open('A0_125m.tif') as d: A=d.read(1); tr=d.transform
L=features.rasterize(((geom,code[gr]) for geom,gr in zip(g.geometry,g.grupo)),out_shape=A.shape,transform=tr,fill=0)
cor=gpd.read_file('lineamento_colatina_corredor.gpkg'); C=features.geometry_mask(cor.geometry,A.shape,tr,invert=True)
rows=[]
for k,i in code.items():
    m=(L==i)&~np.isnan(A)
    if m.sum()<100: continue
    a,ai,ao=A[m],A[m&C],A[m&~C]
    rows.append(dict(grupo=k,area_km2=round(m.sum()*0.015625),mediana=np.median(a),p25=np.percentile(a,25),p75=np.percentile(a,75),
        corr_med=np.median(ai) if ai.size>200 else np.nan, fora_med=np.median(ao) if ao.size>200 else np.nan, frac_no_corredor=round(100*ai.size/a.size,1)))
T=pd.DataFrame(rows); T['dif_%']=100*(T.corr_med/T.fora_med-1)
pd.set_option('display.width',200); print(T.round(2).to_string(index=False))
T.round(3).to_csv('A0_por_litologia.csv',index=False)
# teste de permutação em blocos (dependência espacial): corredor vs fora dentro de G2 e paraderivados
rng=np.random.default_rng(0)
for k in ['G1 pré-colisional','G2 sin-colisional','Complexos paraderivados NP']:
    m=(L==code[k])&~np.isnan(A)
    # amostra em blocos de 10 km (80 células): médias de bloco
    b=80; H,W=A.shape; vals=[];inc=[]
    for r in range(0,H-b,b):
        for c in range(0,W-b,b):
            mm=m[r:r+b,c:c+b]
            if mm.sum()>b*b*0.3: vals.append(np.median(A[r:r+b,c:c+b][mm])); inc.append(C[r:r+b,c:c+b][mm].mean()>0.5)
    vals=np.array(vals); inc=np.array(inc)
    if inc.sum()<3: continue
    obs=np.median(vals[inc])-np.median(vals[~inc]); null=[]
    for _ in range(5000): p=rng.permutation(inc); null.append(np.median(vals[p])-np.median(vals[~p]))
    pv=np.mean(np.array(null)>=obs)
    print(f'{k}: blocos 10 km, n_corr={inc.sum()} n_fora={(~inc).sum()}, Δmediana={obs:.2f} µW/m³, p (permutação, unilateral)={pv:.3f}')
# boxplot em pares
short={'G1 pré-colisional':'G1\npré-colisional','G2 sin-colisional':'G2\nsin-colisional','G5 pós-colisional':'G5\npós-colisional','Complexos paraderivados NP':'Paraderivados\nneoproterozoicos','Embasamento paleoproterozoico':'Embasamento\npaleoproterozoico','Coberturas cenozoicas':'Coberturas\ncenozoicas','Outros':'Outros'}
fig,ax=plt.subplots(figsize=(11,5)); pos=0; ticks=[];tl=[]
for k in order:
    m=(L==code[k])&~np.isnan(A)
    if m.sum()<100: continue
    for j,(mm,c) in enumerate([(m&C,'C1'),(m&~C,'0.75')]):
        if mm.sum()>200:
            b=ax.boxplot([A[mm]],positions=[pos+j*0.42],widths=.36,showfliers=False,patch_artist=True)
            b['boxes'][0].set_facecolor(c); b['medians'][0].set_color('k')
    ticks.append(pos+0.21); tl.append(short[k]); pos+=1.3
ax.set_xticks(ticks); ax.set_xticklabels(tl,fontsize=9); ax.set_ylabel('A0 (µW/m³)')
from matplotlib.patches import Patch
ax.legend([Patch(fc='C1'),Patch(fc='0.75')],['dentro do corredor do Colatina','fora do corredor'],frameon=False,loc='upper right')
ax.set_title('Produção de calor radiogênico por grupo litológico (grade 125 m)',fontsize=11)
ax.spines[['top','right']].set_visible(False); plt.tight_layout(); plt.savefig('A0_por_litologia.png',dpi=150)
