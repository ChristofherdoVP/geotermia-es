from config import *
import os; exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'curie_emag.py'), encoding='utf-8').read().split('out={}')[0])
for L in [200e3,300e3,400e3]:
  kr,pr=spec(window(L),L)
  row=[]
  for beta in [0,1,2,3]:
    p=pr*kr**beta   # correção fractal (magnetização fractal, Maus & Dimri / Bouligand)
    b=kr<=0.20; a=(kr>=0.25)&(kr<=0.6)
    z0=-np.polyfit(kr[b],np.log(np.sqrt(p[b])/kr[b]),1)[0]; zt=-np.polyfit(kr[a],np.log(np.sqrt(p[a])),1)[0]
    row.append(2*z0-zt-4)  # profundidade abaixo da superfície (grade a 4 km)
  print(f'janela {int(L/1e3)} km: Zb (km abaixo da sup.) para beta=0,1,2,3:',np.round(row,1))
