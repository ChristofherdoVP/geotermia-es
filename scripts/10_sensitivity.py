import numpy as np, pandas as pd
P=pd.read_csv('mc_parametros.csv')
A_c=1.77; A_med=1.99
def stats(w,label):
    w=w/w.sum(); ess=1/np.sum(w**2)
    def wq(x,q): i=np.argsort(x); c=np.cumsum(w[i]); return x[i][np.searchsorted(c,q)]
    qr=[wq(P.qr.values,q) for q in (.05,.5,.95)]; hr=wq(P.hr.values,.5)
    q0=P.qr+P.rho*A_med*P.hr
    Tm=P.T0+P.qr*P.zm/P.k+(P.rho*A_med*P.hr**2/P.k)*(1-np.exp(-P.zm/P.hr))
    # T at Moho > 580 ?
    pabove=np.sum(w*(Tm>580))
    print(f'{label:58s} ESS={ess:5.0f} | qr {qr[1]:.1f} ({qr[0]:.1f}-{qr[2]:.1f}) | hr {hr:.1f} | q0 {wq(q0.values,.5):.0f} ({wq(q0.values,.05):.0f}-{wq(q0.values,.95):.0f}) | T_Moho {wq(Tm.values,.5):.0f} ({wq(Tm.values,.05):.0f}-{wq(Tm.values,.95):.0f}) | P(T_Moho>580) {pabove:.3f}')
    return dict(case=label,ess=ess,qr50=qr[1],qr05=qr[0],qr95=qr[2],hr50=hr,q0_50=wq(q0.values,.5),q0_05=wq(q0.values,.05),q0_95=wq(q0.values,.95),TM50=wq(Tm.values,.5),TM05=wq(Tm.values,.05),TM95=wq(Tm.values,.95),p_above=pabove)
L=lambda pred,obs,s: np.exp(-0.5*((pred-obs)/s)**2)
qc=P.qr+P.rho*A_c*P.hr; qm=P.qr+P.rho*A_med*P.hr
rows=[]
rows.append(stats(np.ones(len(P)),'S0 prior (no heat-flow constraint)'))
rows.append(stats(L(qc,26,7.2),'S1 reference: Cachoeiro 26, sigma 7.2'))
rows.append(stats(L(qc,26,12),'S2 Cachoeiro 26, sigma inflated to 12'))
rows.append(stats(L(qc,34,7.2),'S3 Cachoeiro corrected +30% (34 mW/m2)'))
reg=np.array([26,28,33,29,31,39]); rm=reg.mean(); rs=reg.std(ddof=1)
print('regional onshore mean %.1f sd %.1f (n=%d)'%(rm,rs,len(reg)))
rows.append(stats(L(qm,rm,np.hypot(rs,6)),f'S4 regional mean {rm:.0f}±{rs:.0f} (6 sites) at median A0'))
rows.append(stats(L(qm,rm,np.hypot(rs,6))*L(qc,26,7.2),'S5 S1 and S4 combined'))
pd.DataFrame(rows).round(2).to_csv('sensitivity_heatflow.csv',index=False)
