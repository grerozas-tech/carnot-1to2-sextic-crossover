import json,mpmath as mp,pandas as pd
from pathlib import Path
BASE=Path('/mnt/data/carnot_resonance_paper')
mp.mp.dps=50
j=json.load(open(BASE/'data'/'critical_jet_transport_60dps.json'))
def load(d): return {tuple(map(int,k.split(','))):mp.mpf(v) for k,v in d.items()}
Nmap=[load(j['Nmap_nonresonant'][0]),load(j['Nmap_nonresonant'][1])]
alpha=mp.mpf(j['alpha_used']); sa=mp.sqrt(alpha)
def ev(poly,q,p): return sum(v*q**i*p**jj for (i,jj),v in poly.items())
def nf(q,p):
 q=mp.mpf(str(q));p=mp.mpf(str(p)); qr=ev(Nmap[0],q,p); pr=ev(Nmap[1],q,p)
 Q=qr/sa; P=pr*sa; I=(Q*Q+P*P)/2
 return Q,P,I,qr,pr
p=BASE/'data'/'crossover.csv'; df=pd.read_csv(p)
rows=[]
for _,r in df.iterrows():
 Q,P,I,qr,pr=nf(r.q,r.p); Qc,Pc,Ic,qrc,prc=nf(r.Pq,r.Pp)
 odd=mp.sqrt((Q+Qc)**2+(P+Pc)**2)
 rows.append((float(Q),float(P),float(I),float(Qc),float(Pc),float(Ic),float(odd)))
cols=['Q_NF','P_NF','I_NF','Qcomp_NF','Pcomp_NF','Icomp_NF','NF_odd_pair_residual']
for k,c in enumerate(cols): df[c]=[x[k] for x in rows]
A=mp.mpf('203.67235802390137'); rho=mp.mpf('2.0011829182526856'); sigma=mp.mpf('-0.0047330723063350145'); G=mp.mpf('886147.1053518064')
a=A*sigma; delta_x=a*a/G; I0=2*abs(a)/(3*G); J0=-4*abs(a)**3/(27*G**2); theta0=8*rho*I0*mp.sqrt(A*abs(a))
df['u_actual']=df['delta_actual']/float(delta_x)
import numpy as np
u=df['u_actual'].to_numpy(); S=np.sqrt(1+3*u); R=(1+S)/2
df['X_NF_norm']=df['I_NF']/float(I0); df['X_scalar']=R
df['Y_exact_norm']=df['J_reduced']/float(J0); df['Y_scalar']=4*R**3-3*R**2
df['Z_exact_norm']=df['theta_F']/float(theta0); df['Z_scalar']=R*np.sqrt(S)
df['X_rel_dev']=(df['X_NF_norm']-df['X_scalar'])/df['X_scalar']
df['Y_rel_dev']=(df['Y_exact_norm']-df['Y_scalar'])/df['Y_scalar']
df['Z_rel_dev']=(df['Z_exact_norm']-df['Z_scalar'])/df['Z_scalar']
df['I_pair_rel_diff']=abs(df['I_NF']-df['Icomp_NF'])/df['I_NF']
df.to_csv(p,index=False)
print(df[['u_target','I_NF','X_NF_norm','X_scalar','I_pair_rel_diff','NF_odd_pair_residual']].to_string(index=False))
