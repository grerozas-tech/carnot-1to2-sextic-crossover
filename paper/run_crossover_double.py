import sys,csv,numpy as np
from scipy.optimize import root_scalar
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,'/mnt/data/carnot_resonance_paper'); import numerics_double as nd

dx=1.04867934897e-6
alpha=nd.alpha
A=203.67235802390137
mu_star=nd.mu_star

def delta_robust(mu):
    y,T,DP=nd.mother_DP(mu)
    s=np.sqrt(max(0.,-DP[0,1]*DP[1,0])); c=-0.5*np.trace(DP)
    return np.arctan2(s,c),y,T,DP

def mu_for_d(d):
    if d==0:return mu_star
    # linear bracket based on slope ~326.55
    a=mu_star-d/250.; b=mu_star-d/400.
    r=root_scalar(lambda m:delta_robust(m)[0]-d,bracket=[a,b],xtol=5e-16,rtol=1e-14,maxiter=50)
    return r.root

def section_action(mu,y,z,n=20):
    xs,ws=leggauss(n); ss=(xs+1)/2; ww=ws/2; z=np.asarray(z,float); sm=0.
    for s,w in zip(ss,ww):
        x=s*z; F,DF,z1=nd.fmap(mu,y,x,with_dp=True)
        sm += w*(F[0]*(DF@z)[1] - (s*z[0])*z[1])
    return sm/2

us=[0,.01,.03,.1,.3,1,3,10,30,100]
seed=np.array([2.5656e-4,2.8307e-3])
rows=[]
for u in us:
    dtarget=u*dx; mu=mu_for_d(dtarget); d,y,T,DPm=delta_robust(mu)
    z=nd.solve_fixed(mu,y,seed); seed=z.copy()
    F,DF,z1=nd.fmap(mu,y,z,with_dp=True)
    tr=np.trace(DF); det=np.linalg.det(DF); theta=np.arccos(np.clip(tr/2,-1,1)); J=section_action(mu,y,z,8)
    # linearized critical chart amplitude proxy only; not full I_NF
    Q=z[0]/np.sqrt(alpha); P=z[1]*np.sqrt(alpha); Ilin=(Q*Q+P*P)/2
    res=np.linalg.norm(F-z)
    print('u',u,'d',d,'mu',mu,'z',z,'J',J,'theta',theta,'2-tr',2-tr,'deterr',det-1,'Ilin',Ilin)
    rows.append([u,d,mu,z[0],z[1],z1[0],z1[1],Ilin,J,theta,tr,det-1,res])
path='/mnt/data/carnot_resonance_paper/data/crossover.csv'
with open(path,'w',newline='') as f:
    w=csv.writer(f);w.writerow(['u_target','delta_actual','mu','q','p','Pq','Pp','I_linear_proxy','J_reduced','theta_F','trace_DF','det_DF_error','fixed_point_residual']);w.writerows(rows)
print('WROTE',path)
