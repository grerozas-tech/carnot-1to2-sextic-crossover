import sys,numpy as np
from scipy.optimize import root_scalar, root
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,'/mnt/data/carnot_resonance_paper'); import numerics_double as nd

dx=1.04867934897e-6; mustar=nd.mu_star

def delta_robust(mu):
 y,T,DP=nd.mother_DP(mu); s=np.sqrt(max(0.,-DP[0,1]*DP[1,0])); c=-.5*np.trace(DP); return np.arctan2(s,c),y

def mu_for_d(d):
 r=root_scalar(lambda m:delta_robust(m)[0]-d,bracket=[mustar-d/250,mustar-d/400],xtol=1e-15,rtol=1e-14); return r.root

def sa(mu,y,z,n=8):
 xs,ws=leggauss(n); z=np.asarray(z,float); sm=0
 for x,w in zip((xs+1)/2,ws/2):
  F,DF,z1=nd.fmap(mu,y,x*z,with_dp=True); sm += w*(F[0]*(DF@z)[1]-(x*z[0])*z[1])
 return sm/2

prev_u=30.; prev=np.array([0.00057521,0.00670795]); prevprev_u=10.; prevprev=np.array([0.00045826,0.00523521])
for u in [40,50,60,70,80,90,100]:
 mu=mu_for_d(u*dx); d,y=delta_robust(mu)
 # linear predictor in u
 pred=prev + (prev-prevprev)*(u-prev_u)/(prev_u-prevprev_u)
 def fun(z): return nd.fmap(mu,y,z)[0]-z
 r=root(fun,pred,tol=1e-10)
 z=r.x; F,DF,z1=nd.fmap(mu,y,z,with_dp=True); tr=np.trace(DF); theta=np.arccos(np.clip(tr/2,-1,1)); J=sa(mu,y,z)
 print('u',u,'pred',pred,'z',z,'p',z[1],'J',J,'theta',theta,'tr',tr,'res',np.linalg.norm(F-z))
 prevprev_u,prevprev=prev_u,prev; prev_u,prev=u,z
