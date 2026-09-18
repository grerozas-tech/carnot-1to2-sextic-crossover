import sys,numpy as np
from scipy.optimize import root_scalar
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,'/mnt/data/carnot_resonance_paper'); import numerics_double as nd
DX=1.04867934897e-6; mustar=nd.mu_star
z0=np.array([0.00025656,0.00283072])
def drob(mu):
 y,T,DP=nd.mother_DP(mu); s=np.sqrt(max(0.,-DP[0,1]*DP[1,0])); c=-.5*np.trace(DP); return np.arctan2(s,c),y
def mfor(u):
 d=u*DX
 if d==0:return mustar
 return root_scalar(lambda m:drob(m)[0]-d,bracket=[mustar-d/250,mustar-d/400],rtol=1e-14).root
def sa(mu,y,z,n=8):
 xs,ws=leggauss(n); sm=0.; z=np.asarray(z,float)
 for s,w in zip((xs+1)/2,ws/2):
  F,DF,z1=nd.fmap(mu,y,s*z,with_dp=True); sm+=w*(F[0]*(DF@z)[1]-(s*z[0])*z[1])
 return sm/2
def one(u):
 mu=mfor(u); d,y=drob(mu); R=(1+np.sqrt(1+3*u))/2; seed=z0*np.sqrt(R)
 z=nd.solve_fixed(mu,y,seed); F,DF,z1=nd.fmap(mu,y,z,with_dp=True); tr=np.trace(DF); th=np.arccos(np.clip(tr/2,-1,1)); J=sa(mu,y,z)
 return [u,d,mu,*z,*z1,J,th,tr,np.linalg.det(DF)-1,np.linalg.norm(F-z)]
if __name__=='__main__':
 import sys
 for u in map(float,sys.argv[1:]):
  r=one(u); print(','.join(repr(x) for x in r),flush=True)
