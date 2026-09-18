import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, root_scalar

A=203.67235802390137
alpha=0.0946271414514154
mu_star=0.6628690069305101


def f(t,X,mu):
    x,y,z=X
    return np.array([np.cos(z),np.sin(z),-0.5*x*x+0.5*y*y-mu])

def jac3(X):
    x,y,z=X
    return np.array([[0,0,-np.sin(z)],[0,0,np.cos(z)],[-x,y,0.]])

def shoot_res(v,mu):
    a,T=v
    sol=solve_ivp(lambda t,X:f(t,X,mu),(0,T),(0,a,0),method='DOP853',rtol=2e-12,atol=2e-14,max_step=.12)
    X=sol.y[:,-1]
    return np.array([X[1],X[2]-np.pi/2])

def mother(mu,guess=None):
    if guess is None:
        dm=mu-mu_star
        guess=(-2.57258930747695 + 2.2*dm, 3.01795545336086 + 7.5*dm)
    r=root(lambda v:shoot_res(v,mu),guess,tol=1e-11)
    if np.linalg.norm(r.fun)>1e-9: raise RuntimeError(('mother fail',mu,r.fun,r.message))
    return r.x

def mother_DP(mu,guess=None):
    a,Tq=mother(mu,guess)
    Y0=np.r_[[0.,a,0.],np.eye(3).ravel()]
    def rhs(t,Y):
        X=Y[:3]; P=Y[3:].reshape(3,3)
        return np.r_[f(t,X,mu),(jac3(X)@P).ravel()]
    sol=solve_ivp(rhs,(0,4*Tq),Y0,method='DOP853',rtol=2e-12,atol=2e-14,max_step=.08)
    Y=sol.y[:,-1]; X=Y[:3]; Phi=Y[3:].reshape(3,3)
    ff=f(4*Tq,X,mu); Proj=np.eye(3)-np.outer(ff,[1.,0,0])/ff[0]
    DE=np.array([[0,0],[1,0],[0,1.]])
    DC=np.array([[0,1,0],[0,0,np.cos(X[2])]])
    DP=DC@Proj@Phi@DE
    return a,Tq,DP

def delta_of_mu(mu):
    a,T,DP=mother_DP(mu)
    c=np.clip(-np.trace(DP)/2,-1,1)
    return np.arccos(c),a,T,DP

def mu_for_delta(delta,side='below'):
    # solve using b/alpha as robust small-delta signed proxy, otherwise exact acos
    if delta==0: return mu_star
    # k around 327; broaden bracket
    lo=mu_star-1e-3; hi=mu_star-1e-12 if side=='below' else mu_star+1e-12
    if side!='below': lo,hi=mu_star+1e-12,mu_star+1e-3
    def fun(mu): return delta_of_mu(mu)[0]-delta
    # find bracket adaptively near star
    if side=='below':
        a=mu_star-1e-10; b=mu_star-1e-3
    else:
        a=mu_star+1e-10; b=mu_star+1e-3
    fa=fun(a); fb=fun(b)
    if fa*fb>0: raise RuntimeError(('no bracket',delta,fa,fb))
    r=root_scalar(fun,bracket=sorted([a,b]),xtol=2e-13,rtol=2e-13,maxiter=40)
    return r.root

def _event():
    def ev(t,Y): return Y[0]
    ev.direction=1; ev.terminal=False
    return ev

def pmap(mu,ystar,z0,tmax=15.0,radius=1.0,with_dp=False,with_action=False):
    q,p=z0
    X0=np.array([0.,ystar+q,np.arcsin(p)])
    if with_dp or with_action:
        parts=[X0]
        if with_dp: parts.append(np.eye(3).ravel())
        if with_action: parts.append(np.array([0.]))
        Y0=np.concatenate(parts)
        def rhs(t,Y):
            X=Y[:3]; out=[*f(t,X,mu)]
            off=3
            if with_dp:
                Phi=Y[off:off+9].reshape(3,3); out.extend((jac3(X)@Phi).ravel()); off+=9
            if with_action:
                x,y,z=X
                af=-1 + x*np.sin(z)*(0.5*y*y-mu-x*x/6)
                out.append(af)
            return np.array(out)
    else:
        Y0=X0
        rhs=lambda t,Y:f(t,Y,mu)
    sol=solve_ivp(rhs,(0,tmax),Y0,events=_event(),method='DOP853',rtol=2e-12,atol=2e-14,max_step=.06)
    for t,Y in zip(sol.t_events[0],sol.y_events[0]):
        if t<1e-6: continue
        X=Y[:3]
        out=np.array([X[1]-ystar,np.sin(X[2])])
        dist=np.linalg.norm(out)
        if np.cos(X[2])>0 and dist<radius:
            ans=[out]
            off=3
            if with_dp:
                Phi=Y[off:off+9].reshape(3,3); off+=9
                ff=f(t,X,mu); Proj=np.eye(3)-np.outer(ff,[1.,0,0])/ff[0]
                DE=np.array([[0,0],[1,0],[0,1/np.sqrt(1-p*p)]])
                DC=np.array([[0,1,0],[0,0,np.cos(X[2])]])
                DP=DC@Proj@Phi@DE
                ans.append(DP)
            if with_action: ans.append(Y[off])
            ans.extend([t,dist,X])
            return tuple(ans)
    raise RuntimeError(('no local return',mu,z0,[(t, np.linalg.norm([Y[1]-ystar,np.sin(Y[2])])) for t,Y in zip(sol.t_events[0],sol.y_events[0])]))

def fmap(mu,ystar,z,with_dp=False,with_action=False):
    r1=pmap(mu,ystar,z,with_dp=with_dp,with_action=with_action)
    z1=r1[0]; idx=1
    DP1=None; A1=None
    if with_dp: DP1=r1[idx]; idx+=1
    if with_action: A1=r1[idx]; idx+=1
    r2=pmap(mu,ystar,z1,with_dp=with_dp,with_action=with_action)
    z2=r2[0]; idx=1
    DP2=None; A2=None
    if with_dp: DP2=r2[idx]; idx+=1
    if with_action: A2=r2[idx]; idx+=1
    ans=[z2]
    if with_dp: ans.append(DP2@DP1)
    if with_action: ans.append(A1+A2)
    ans.append(z1)
    return tuple(ans)

def solve_fixed(mu,ystar,seed):
    def fun(z): return fmap(mu,ystar,z)[0]-z
    r=root(fun,seed,tol=1e-10)
    if np.linalg.norm(r.fun)>1e-8: raise RuntimeError(('fixed fail',mu,seed,r.x,r.fun,r.message))
    return r.x

if __name__=='__main__':
    for delta in [0.02,0.01,0.005,0.001]:
        mu=mu_for_delta(delta)
        y,T,DPm=mother_DP(mu)
        print('\ndelta',delta,'mu',mu,'mother tr',np.trace(DPm),'y',y)
        amp=np.sqrt(delta/(4*A))
        seeds={'+Q':np.array([np.sqrt(alpha)*amp,0.]),'-Q':np.array([-np.sqrt(alpha)*amp,0.]),'+P':np.array([0.,amp/np.sqrt(alpha)]),'-P':np.array([0.,-amp/np.sqrt(alpha)])}
        for name,s in seeds.items():
            try:
                z=solve_fixed(mu,y,s)
                F,DF,z1=fmap(mu,y,z,with_dp=True)
                print(name,'seed',s,'->',z,'P->',z1,'res',np.linalg.norm(F-z),'trF',np.trace(DF),'det',np.linalg.det(DF))
            except Exception as e: print(name,'FAIL',e)
