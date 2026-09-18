import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, root_scalar


def f(t,X,mu):
    x,y,z=X
    return np.array([np.cos(z), np.sin(z), -0.5*x*x+0.5*y*y-mu])

def jac(X):
    x,y,z=X
    return np.array([[0,0,-np.sin(z)],[0,0,np.cos(z)],[-x,y,0.]],float)

def flow_with_var(mu,X0,t_end,rtol=2e-12,atol=2e-14,dense=False,events=None,max_step=np.inf):
    Y0=np.r_[X0,np.eye(3).ravel()]
    def rhs(t,Y):
        X=Y[:3]; Phi=Y[3:].reshape(3,3)
        return np.r_[f(t,X,mu),(jac(X)@Phi).ravel()]
    return solve_ivp(rhs,[0,t_end],Y0,rtol=rtol,atol=atol,dense_output=dense,events=events,max_step=max_step,method='DOP853')

def shoot_res(v,mu):
    a,Tq=v
    sol=solve_ivp(lambda t,X:f(t,X,mu),[0,Tq],[0,a,0],rtol=3e-12,atol=3e-14,method='DOP853',max_step=0.15)
    X=sol.y[:,-1]
    return [X[1], X[2]-np.pi/2]

def mother(mu,guess=(-2.5726,3.018)):
    r=root(lambda v:shoot_res(v,mu),guess,tol=1e-11)
    if not r.success:
        print('root fail',r.message,r.fun)
    return r.x

def local_return(mu,y_star,q=0,p=0,tmax=25):
    X0=np.array([0., y_star+q, np.arcsin(p)])
    # collect direction +1 x crossings; skip t=0 and select closest section point after >1
    def ev(t,Y): return Y[0]
    ev.direction=1; ev.terminal=False
    sol=flow_with_var(mu,X0,tmax,events=ev,max_step=.08)
    ts=sol.t_events[0]; Ys=sol.y_events[0]
    c=[]
    for t,Y in zip(ts,Ys):
        if t<1e-6: continue
        X=Y[:3]
        if np.cos(X[2])<=0: continue
        dist=np.hypot(X[1]-y_star,np.sin(X[2]))
        c.append((dist,t,Y))
    c.sort(key=lambda a:a[0])
    if not c: raise RuntimeError('no return')
    dist,t,Y=c[0]; Xf=Y[:3]; Phi=Y[3:].reshape(3,3)
    ff=f(t,Xf,mu)
    n=np.array([1.,0,0])
    Proj=np.eye(3)-np.outer(ff,n)/(n@ff)
    DE=np.array([[0,0],[1,0],[0,1/np.sqrt(1-p*p)]],float)
    DC=np.array([[0,1,0],[0,0,np.cos(Xf[2])]],float)
    DP=DC@Proj@Phi@DE
    out=np.array([Xf[1]-y_star,np.sin(Xf[2])])
    return out,DP,t,dist,Xf,[(a,b) for a,b,_ in c[:6]]

if __name__=='__main__':
    mu=0.66286900694
    a,Tq=mother(mu)
    print('mother',mu,a,Tq, 'res',shoot_res((a,Tq),mu))
    out,DP,t,d,Xf,cands=local_return(mu,a)
    print('return t',t,'out',out,'dist',d,'Xf',Xf)
    print('cands',cands)
    print('DP',DP,'tr',np.trace(DP),'det',np.linalg.det(DP),'DP+I',DP+np.eye(2),'norm',np.linalg.norm(DP+np.eye(2)))
