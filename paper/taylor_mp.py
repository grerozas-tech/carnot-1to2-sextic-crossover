import mpmath as mp


def taylor_coeffs(mu, state, Phi=None, order=24):
    # coefficients t^n about current point for state and optionally 3x3 variational matrix
    x=[mp.mpf('0')]*(order+1); y=x.copy(); z=x.copy(); S=x.copy(); C=x.copy()
    x[0],y[0],z[0]=map(mp.mpf,state)
    S[0]=mp.sin(z[0]); C[0]=mp.cos(z[0])
    ph=None
    if Phi is not None:
        ph=[[[mp.mpf('0')]*(order+1) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3): ph[i][j][0]=mp.mpf(Phi[i][j])
    for n in range(order):
        # state coefficient n+1
        x[n+1]=C[n]/(n+1)
        y[n+1]=S[n]/(n+1)
        xx=sum(x[k]*x[n-k] for k in range(n+1))
        yy=sum(y[k]*y[n-k] for k in range(n+1))
        rhsz=(-xx+yy)/2
        if n==0: rhsz -= mu
        z[n+1]=rhsz/(n+1)
        # sin/cos composition coefficients n+1 using S'=C z', C'=-S z'
        convS=mp.mpf('0'); convC=mp.mpf('0')
        for k in range(n+1):
            zd=(n-k+1)*z[n-k+1]
            convS += C[k]*zd
            convC -= S[k]*zd
        S[n+1]=convS/(n+1)
        C[n+1]=convC/(n+1)
        if ph is not None:
            # J coefficients at m: rows: [0,0,-S]; [0,0,C]; [-x,y,0]
            for i in range(3):
                for j in range(3):
                    cv=mp.mpf('0')
                    for k in range(n+1):
                        m=k; r=n-k
                        if i==0:
                            cv += (-S[m])*ph[2][j][r]
                        elif i==1:
                            cv += C[m]*ph[2][j][r]
                        else:
                            cv += (-x[m])*ph[0][j][r] + y[m]*ph[1][j][r]
                    ph[i][j][n+1]=cv/(n+1)
    return x,y,z,S,C,ph

def peval(co,h):
    v=co[-1]
    for a in reversed(co[:-1]): v=v*h+a
    return v

def step(mu,state,h,Phi=None,order=24,return_coeffs=False):
    co=taylor_coeffs(mu,state,Phi,order)
    x,y,z,S,C,ph=co
    out=[peval(x,h),peval(y,h),peval(z,h)]
    Pout=None
    if Phi is not None:
        Pout=[[peval(ph[i][j],h) for j in range(3)] for i in range(3)]
    if return_coeffs: return out,Pout,co
    return out,Pout

def integrate(mu,state,T,Phi=None,order=24,hmax=mp.mpf('0.15')):
    T=mp.mpf(T); n=max(1,int(mp.ceil(abs(T)/hmax))); h=T/n
    st=list(map(mp.mpf,state))
    P=[[mp.mpf(Phi[i][j]) for j in range(3)] for i in range(3)] if Phi is not None else None
    for _ in range(n): st,P=step(mu,st,h,P,order)
    return st,P

def mother(mu,a0,T0,order=24,hmax=mp.mpf('0.15'),tol=None,maxit=10):
    a=mp.mpf(a0); T=mp.mpf(T0)
    if tol is None: tol=mp.mpf(10)**(-(mp.mp.dps-12))
    I=[[mp.mpf(1 if i==j else 0) for j in range(3)] for i in range(3)]
    for it in range(maxit):
        st,P=integrate(mu,[0,a,0],T,I,order,hmax)
        r1=st[1]; r2=st[2]-mp.pi/2
        if max(abs(r1),abs(r2))<tol: return a,T,it+1
        # Jacobian columns: derivative wrt a = Phi[:,1], wrt T = f(final)
        x,y,z=st
        fy=mp.sin(z); fz=-x*x/2+y*y/2-mu
        J=mp.matrix([[P[1][1],fy],[P[2][1],fz]])
        r=mp.matrix([r1,r2]); d=mp.lu_solve(J,-r)
        a += d[0]; T += d[1]
    return a,T,maxit

def full_DP(mu,a,Tq,order=24,hmax=mp.mpf('0.15')):
    I=[[mp.mpf(1 if i==j else 0) for j in range(3)] for i in range(3)]
    st,Phi=integrate(mu,[0,a,0],4*Tq,I,order,hmax)
    x,y,z=st
    ff=mp.matrix([mp.cos(z),mp.sin(z),-x*x/2+y*y/2-mu])
    # projection I - f n^T/(n^T f), n=(1,0,0)
    Proj=mp.eye(3)
    for i in range(3): Proj[i,0]-=ff[i]/ff[0]
    PM=mp.matrix(Phi)
    DE=mp.matrix([[0,0],[1,0],[0,1]])
    DC=mp.matrix([[0,1,0],[0,0,mp.cos(z)]])
    DP=DC*Proj*PM*DE
    return st,DP

if __name__=='__main__':
    mp.mp.dps=50
    mu=mp.mpf('0.66286900693055')
    for h in ['0.20','0.15','0.10']:
        a,T,it=mother(mu,'-2.57258930748','3.01795545336',order=26,hmax=mp.mpf(h),maxit=8)
        st,DP=full_DP(mu,a,T,order=26,hmax=mp.mpf(h))
        print('h',h,'a',mp.nstr(a,30),'T',mp.nstr(T,30),'it',it)
        print('state',*[mp.nstr(v,20) for v in st])
        print('DP',[[mp.nstr(DP[i,j],25) for j in range(2)] for i in range(2)])
        print('tr+2',mp.nstr(DP[0,0]+DP[1,1]+2,20),'det-1',mp.nstr(mp.det(DP)-1,20))

# ---- local-return utilities with optional action ----
def _conv_at(a,b,n):
    return sum(a[k]*b[n-k] for k in range(n+1))

def _conv3_at(a,b,c,n):
    s=mp.mpf('0')
    for i in range(n+1):
        for j in range(n-i+1):
            s += a[i]*b[j]*c[n-i-j]
    return s

def taylor_coeffs_ext(mu,state,Phi=None,action0=None,order=24):
    x=[mp.mpf('0')]*(order+1); y=x.copy(); z=x.copy(); S=x.copy(); C=x.copy(); Ac=x.copy() if action0 is not None else None
    x[0],y[0],z[0]=map(mp.mpf,state); S[0]=mp.sin(z[0]); C[0]=mp.cos(z[0])
    if Ac is not None: Ac[0]=mp.mpf(action0)
    ph=None
    if Phi is not None:
        ph=[[[mp.mpf('0')]*(order+1) for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3): ph[i][j][0]=mp.mpf(Phi[i][j])
    for n in range(order):
        x[n+1]=C[n]/(n+1); y[n+1]=S[n]/(n+1)
        xx=_conv_at(x,x,n); yy=_conv_at(y,y,n)
        rz=(-xx+yy)/2 - (mu if n==0 else 0)
        z[n+1]=rz/(n+1)
        cs=mp.mpf('0'); cc=mp.mpf('0')
        for k in range(n+1):
            zd=(n-k+1)*z[n-k+1]; cs += C[k]*zd; cc -= S[k]*zd
        S[n+1]=cs/(n+1); C[n+1]=cc/(n+1)
        if ph is not None:
            for i in range(3):
                for j in range(3):
                    cv=mp.mpf('0')
                    for k in range(n+1):
                        r=n-k
                        if i==0: cv += (-S[k])*ph[2][j][r]
                        elif i==1: cv += C[k]*ph[2][j][r]
                        else: cv += (-x[k])*ph[0][j][r]+y[k]*ph[1][j][r]
                    ph[i][j][n+1]=cv/(n+1)
        if Ac is not None:
            # af = -1 + x*S*( y^2/2 - mu - x^2/6 )
            y2=[_conv_at(y,y,m) for m in range(n+1)]
            x2=[_conv_at(x,x,m) for m in range(n+1)]
            B=[y2[m]/2 - x2[m]/6 - (mu if m==0 else 0) for m in range(n+1)]
            af=(-1 if n==0 else 0) + _conv3_at(x,S,B,n)
            Ac[n+1]=af/(n+1)
    return x,y,z,S,C,ph,Ac

def step_ext(mu,state,h,Phi=None,action0=None,order=24,return_coeffs=False):
    co=taylor_coeffs_ext(mu,state,Phi,action0,order)
    x,y,z,S,C,ph,Ac=co
    st=[peval(x,h),peval(y,h),peval(z,h)]
    Pout=None if Phi is None else [[peval(ph[i][j],h) for j in range(3)] for i in range(3)]
    Aout=None if Ac is None else peval(Ac,h)
    return (st,Pout,Aout,co) if return_coeffs else (st,Pout,Aout)

def _root_poly_bisect(co,h,its=180):
    lo=mp.mpf('0'); hi=mp.mpf(h); flo=peval(co,lo); fhi=peval(co,hi)
    if flo==0: lo=mp.mpf('1e-40')*h; flo=peval(co,lo)
    if flo*fhi>0: return None
    for _ in range(its):
        mid=(lo+hi)/2; fm=peval(co,mid)
        if flo*fm<=0: hi=mid; fhi=fm
        else: lo=mid; flo=fm
    return (lo+hi)/2

def local_return(mu,ystar,qp,order=26,hmax=mp.mpf('.08'),radius=mp.mpf('1.0'),with_dp=True,with_action=True,tmax=mp.mpf('15')):
    q,p=map(mp.mpf,qp); st=[mp.mpf('0'),ystar+q,mp.asin(p)]
    Phi=[[mp.mpf(1 if i==j else 0) for j in range(3)] for i in range(3)] if with_dp else None
    act=mp.mpf('0') if with_action else None
    t=mp.mpf('0'); first=True
    while t < tmax:
        h=min(hmax,tmax-t)
        st_end,P_end,A_end,co=step_ext(mu,st,h,Phi,act,order,True)
        xco,yco,zco,Sco,Cco,phco,Aco=co
        # detect negative-to-positive crossing inside step, excluding initial t=0
        x0=st[0]; x1=st_end[0]
        if (x0 < 0 and x1 >= 0) or (not first and x0==0 and x1>0):
            hr=_root_poly_bisect(xco,h)
            if hr is not None and hr>mp.mpf('1e-25'):
                sr=[peval(xco,hr),peval(yco,hr),peval(zco,hr)]
                qr=sr[1]-ystar; pr=mp.sin(sr[2]); dist=mp.sqrt(qr*qr+pr*pr)
                if mp.cos(sr[2])>0 and dist<radius and t+hr>mp.mpf('1e-8'):
                    P_r=None
                    if with_dp:
                        P_r=[[peval(phco[i][j],hr) for j in range(3)] for i in range(3)]
                        ff=mp.matrix([mp.cos(sr[2]),mp.sin(sr[2]),-sr[0]*sr[0]/2+sr[1]*sr[1]/2-mu])
                        Proj=mp.eye(3)
                        for i in range(3): Proj[i,0]-=ff[i]/ff[0]
                        DE=mp.matrix([[0,0],[1,0],[0,1/mp.sqrt(1-p*p)]])
                        DC=mp.matrix([[0,1,0],[0,0,mp.cos(sr[2])]])
                        DP=DC*Proj*mp.matrix(P_r)*DE
                    else: DP=None
                    ar=peval(Aco,hr) if with_action else None
                    return [qr,pr],DP,ar,t+hr,sr
        st,Phi,act=st_end,P_end,A_end; t+=h; first=False
    raise RuntimeError('local return not found')

def fmap_local(mu,ystar,qp,order=26,hmax=mp.mpf('.08'),with_dp=True,with_action=True):
    z1,D1,A1,t1,s1=local_return(mu,ystar,qp,order,hmax,with_dp=with_dp,with_action=with_action)
    z2,D2,A2,t2,s2=local_return(mu,ystar,z1,order,hmax,with_dp=with_dp,with_action=with_action)
    DF=D2*D1 if with_dp else None
    AF=A1+A2 if with_action else None
    return z2,DF,AF,z1,(t1,t2)

def fixed_newton(mu,ystar,seed,order=26,hmax=mp.mpf('.08'),tol=None,maxit=8):
    z=mp.matrix([mp.mpf(seed[0]),mp.mpf(seed[1])])
    if tol is None: tol=mp.mpf(10)**(-(mp.mp.dps-15))
    for it in range(maxit):
        F,DF,AF,z1,ts=fmap_local(mu,ystar,[z[0],z[1]],order,hmax,True,False)
        r=mp.matrix([F[0]-z[0],F[1]-z[1]])
        if max(abs(r[0]),abs(r[1]))<tol: return [z[0],z[1]],it+1
        J=DF-mp.eye(2); dz=mp.lu_solve(J,-r); z += dz
    return [z[0],z[1]],maxit
