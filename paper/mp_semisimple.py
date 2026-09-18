import mpmath as mp, time

def state_quarter(mu,a,T,dps=50,tol_exp=38,degree=36):
    def ff(t,Y):
        x,y,z=Y
        return (mp.cos(z),mp.sin(z),-x*x/2+y*y/2-mu)
    sol=mp.odefun(ff, mp.mpf('0'), (mp.mpf('0'),a,mp.mpf('0')), tol=mp.mpf(10)**(-tol_exp), degree=degree)
    return sol(T)

def mother(mu,a0,T0,dps=50,tol_exp=38,degree=36):
    mp.mp.dps=dps
    def F(a,T):
        Y=state_quarter(mu,a,T,dps,tol_exp,degree)
        return Y[1], Y[2]-mp.pi/2
    a,T=mp.findroot(F,(mp.mpf(a0),mp.mpf(T0)),tol=mp.mpf(10)**(-(dps-10)),maxsteps=12,solver='mdnewton')
    return a,T

def full_DP(mu,a,Tq,dps=50,tol_exp=36,degree=32):
    mp.mp.dps=dps
    def rhs(t,Y):
        x,y,z=Y[0],Y[1],Y[2]
        J=((mp.mpf(0),mp.mpf(0),-mp.sin(z)),(mp.mpf(0),mp.mpf(0),mp.cos(z)),(-x,y,mp.mpf(0)))
        ph=[Y[3+i] for i in range(9)]
        out=[mp.cos(z),mp.sin(z),-x*x/2+y*y/2-mu]
        for i in range(3):
            for j in range(3):
                out.append(sum(J[i][k]*ph[3*k+j] for k in range(3)))
        return tuple(out)
    Y0=[mp.mpf('0'),a,mp.mpf('0')]+[mp.mpf(1 if i==j else 0) for i in range(3) for j in range(3)]
    sol=mp.odefun(rhs,mp.mpf('0'),tuple(Y0),tol=mp.mpf(10)**(-tol_exp),degree=degree)
    Y=sol(4*Tq)
    x,y,z=Y[0],Y[1],Y[2]
    Phi=mp.matrix(3,3)
    for i in range(3):
        for j in range(3): Phi[i,j]=Y[3+3*i+j]
    ff=mp.matrix([mp.cos(z),mp.sin(z),-x*x/2+y*y/2-mu])
    n=mp.matrix([[1,0,0]])
    Proj=mp.eye(3)-(ff*n)/(ff[0])
    DE=mp.matrix([[0,0],[1,0],[0,1]])
    DC=mp.matrix([[0,1,0],[0,0,mp.cos(z)]])
    DP=DC*Proj*Phi*DE
    return Y[:3],DP

if __name__=='__main__':
    mp.mp.dps=45
    mu=mp.mpf('0.66286900693055')
    t=time.time(); a,T=mother(mu,'-2.572589307477','3.01795545336',dps=45,tol_exp=34,degree=32); print('mother time',time.time()-t); print(mp.nstr(a,30),mp.nstr(T,30))
    t=time.time(); X,DP=full_DP(mu,a,T,dps=45,tol_exp=32,degree=28); print('DP time',time.time()-t); print([mp.nstr(v,30) for v in X]); print([[mp.nstr(DP[i,j],30) for j in range(2)] for i in range(2)]); print('tr+2',mp.nstr(DP[0,0]+DP[1,1]+2,30)); print('b',mp.nstr(DP[0,1],30),'c',mp.nstr(DP[1,0],30)); print('det-1',mp.nstr(mp.det(DP)-1,30))
