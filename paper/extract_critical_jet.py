from __future__ import annotations
import mpmath as mp
import math, json, os, sys, time
from pathlib import Path
from multiprocessing import Pool

BASE=Path('/mnt/data/carnot_resonance_paper')
sys.path.insert(0,str(BASE))
from taylor_mp import local_return

MU='0.6628690069305101280936908638603636487'
YSTAR='-2.5725893074769489665217743248705277'

# ---------- polynomial algebra, total-degree truncated ----------
def clean(p, tol=None):
    if tol is None: return {k:v for k,v in p.items() if v != 0}
    return {k:v for k,v in p.items() if abs(v)>tol}

def padd(a,b,N=7):
    c=dict(a)
    for k,v in b.items():
        if sum(k)<=N: c[k]=c.get(k,mp.mpf('0'))+v
    return clean(c)

def pscale(a,s): return clean({k:s*v for k,v in a.items()})
def pmul(a,b,N=7):
    c={}
    for (i,j),u in a.items():
        for (k,l),v in b.items():
            d=i+j+k+l
            if d<=N:
                key=(i+k,j+l); c[key]=c.get(key,mp.mpf('0'))+u*v
    return clean(c)
def ppow(a,n,N=7):
    r={(0,0):mp.mpf(1)}
    for _ in range(n): r=pmul(r,a,N)
    return r
def pder(a,var):
    c={}
    for (i,j),v in a.items():
        if var==0 and i: c[(i-1,j)]=v*i
        if var==1 and j: c[(i,j-1)]=v*j
    return c
def pint_p(a):
    return {(i,j+1):v/(j+1) for (i,j),v in a.items()}
def pcompose(poly, mpair, N=7):
    q,p=mpair; out={}
    cacheq={0:{(0,0):mp.mpf(1)}}; cachep={0:{(0,0):mp.mpf(1)}}
    for (i,j),v in poly.items():
        if i not in cacheq: cacheq[i]=ppow(q,i,N)
        if j not in cachep: cachep[j]=ppow(p,j,N)
        term=pscale(pmul(cacheq[i],cachep[j],N),v)
        out=padd(out,term,N)
    return out

def map_comp(A,B,N=7): # A o B
    return [pcompose(A[0],B,N),pcompose(A[1],B,N)]
def map_add(A,B,N=7): return [padd(A[0],B[0],N),padd(A[1],B[1],N)]
def map_scale(A,s): return [pscale(A[0],s),pscale(A[1],s)]
ID=[{(1,0):mp.mpf(1)},{(0,1):mp.mpf(1)}]

def map_disp(A,N=7): return [padd(A[0],pscale(ID[0],-1),N),padd(A[1],pscale(ID[1],-1),N)]
def map_jac_vec(A,V,N=7):
    return [padd(pmul(pder(A[i],0),V[0],N),pmul(pder(A[i],1),V[1],N),N) for i in range(2)]

def homogeneous(p,d): return {k:v for k,v in p.items() if sum(k)==d}
def map_hom(A,d): return [homogeneous(A[0],d),homogeneous(A[1],d)]

def map_inverse(A,N=7):
    # near-identity inverse by degree correction
    B=[dict(ID[0]),dict(ID[1])]
    for d in range(2,N+1):
        C=map_comp(A,B,N)
        err=map_hom(map_disp(C,N),d)
        B=[padd(B[i],pscale(err[i],-1),N) for i in range(2)]
    return B

def det_jac(A,N=7):
    a=pder(A[0],0); b=pder(A[0],1); c=pder(A[1],0); d=pder(A[1],1)
    return padd(pmul(a,d,N),pscale(pmul(b,c,N),-1),N)

def peval2(p,q,r):
    return sum(v*(q**i)*(r**j) for (i,j),v in p.items())

def norm_coeff(p, degrees=None):
    vals=[]
    for k,v in p.items():
        if degrees is None or sum(k) in degrees: vals.append(abs(v))
    return max(vals) if vals else mp.mpf('0')

# ---------- sampling ----------
def worker(args):
    qstr,pstr,dps,order,hmax=args
    import mpmath as _mp, sys as _sys
    _mp.mp.dps=dps
    _sys.path.insert(0,str(BASE))
    from taylor_mp import local_return as _lr
    mu=_mp.mpf(MU); ys=_mp.mpf(YSTAR)
    q=_mp.mpf(qstr); p=_mp.mpf(pstr)
    out=_lr(mu,ys,[q,p],order=order,hmax=_mp.mpf(hmax),radius=_mp.mpf('.2'),with_dp=False,with_action=False,tmax=_mp.mpf('14'))[0]
    return (qstr,pstr,_mp.nstr(out[0],dps),_mp.nstr(out[1],dps))

def sample_radius(radius, nnode=7,dps=70,order=34,hmax='.06',procs=6):
    # Chebyshev-Lobatto tensor grid; full return, including origin.
    nodes=[math.cos(math.pi*k/(nnode-1)) for k in range(nnode)]
    tasks=[]
    rr=mp.mpf(radius)
    for x in nodes:
        for y in nodes:
            q=rr*mp.mpf(str(x)); p=rr*mp.mpf(str(y))
            tasks.append((mp.nstr(q,dps),mp.nstr(p,dps),dps,order,hmax))
    t=time.time()
    with Pool(processes=procs) as pool:
        vals=list(pool.imap_unordered(worker,tasks,chunksize=1))
    vals.sort(key=lambda z:(mp.mpf(z[0]),mp.mpf(z[1])))
    print(f'sampled radius={radius} n={len(vals)} in {time.time()-t:.1f}s',flush=True)
    return vals

def fit_jet(vals,radius,maxdeg=7,mindeg=2,dps=70):
    mp.mp.dps=dps
    mons=[(i,d-i) for d in range(mindeg,maxdeg+1) for i in range(d+1)]
    m=len(vals); n=len(mons)
    A=mp.matrix(m,n); bq=mp.matrix(m,1); bp=mp.matrix(m,1)
    rr=mp.mpf(radius)
    for row,(qs,ps,oqs,ops) in enumerate(vals):
        q=mp.mpf(qs); p=mp.mpf(ps); oq=mp.mpf(oqs); op=mp.mpf(ops)
        x=q/rr; y=p/rr
        for col,(i,j) in enumerate(mons): A[row,col]=x**i*y**j
        bq[row]=oq+q; bp[row]=op+p
    cq,resq=mp.qr_solve(A,bq); cp,resp=mp.qr_solve(A,bp)
    P=[{(1,0):mp.mpf(-1)},{(0,1):mp.mpf(-1)}]
    for col,(i,j) in enumerate(mons):
        d=i+j
        P[0][(i,j)]=cq[col]/rr**d
        P[1][(i,j)]=cp[col]/rr**d
    return P, {'resq':resq,'resp':resp,'mons':mons}

# ---------- normal-form construction ----------
def build_normal_form(P,N=5):
    # F=P^2 = Id + f3+f4+f5 (quadratic cancels formally).
    F=map_comp(P,P,N)
    V=map_disp(F,N)
    V3=map_hom(V,3); V4=map_hom(V,4); V5=map_hom(V,5)
    # log F through degree 5 (if V starts degree 3): X3=V3, X4=V4, X5=V5-1/2 DV3 V3
    DV3V3=map_jac_vec(V3,V3,N)
    X=[padd(padd(V3[0],V4[0],N),padd(V5[0],pscale(DV3V3[0],-mp.mpf('.5')),N),N),
       padd(padd(V3[1],V4[1],N),padd(V5[1],pscale(DV3V3[1],-mp.mpf('.5')),N),N)]
    # M=exp(X/2) through deg5: Id + X/2 + 1/8 D X3 X3
    M=[dict(ID[0]),dict(ID[1])]
    for i in range(2):
        M[i]=padd(M[i],pscale(X[i],mp.mpf('.5')),N)
        M[i]=padd(M[i],pscale(DV3V3[i],mp.mpf('.125')),N)
    Minv=map_inverse(M,N)
    R=map_comp(P,Minv,N) # involution, near -Id
    # u=(w-Rw)/2, derivative Id
    U=[pscale(padd(ID[i],pscale(R[i],-1),N),mp.mpf('.5')) for i in range(2)]
    Uinv=map_inverse(U,N)
    # omega in u coords: f(u)=det D(U^{-1})(u)
    dens=det_jac(Uinv,N-1)
    # Darboux odd correction C(u)=(u1, integral_0^u2 dens(u1,s) ds)
    C=[dict(ID[0]),pint_p(dens)]
    # Trim C to N; should be odd if dens even.
    C=[{k:v for k,v in ci.items() if sum(k)<=N} for ci in C]
    Nmap=map_comp(C,U,N)  # original w -> canonical odd coordinates W
    Ninv=map_inverse(Nmap,N)
    Ptilde=map_comp(map_comp(Nmap,P,N),Ninv,N)
    Ftilde=map_comp(Ptilde,Ptilde,N)
    Vf=map_disp(Ftilde,N); f3=map_hom(Vf,3); f4=map_hom(Vf,4); f5=map_hom(Vf,5)
    Df3f3=map_jac_vec(f3,f3,N)
    # H generator for Mtilde = sqrt(Ftilde): vector field Y=1/2 log Ftilde
    # Y3=.5f3; Y4=.5f4; Y5=.5(f5-.5 Df3 f3)
    Y3=map_scale(f3,mp.mpf('.5'))
    Y4=map_scale(f4,mp.mpf('.5'))
    Y5=[pscale(padd(f5[i],pscale(Df3f3[i],-mp.mpf('.5')),N),mp.mpf('.5')) for i in range(2)]
    return {'F':F,'XlogF':X,'M':M,'R':R,'U':U,'dens':dens,'C':C,'Nmap':Nmap,'Ninv':Ninv,'Ptilde':Ptilde,'Ftilde':Ftilde,'Y3':Y3,'Y4':Y4,'Y5':Y5}

def integrate_hamiltonian_from_vector(Y, Hdeg):
    # Y=(H_p,-H_q). Integrate H_p=Yq wrt p; correct q-only term using -H_q=Yp.
    H=pint_p(Y[0])
    Hq=pder(H,0)
    resid=padd(Y[1],Hq,7)  # Yp + H_q should be 0; any q-only part correction needed? desired H_q=-Yp
    # correction C'(q) = -Yp - Hq = -resid; integrate terms with p^0
    corr={}
    for (i,j),v in resid.items():
        if j==0:
            corr[(i+1,0)]=-v/(i+1)
    H=padd(H,corr,7)
    return {k:v for k,v in H.items() if sum(k)==Hdeg}

def linear_scale_H(H,alpha,N=7):
    # old section canonical coords w=(q,p); q=sqrt(alpha) Q, p=P/sqrt(alpha)
    out={}; sa=mp.sqrt(alpha)
    for (i,j),v in H.items(): out[(i,j)]=v*(sa**i)*(sa**(-j))
    return out

def dump_poly(p):
    return {f'{i},{j}':mp.nstr(v,50) for (i,j),v in sorted(p.items(),key=lambda kv:(sum(kv[0]),kv[0]))}

def run(radius='0.0006', dps=70, nnode=7, procs=6):
    mp.mp.dps=dps
    cache=BASE/'data'/f'jet_samples_r{radius.replace(".","p")}.json'
    if cache.exists(): vals=json.loads(cache.read_text())
    else:
        vals=sample_radius(radius,nnode,dps,34,'.06',procs)
        cache.write_text(json.dumps(vals))
    P,fit=fit_jet(vals,radius,7,2,dps)
    nf=build_normal_form(P,5)
    H4=integrate_hamiltonian_from_vector(nf['Y3'],4)
    H5=integrate_hamiltonian_from_vector(nf['Y4'],5)
    H6=integrate_hamiltonian_from_vector(nf['Y5'],6)
    # diagnostics in raw section coordinates
    ev2=map_hom(nf['Ptilde'],2); ev4=map_hom(nf['Ptilde'],4); od3=map_hom(nf['Ptilde'],3); od5=map_hom(nf['Ptilde'],5)
    evenPtilde=max(norm_coeff(ev2[0]),norm_coeff(ev2[1]),norm_coeff(ev4[0]),norm_coeff(ev4[1]))
    oddscale=max(norm_coeff(od3[0]),norm_coeff(od3[1]),norm_coeff(od5[0]),norm_coeff(od5[1]))
    dens_odd=max([abs(v) for k,v in nf['dens'].items() if sum(k)%2==1] or [mp.mpf('0')])
    divY3=norm_coeff(padd(pder(nf['Y3'][0],0),pder(nf['Y3'][1],1),7))
    divY5=norm_coeff(padd(pder(nf['Y5'][0],0),pder(nf['Y5'][1],1),7))
    out={
      'radius':radius,'dps':dps,'fit_resq':mp.nstr(fit['resq'],30),'fit_resp':mp.nstr(fit['resp'],30),
      'P': [dump_poly(P[0]),dump_poly(P[1])],
      'F2norm':mp.nstr(max(norm_coeff(map_hom(map_disp(nf['F'],5),2)[0]),norm_coeff(map_hom(map_disp(nf['F'],5),2)[1])),30),
      'R_involution_error': None,
      'dens_odd':mp.nstr(dens_odd,30),'oddification_even_ratio':mp.nstr(evenPtilde/oddscale,30),
      'divY3':mp.nstr(divY3,30),'divY5':mp.nstr(divY5,30),
      'H4_raw':dump_poly(H4),'H5_raw':dump_poly(H5),'H6_raw':dump_poly(H6),
      'Nmap': [dump_poly(nf['Nmap'][0]),dump_poly(nf['Nmap'][1])],
      'Ninv': [dump_poly(nf['Ninv'][0]),dump_poly(nf['Ninv'][1])],
      'Ptilde':[dump_poly(nf['Ptilde'][0]),dump_poly(nf['Ptilde'][1])]
    }
    RR=map_comp(nf['R'],nf['R'],5); rrdisp=map_disp(RR,5)
    out['R_involution_error']=mp.nstr(max(norm_coeff(rrdisp[0]),norm_coeff(rrdisp[1])),30)
    outfile=BASE/'data'/f'critical_jet_r{radius.replace(".","p")}.json'
    outfile.write_text(json.dumps(out,indent=2))
    print('radius',radius,'fit',fit['resq'],fit['resp'])
    print('F2',out['F2norm'],'R2',out['R_involution_error'],'densodd',out['dens_odd'],'evenratio',out['oddification_even_ratio'])
    print('div',out['divY3'],out['divY5'])
    print('H4raw',out['H4_raw'])
    print('H5max',mp.nstr(norm_coeff(H5),20))
    print('H6raw',out['H6_raw'])
    return out

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument('--radius',default='0.0006'); ap.add_argument('--dps',type=int,default=70); ap.add_argument('--nodes',type=int,default=7); ap.add_argument('--procs',type=int,default=6)
    args=ap.parse_args(); run(args.radius,args.dps,args.nodes,args.procs)

# ---------- nonresonant Lie-transform oddification ----------
def vector_flow_map(Y,N=5):
    """Time-one flow map exp(Y) for polynomial vector field Y, truncated."""
    out=[dict(ID[0]),dict(ID[1])]
    V=[dict(Y[0]),dict(Y[1])]
    fact=mp.mpf(1)
    for n in range(1,8):
        fact *= n
        out=[padd(out[i],pscale(V[i],1/fact),N) for i in range(2)]
        Vnext=map_jac_vec(V,Y,N)
        if max(norm_coeff(Vnext[0]),norm_coeff(Vnext[1]))==0: break
        V=Vnext
    return out

def oddify_nonresonant(P,N=5):
    # First eliminate degree-2 term with canonical homogeneous Y2=-P2/2.
    P2=map_hom(P,2)
    Y2=map_scale(P2,mp.mpf('-.5'))
    # Hamiltonian diagnostic
    div2=padd(pder(Y2[0],0),pder(Y2[1],1),N)
    H3=integrate_hamiltonian_from_vector(Y2,3)
    H2map=vector_flow_map(Y2,N); H2inv=vector_flow_map(map_scale(Y2,-1),N)
    P1=map_comp(map_comp(H2map,P,N),H2inv,N)
    # Eliminate resulting degree-4 term with Y4=-P1_4/2.
    P14=map_hom(P1,4)
    Y4=map_scale(P14,mp.mpf('-.5'))
    div4=padd(pder(Y4[0],0),pder(Y4[1],1),N)
    H5gen=integrate_hamiltonian_from_vector(Y4,5)
    H4map=vector_flow_map(Y4,N); H4inv=vector_flow_map(map_scale(Y4,-1),N)
    Pt=map_comp(map_comp(H4map,P1,N),H4inv,N)
    # overall coordinate map old->new is H4map o H2map (since each conjugation h P h^-1)
    Nmap=map_comp(H4map,H2map,N); Ninv=map_inverse(Nmap,N)
    return {'Ptilde':Pt,'Nmap':Nmap,'Ninv':Ninv,'Y2':Y2,'Y4':Y4,'G3':H3,'G5':H5gen,'div2':div2,'div4':div4}

def generator_from_odd_P(Pt,N=5):
    # M=-Pt = Id + M3+M5; log directly.
    M=[pscale(Pt[0],-1),pscale(Pt[1],-1)]
    V=map_disp(M,N); M3=map_hom(V,3); M4=map_hom(V,4); M5=map_hom(V,5)
    DM3M3=map_jac_vec(M3,M3,N)
    Y3=M3
    Y5=[padd(M5[i],pscale(DM3M3[i],mp.mpf('-.5')),N) for i in range(2)]
    H4=integrate_hamiltonian_from_vector(Y3,4); H6=integrate_hamiltonian_from_vector(Y5,6)
    return {'M':M,'M3':M3,'M4':M4,'M5':M5,'Y3':Y3,'Y5':Y5,'H4':H4,'H6':H6}
