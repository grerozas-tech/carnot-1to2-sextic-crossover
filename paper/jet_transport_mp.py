from __future__ import annotations
import mpmath as mp
from pathlib import Path
import sys,json,time
BASE=Path('/mnt/data/carnot_resonance_paper')
sys.path.insert(0,str(BASE))
import extract_critical_jet as alg

MU=mp.mpf('0.6628690069305101280936908638603636487')
YSTAR=mp.mpf('-2.5725893074769489665217743248705277')
TQ=mp.mpf('3.0179554533608596100474613978029132')

N=5
ZERO={}
ONE={(0,0):mp.mpf(1)}
Q={(1,0):mp.mpf(1)}
P={(0,1):mp.mpf(1)}

def pc(a): return dict(a)
def const(c): return {(0,0):mp.mpf(c)} if c else {}
def add(a,b): return alg.padd(a,b,N)
def scale(a,s): return alg.pscale(a,s)
def mul(a,b): return alg.pmul(a,b,N)
def power(a,n): return alg.ppow(a,n,N)
def hpart(a,d): return alg.homogeneous(a,d)
def peval_time(arr,h):
    out={}
    hp=mp.mpf(1)
    for a in arr:
        out=add(out,scale(a,hp)); hp*=h
    return out

def pnorm(a): return alg.norm_coeff(a)

def sincos_of_poly(z):
    c=z.get((0,0),mp.mpf(0)); dz=dict(z); dz[(0,0)]=dz.get((0,0),0)-c; dz=alg.clean(dz)
    sind={}; cosd=dict(ONE)
    # sin dz, cos dz through degree N
    for k in range(0,N+1):
        n=2*k+1
        if n<=N: sind=add(sind,scale(power(dz,n), (-1)**k/mp.factorial(n)))
        n=2*k
        if n>=2 and n<=N: cosd=add(cosd,scale(power(dz,n), (-1)**k/mp.factorial(n)))
    S=add(scale(cosd,mp.sin(c)),scale(sind,mp.cos(c)))
    C=add(scale(cosd,mp.cos(c)),scale(sind,-mp.sin(c)))
    return S,C

def taylor_coeffs_jet(mu,state,S0=None,C0=None,order=24):
    x=[{} for _ in range(order+1)]; y=[{} for _ in range(order+1)]; z=[{} for _ in range(order+1)]
    S=[{} for _ in range(order+1)]; C=[{} for _ in range(order+1)]
    x[0],y[0],z[0]=map(pc,state)
    if S0 is None or C0 is None: S[0],C[0]=sincos_of_poly(z[0])
    else: S[0],C[0]=pc(S0),pc(C0)
    for n in range(order):
        x[n+1]=scale(C[n],mp.mpf(1)/(n+1))
        y[n+1]=scale(S[n],mp.mpf(1)/(n+1))
        xx={}; yy={}
        for k in range(n+1):
            xx=add(xx,mul(x[k],x[n-k])); yy=add(yy,mul(y[k],y[n-k]))
        rz=scale(add(scale(xx,-1),yy),mp.mpf('.5'))
        if n==0: rz=add(rz,const(-mu))
        z[n+1]=scale(rz,mp.mpf(1)/(n+1))
        cs={}; cc={}
        for k in range(n+1):
            zd=scale(z[n-k+1],n-k+1)
            cs=add(cs,mul(C[k],zd)); cc=add(cc,scale(mul(S[k],zd),-1))
        S[n+1]=scale(cs,mp.mpf(1)/(n+1)); C[n+1]=scale(cc,mp.mpf(1)/(n+1))
    return x,y,z,S,C

def step(mu,state,S,C,h,order=24):
    co=taylor_coeffs_jet(mu,state,S,C,order)
    x,y,z,ss,cc=co
    return [peval_time(x,h),peval_time(y,h),peval_time(z,h)],peval_time(ss,h),peval_time(cc,h)

def integrate(mu,state,S,C,T,order=24,hmax=mp.mpf('.08')):
    n=max(1,int(mp.ceil(abs(T)/hmax))); h=T/n
    st=list(map(pc,state)); ss=pc(S); cc=pc(C)
    for k in range(n):
        st,ss,cc=step(mu,st,ss,cc,h,order)
        if (k+1)%50==0: print(' step',k+1,'/',n,flush=True)
    return st,ss,cc

def compose_time_series(arr,tau):
    out={}; tp=dict(ONE)
    for a in arr:
        out=add(out,mul(a,tp)); tp=mul(tp,tau)
    return out

def return_jet(dps=70,order=24,hmax='.08'):
    mp.mp.dps=dps
    mu=mp.mpf(str(MU)); ys=mp.mpf(str(YSTAR)); T=4*mp.mpf(str(TQ))
    # initial section embedding: y=ys+q, z=asin p; sin z=p, cos z=sqrt(1-p^2)
    z0=add(P,add(scale(power(P,3),mp.mpf(1)/6),scale(power(P,5),mp.mpf(3)/40)))
    y0=add(const(ys),Q); x0={}
    S0=dict(P)
    C0=add(ONE,add(scale(power(P,2),mp.mpf(-1)/2),scale(power(P,4),mp.mpf(-1)/8)))
    t=time.time(); st,Sf,Cf=integrate(mu,[x0,y0,z0],S0,C0,T,order,mp.mpf(hmax)); print('integrated',time.time()-t,'s')
    # center constant endpoint exactly on mother section; z returns to 2*pi for this branch
    st[0][(0,0)]=mp.mpf(0); st[1][(0,0)]=ys; st[2][(0,0)]=2*mp.pi
    Sf[(0,0)]=mp.mpf(0); Cf[(0,0)]=mp.mpf(1)
    # local time Taylor around reference return, only order N needed
    tx,ty,tz,tS,tC=taylor_coeffs_jet(mu,st,Sf,Cf,order=N)
    tau={}
    for d in range(1,N+1):
        xr=compose_time_series(tx,tau)
        rd=hpart(xr,d)
        tau=add(tau,scale(rd,-1))  # x_t(0)=1
    xev=compose_time_series(tx,tau); yev=compose_time_series(ty,tau); zev=compose_time_series(tz,tau)
    qout=add(yev,const(-ys))
    pout,_=sincos_of_poly(zev)
    qout[(0,0)]=mp.mpf(0); pout[(0,0)]=mp.mpf(0)
    Pmap=[alg.clean(qout),alg.clean(pout)]
    return Pmap,tau,{'fixed_state':st,'xevent':xev}

def summarize(dps=70,order=24,hmax='.08'):
    mp.mp.dps=dps
    Pmap,tau,extra=return_jet(dps,order,hmax)
    # force tiny linear deviations? report first
    print('linear qout', {k:mp.nstr(v,30) for k,v in Pmap[0].items() if sum(k)<=1})
    print('linear pout', {k:mp.nstr(v,30) for k,v in Pmap[1].items() if sum(k)<=1})
    # Use actual linear terms; normal-form algebra assumes exactly -I, set to exact after diagnostics
    linerr=max(abs(Pmap[0].get((1,0),0)+1),abs(Pmap[0].get((0,1),0)),abs(Pmap[1].get((1,0),0)),abs(Pmap[1].get((0,1),0)+1))
    Pmap[0][(1,0)]=mp.mpf(-1); Pmap[0].pop((0,1),None); Pmap[1][(0,1)]=mp.mpf(-1); Pmap[1].pop((1,0),None)
    det=alg.det_jac(Pmap,5); det[(0,0)]=det.get((0,0),0)-1
    sym={str(d):mp.nstr(alg.norm_coeff(alg.homogeneous(det,d)),30) for d in range(1,5)}
    od=alg.oddify_nonresonant(Pmap,5)
    gen=alg.generator_from_odd_P(od['Ptilde'],5)
    H4=gen['H4']; H6=gen['H6']; H5={}
    ev2=alg.map_hom(od['Ptilde'],2); ev4=alg.map_hom(od['Ptilde'],4); od3=alg.map_hom(od['Ptilde'],3); od5=alg.map_hom(od['Ptilde'],5)
    even=max(alg.norm_coeff(ev2[0]),alg.norm_coeff(ev2[1]),alg.norm_coeff(ev4[0]),alg.norm_coeff(ev4[1])); odd=max(alg.norm_coeff(od3[0]),alg.norm_coeff(od3[1]),alg.norm_coeff(od5[0]),alg.norm_coeff(od5[1]))
    div3=alg.norm_coeff(alg.padd(alg.pder(gen['Y3'][0],0),alg.pder(gen['Y3'][1],1),7)); div5=alg.norm_coeff(alg.padd(alg.pder(gen['Y5'][0],0),alg.pder(gen['Y5'][1],1),7))
    alpha=mp.mpf('0.0946271414514154')
    H4s=alg.linear_scale_H(H4,alpha); H6s=alg.linear_scale_H(H6,alpha)
    alpha_bal=(H4.get((0,4),mp.mpf(0))/H4.get((4,0),mp.mpf(1)))**(mp.mpf(1)/4)
    out={'dps':dps,'order':order,'hmax':hmax,'linear_error_before_forcing':mp.nstr(linerr,30),'symplectic_residual_by_degree':sym,'P':[alg.dump_poly(Pmap[0]),alg.dump_poly(Pmap[1])],'tau':alg.dump_poly(tau),'H4_raw':alg.dump_poly(H4),'H6_raw':alg.dump_poly(H6),'H4_scaled':alg.dump_poly(H4s),'H6_scaled':alg.dump_poly(H6s),'alpha_used':mp.nstr(alpha,30),'alpha_bal':mp.nstr(alpha_bal,30),'Nmap_nonresonant':[alg.dump_poly(od['Nmap'][0]),alg.dump_poly(od['Nmap'][1])],'Ninv_nonresonant':[alg.dump_poly(od['Ninv'][0]),alg.dump_poly(od['Ninv'][1])],'oddification_even_ratio':mp.nstr(even/odd,30),'divY3':mp.nstr(div3,30),'divY5':mp.nstr(div5,30),'div_generator2':mp.nstr(alg.norm_coeff(od['div2']),30),'div_generator4':mp.nstr(alg.norm_coeff(od['div4']),30)}
    RR=alg.map_comp(od['Ptilde'],od['Ptilde'],5)
    path=BASE/'data'/f'critical_jet_transport_{dps}dps.json'; path.write_text(json.dumps(out,indent=2))
    print('linerr',out['linear_error_before_forcing']); print('sym',sym); print('evenratio',out['oddification_even_ratio'],'div',out['divY3'],out['divY5'],'gendiv',out['div_generator2'],out['div_generator4']); print('H4',out['H4_raw']); print('H6',out['H6_raw']); print(path)
    return out

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser(); ap.add_argument('--dps',type=int,default=70); ap.add_argument('--order',type=int,default=24); ap.add_argument('--hmax',default='.08')
 a=ap.parse_args(); summarize(a.dps,a.order,a.hmax)
