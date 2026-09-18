from __future__ import annotations
import mpmath as mp
from multiprocessing import Pool
from pathlib import Path
import sys,json,time,math
BASE=Path('/mnt/data/carnot_resonance_paper')
sys.path.insert(0,str(BASE))
import extract_critical_jet as alg
MU='0.6628690069305101280936908638603636487'
YS='-2.5725893074769489665217743248705277'

def worker(task):
    comp,i,j,dps,addprec,order,hmax=task
    import mpmath as _mp,sys as _sys
    _mp.mp.dps=dps
    _sys.path.insert(0,str(BASE))
    from taylor_mp import local_return
    mu=_mp.mpf(MU); ys=_mp.mpf(YS)
    def fun(q,p):
        return local_return(mu,ys,[q,p],order=order,hmax=_mp.mpf(hmax),radius=_mp.mpf('.2'),with_dp=False,with_action=False,tmax=_mp.mpf('14'))[0][comp]
    der=_mp.diff(fun,(_mp.mpf('0'),_mp.mpf('0')),(i,j),addprec=addprec)
    coeff=der/(_mp.factorial(i)*_mp.factorial(j))
    return comp,i,j,_mp.nstr(coeff,dps)

def extract(dps=45,addprec=18,procs=8,order=30,hmax='.06'):
    tasks=[]
    for comp in [0,1]:
      for d in range(2,6):
        for i in range(d+1): tasks.append((comp,i,d-i,dps,addprec,order,hmax))
    t=time.time()
    with Pool(processes=procs) as pool:
        vals=list(pool.imap_unordered(worker,tasks,chunksize=1))
    print('derivatives',len(vals),'time',time.time()-t,flush=True)
    mp.mp.dps=dps
    P=[{(1,0):mp.mpf(-1)},{(0,1):mp.mpf(-1)}]
    for comp,i,j,s in vals: P[comp][(i,j)]=mp.mpf(s)
    # diagnostics raw symplectic
    det=alg.det_jac(P,5); det[(0,0)]=det.get((0,0),0)-1
    sym={str(d):mp.nstr(alg.norm_coeff(alg.homogeneous(det,d)),30) for d in range(1,5)}
    nf=alg.build_normal_form(P,5)
    H4=alg.integrate_hamiltonian_from_vector(nf['Y3'],4)
    H5=alg.integrate_hamiltonian_from_vector(nf['Y4'],5)
    H6=alg.integrate_hamiltonian_from_vector(nf['Y5'],6)
    ev2=alg.map_hom(nf['Ptilde'],2); ev4=alg.map_hom(nf['Ptilde'],4); od3=alg.map_hom(nf['Ptilde'],3); od5=alg.map_hom(nf['Ptilde'],5)
    even=max(alg.norm_coeff(ev2[0]),alg.norm_coeff(ev2[1]),alg.norm_coeff(ev4[0]),alg.norm_coeff(ev4[1])); odd=max(alg.norm_coeff(od3[0]),alg.norm_coeff(od3[1]),alg.norm_coeff(od5[0]),alg.norm_coeff(od5[1]))
    div3=alg.norm_coeff(alg.padd(alg.pder(nf['Y3'][0],0),alg.pder(nf['Y3'][1],1),7))
    div5=alg.norm_coeff(alg.padd(alg.pder(nf['Y5'][0],0),alg.pder(nf['Y5'][1],1),7))
    out={'dps':dps,'addprec':addprec,'P':[alg.dump_poly(P[0]),alg.dump_poly(P[1])],'symplectic_residual_by_degree':sym,'H4_raw':alg.dump_poly(H4),'H5_raw':alg.dump_poly(H5),'H6_raw':alg.dump_poly(H6),'Nmap':[alg.dump_poly(nf['Nmap'][0]),alg.dump_poly(nf['Nmap'][1])],'Ninv':[alg.dump_poly(nf['Ninv'][0]),alg.dump_poly(nf['Ninv'][1])],'oddification_even_ratio':mp.nstr(even/odd,30),'divY3':mp.nstr(div3,30),'divY5':mp.nstr(div5,30),'dens_odd':mp.nstr(max([abs(v) for k,v in nf['dens'].items() if sum(k)%2==1] or [mp.mpf(0)]),30)}
    RR=alg.map_comp(nf['R'],nf['R'],5); rr=alg.map_disp(RR,5)
    out['R_involution_error']=mp.nstr(max(alg.norm_coeff(rr[0]),alg.norm_coeff(rr[1])),30)
    path=BASE/'data'/f'critical_jet_diff_{dps}dps.json'; path.write_text(json.dumps(out,indent=2))
    print('sym',sym); print('evenratio',out['oddification_even_ratio'],'densodd',out['dens_odd'],'div',out['divY3'],out['divY5']); print('H4',out['H4_raw']); print('H5max',mp.nstr(alg.norm_coeff(H5),20)); print('H6',out['H6_raw']); print(path)
    return out

if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser(); ap.add_argument('--dps',type=int,default=45); ap.add_argument('--addprec',type=int,default=18); ap.add_argument('--procs',type=int,default=8); ap.add_argument('--order',type=int,default=30); ap.add_argument('--hmax',default='.06')
 a=ap.parse_args(); extract(a.dps,a.addprec,a.procs,a.order,a.hmax)
