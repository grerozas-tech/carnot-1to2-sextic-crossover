import sys, csv
import mpmath as mp
sys.path.insert(0,'/mnt/data/carnot_resonance_paper')
import taylor_mp as tm

mp.mp.dps=50
A=mp.mpf('203.67235802390137'); rho=mp.mpf('2.0011829182526856'); alpha=mp.mpf('0.0946271414514154')
rows=[]
cases=[('A',mp.mpf('.020'),mp.mpf('0.6628074587750346')),('A',mp.mpf('.010'),mp.mpf('0.6628383082964905')),('A',mp.mpf('.005'),mp.mpf('0.662853676522959')),('B',mp.mpf('.016'),mp.mpf('0.6628198166401073')),('B',mp.mpf('.008'),mp.mpf('0.662844460119181')),('B',mp.mpf('.004'),mp.mpf('0.6628567456330406'))]
a0=mp.mpf('-2.5725893'); T0=mp.mpf('3.017955')
for ladder,dnom,mu in cases:
    a,T,it=tm.mother(mu,a0,T0,order=24,hmax=mp.mpf('.1'),tol=mp.mpf('1e-38')); a0,T0=a,T
    st,DPm=tm.full_DP(mu,a,T,order=24,hmax=mp.mpf('.1'))
    trm=DPm[0,0]+DPm[1,1]; delta=mp.acos(-trm/2)
    F0,D0,AF0,z10,t0=tm.fmap_local(mu,a,[0,0],order=24,hmax=mp.mpf('.1'),with_dp=True,with_action=True)
    amp=mp.sqrt(delta/(4*A))
    seeds={'Q':[mp.sqrt(alpha)*amp,0], 'P':[0,amp/mp.sqrt(alpha)]}
    out={}
    for br,seed in seeds.items():
        z,nit=tm.fixed_newton(mu,a,seed,order=24,hmax=mp.mpf('.1'),tol=mp.mpf('1e-32'),maxit=8)
        F,DF,AF,z1,ts=tm.fmap_local(mu,a,z,order=24,hmax=mp.mpf('.1'),with_dp=True,with_action=True)
        tr=DF[0,0]+DF[1,1]; det=mp.det(DF); k=mp.acosh(tr/2); J=(AF-AF0)/2
        res=max(abs(F[0]-z[0]),abs(F[1]-z[1]))
        out[br]=(z,z1,J,k,tr,det,res)
        Y=(J+delta**2/(16*A))/delta**3
        Lam=(k-2*rho*delta)/delta**2
        rows.append([ladder,mp.nstr(dnom,20),mp.nstr(delta,30),mp.nstr(mu,30),br,mp.nstr(z[0],30),mp.nstr(z[1],30),mp.nstr(z1[0],30),mp.nstr(z1[1],30),mp.nstr(J,30),mp.nstr(k,30),mp.nstr(tr,30),mp.nstr(det-1,12),mp.nstr(res,12),mp.nstr(Y,25),mp.nstr(Lam,25)])
    Jm=(out['Q'][2]+out['P'][2])/2; km=(out['Q'][3]+out['P'][3])/2
    Aorb=-delta**2/(16*Jm); rhoorb=km/(2*delta)
    Yq=(out['Q'][2]+delta**2/(16*A))/delta**3; Yp=(out['P'][2]+delta**2/(16*A))/delta**3
    Lq=(out['Q'][3]-2*rho*delta)/delta**2; Lp=(out['P'][3]-2*rho*delta)/delta**2
    Yplus=(Yq+Yp)/2; Yminus=(Yq-Yp)/2; Lplus=(Lq+Lp)/2; Lminus=(Lq-Lp)/2
    Geff=512*A**3*Yplus-8*A**2*rho*Lplus
    Deff=512*A**3*Yminus-8*A**2*rho*Lminus
    print('ladder',ladder,'dnom',dnom,'delta',mp.nstr(delta,18),'mu',mp.nstr(mu,18))
    print('  Aorb',mp.nstr(Aorb,18),'rhoorb',mp.nstr(rhoorb,18),'rho-2',mp.nstr(rhoorb-2,12))
    print('  JQ JP',mp.nstr(out['Q'][2],18),mp.nstr(out['P'][2],18))
    print('  kQ kP',mp.nstr(out['Q'][3],18),mp.nstr(out['P'][3],18),'diff',mp.nstr(out['Q'][3]-out['P'][3],12))
    print('  Y+',mp.nstr(Yplus,14),'L+',mp.nstr(Lplus,14),'Geff',mp.nstr(Geff,14),'Deff',mp.nstr(Deff,14))
    rows.append([ladder,mp.nstr(dnom,20),mp.nstr(delta,30),mp.nstr(mu,30),'MEAN','','','','',mp.nstr(Jm,30),mp.nstr(km,30),'','','',mp.nstr(Aorb,25),mp.nstr(rhoorb,25),mp.nstr(Geff,25),mp.nstr(Deff,25)])

path='/mnt/data/carnot_resonance_paper/data/axial_reconstruction_mp.csv'
with open(path,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['ladder','delta_nominal','delta_actual','mu','branch','q','p','Pq','Pp','J_reduced','kappa_F','trace_DF','det_DF_error','fixed_residual','Y_or_Aorb','Lambda_or_rhoorb','G_effective','Delta6_effective']); w.writerows(rows)
print('WROTE',path)
