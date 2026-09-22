import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40
def setup(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); w=(v+s)/2; a=1+v
    return v,c,s,w,a
print("="*80)
print("EM TERM BREAKDOWN (is there cancellation?)")
print("="*80)
for (x,y) in [(10000,0.001),(0.001,10000),(1,10),(100,1),(2,7),(1,2),(1000,999)]:
    v,c,s,w,a=setup(x,y)
    K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
    L=lambda t: 1/(t-w)**2
    G=lambda t: K(t)-L(t)
    I=mp.quad(G,[a,mp.inf])
    t1=G(a)/2; t2=-mp.diff(G,a,1)/12; t3=mp.diff(G,a,3)/720; t4=-mp.diff(G,a,5)/30240
    true=sum(G(n+v) for n in range(1,4001))+mp.quad(G,[4000+v,mp.inf])
    print(f" (x,y)=({x},{y})  rho={mp.nstr(c/a,6)}  w-rho_A... ")
    print(f"    I={mp.nstr(I,8):>14}  G(a)/2={mp.nstr(t1,8):>14}  -G'(a)/12={mp.nstr(t2,8):>14}")
    print(f"    G'''(a)/720={mp.nstr(t3,8):>14}  -G5/30240={mp.nstr(t4,8):>14}  true={mp.nstr(true,8):>14}")
    print(f"    |I|={mp.nstr(abs(I),4)}  sum|terms|={mp.nstr(abs(I)+abs(t1)+abs(t2)+abs(t3)+abs(t4),4)}   ratio={mp.nstr((abs(I)+abs(t1)+abs(t2)+abs(t3)+abs(t4))/abs(I),5)}")
print()
print("="*80)
print("c^4 COEFFICIENT (corrected: zeta_4^3, not zeta_4^2)")
print("="*80)
def tval(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    u=(x-y)/2; v=(x+y)/2; s=mp.sqrt(x*y)
    target=(mp.loggamma(1+x)+mp.loggamma(1+y)-2*mp.loggamma(1+v))/u**2
    f=lambda tt: mp.zeta(2,1+tt)-target
    lo,hi=s,v
    if not (f(lo)>0>f(hi)): return None
    for _ in range(220):
        mid=(lo+hi)/2
        if f(mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2
for (x,y) in [(1,2),(0.5,3),(2,7),(1.5,1.6),(10,11),(100,101)]:
    v,c,s,w,a=setup(x,y)
    t=tval(x,y)
    d=v-t
    z3=mp.zeta(3,a); z4=mp.zeta(4,a); z6=mp.zeta(6,a)
    lead=c**2*z4/(4*z3)
    corr=c**4/(2*z3)*(z6/3-3*z4**3/(16*z3**2))
    print(f"  (x,y)=({x},{y}): v-t={mp.nstr(d,12):>16} lead={mp.nstr(lead,12):>16} lead+corr={mp.nstr(lead+corr,12):>16} resid={mp.nstr(lead+corr-d,4)}")
print()
print("="*80)
print("CONJECTURE: infimum path")
print("="*80)
best=(mp.inf,None)
for V in ['1e2','1e3','1e4','1e5','1e6','1e7']:
    Vm=mp.mpf(V)
    for ei in ['1e-1','3e-2','1e-2','3e-3','1e-3','3e-4','1e-4','3e-5','1e-5']:
        e=mp.mpf(ei); x=Vm*(1+e); y=Vm*(1-e)
        t=tval(x,y)
        if t is None: continue
        s=mp.sqrt(x*y); th=(t-s)/(Vm-s)
        if th<best[0]: best=(th,(V,ei))
        print(f"   V={V:>8} eps={ei:>7}  theta={mp.nstr(th,12)}")
print("  BEST:", mp.nstr(best[0],12), "at", best[1], "  2/3 =", mp.nstr(mp.mpf(2)/3,12))
