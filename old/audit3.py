import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40

def setup(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); w=(v+s)/2; a=1+v
    return v,c,s,w,a

def test(x,y):
    v,c,s,w,a = setup(x,y)
    K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
    L=lambda t: 1/(t-w)**2
    G=lambda t: K(t)-L(t)
    N=4000
    true=sum(G(n+v) for n in range(1,N+1))+mp.quad(G,[N+v,mp.inf])
    I=mp.quad(G,[a,mp.inf])
    em=I+G(a)/2
    terms=[I,G(a)/2]
    for j,co in [(1,mp.mpf(-1)/12),(3,mp.mpf(1)/720),(5,mp.mpf(-1)/30240)]:
        val=co*mp.diff(G,a,j); em+=val; terms.append(val)
    return float(true),float(em),float(I)

print("="*76)
print("EULER-MACLAURIN VIABILITY:  sum_{n>=1}G(n+v) < 0 ?   G=K-L (midpoint target)")
print("="*76)
grid=[0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000]
worst=[];nbad=0;tot=0
for xv in grid:
    for yv in grid:
        if xv==yv: continue
        try: true,em,I=test(xv,yv)
        except Exception as e: continue
        tot+=1
        if em>=0: nbad+=1
        worst.append((em,true,xv,yv))
worst.sort(reverse=True)
print(f"  tested {tot}; EM partial sum >= 0 in {nbad} cases")
print("  worst (largest) EM partial sums:")
for w in worst[:10]:
    print(f"     EM={w[0]:+.5e}  true={w[1]:+.5e}   (x,y)=({w[2]},{w[3]})")

print()
print("="*76)
print("CONJECTURE fine scan: targeted extremal family  x=v(1+e), y=v(1-e)")
print("="*76)
def tval(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    if x==y: return None
    u=(x-y)/2; v=(x+y)/2; s=mp.sqrt(x*y)
    target=(mp.loggamma(1+x)+mp.loggamma(1+y)-2*mp.loggamma(1+v))/u**2
    f=lambda tt: mp.zeta(2,1+tt)-target
    lo,hi=s,v
    if not (f(lo)>0>f(hi)): return None
    for _ in range(200):
        mid=(lo+hi)/2
        if f(mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2
mn=mp.inf;arg=None
for V in ['10','100','1000','10000','100000']:
    V=mp.mpf(V)
    for ei in ['0.001','0.005','0.01','0.05','0.1','0.3','0.5','0.7','0.9','0.99','0.999']:
        e=mp.mpf(ei); x=V*(1+e); y=V*(1-e)
        t=tval(x,y)
        if t is None: continue
        v=V; s=mp.sqrt(x*y)
        th=(t-s)/(v-s)
        if th<mn: mn=th; arg=(str(V),ei)
print("  minimum theta on this family:", mp.nstr(mn,12), " at (V,eps)=", arg)
print("  2/3 =", mp.nstr(mp.mpf(2)/3,12), "  holds:", mn>mp.mpf(2)/3)

print()
print("="*76)
print("CHECK new c^4 coefficient in the expansion of v-t")
print("="*76)
for (x,y) in [(1,2),(0.5,3),(2,7),(1.5,1.6),(10,11)]:
    v,c,s,w,a=setup(x,y)
    c=abs((mp.mpf(x)-mp.mpf(y))/2)
    t=tval(x,y); 
    if t is None: continue
    d=v-t
    z3=mp.zeta(3,a); z4=mp.zeta(4,a); z6=mp.zeta(6,a)
    lead=c**2*z4/(4*z3)
    corr=c**4/(2*z3)*(z6/3-3*z4**2/(16*z3**2))
    print(f"  (x,y)=({x},{y}): v-t={mp.nstr(d,12):>16}  lead={mp.nstr(lead,12):>16}  lead+corr={mp.nstr(lead+corr,12):>16}  lead+corr-(v-t)={mp.nstr(lead+corr-d,4)}")
