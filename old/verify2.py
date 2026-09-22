import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 60
def logGs(x,y):
    v=(x+y)/2
    return mp.loggamma(1+x)+mp.loggamma(1+y)-2*mp.loggamma(1+v)
def tval(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    if x==y: return None
    u=(x-y)/2; v=(x+y)/2; s=mp.sqrt(x*y)
    target=logGs(x,y)/u**2
    f=lambda tt: mp.zeta(2,1+tt)-target
    lo,hi=s,v
    if not (f(lo)>0>f(hi)): return None
    for _ in range(240):
        mid=(lo+hi)/2
        if f(mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2

print("="*76)
print("CHECK T6(ii): bracket  lo <= v-t <= hi   (with rho, not eta)")
print("="*76)
allok=True
for (x,y) in [(1,2),(1,10),(0.5,3),(2,7),(0.3,3),(20,100),(0.1,199.9),(5,20),(100,101),(3,mp.mpf('3.0001'))]:
    x=mp.mpf(x); y=mp.mpf(y)
    u=(x-y)/2; v=(x+y)/2; c=abs(u); A=1+v
    rho=abs(x-y)/(2+x+y)
    t=tval(x,y); d=v-t
    lo=c**2*mp.zeta(4,A)/(4*mp.zeta(3,A)+6*v*mp.zeta(4,A))
    hi=c**2*mp.zeta(4,A)/(4*mp.zeta(3,A))*(1+2*rho**2/(3*(1-rho**2)))
    pred=c**2*mp.zeta(4,A)/(4*mp.zeta(3,A))
    good = lo<=d<=hi
    allok &= good
    print(f"  ({mp.nstr(x,7)},{mp.nstr(y,9)}) rho={mp.nstr(rho,6):>9}  v-t={mp.nstr(d,10):>16}  lo={mp.nstr(lo,10):>16}  pred={mp.nstr(pred,10):>16}  hi={mp.nstr(hi,10):>16} {'ok' if good else 'FAIL'}")
    print(f"        (v-t)/(v-s) = {mp.nstr(d/(v-mp.sqrt(x*y)),10)}   pred/(v-s) = {mp.nstr(pred/(v-mp.sqrt(x*y)),10)}")
print("  ALL BRACKETS HOLD:", allok)

print()
print("="*76)
print("CHECK: (x,y)=(100,1) has Q=4.5>1 but rho<1 and series converges")
print("="*76)
x=mp.mpf(100); y=mp.mpf(1); u=(x-y)/2; v=(x+y)/2
print("  Q =",mp.nstr(abs(x-y)/(2*(1+mp.sqrt(x*y))),6),"  rho =",mp.nstr(abs(x-y)/(2+x+y),8))
for m in [20,100,1000]:
    R=sum(mp.mpf(u**(2*k))/k*mp.zeta(2*k,1+v) for k in range(m,m+2500))
    print(f"   R_{m} = {mp.nstr(R,6)}")

print()
print("="*76)
print("CHECK T2 optimality: (1/m) log R_m -> 2 log rho   at (x,y)=(2,7)")
print("="*76)
x=mp.mpf(2); y=mp.mpf(7); u=(x-y)/2; v=(x+y)/2; rho=abs(x-y)/(2+x+y)
print("  rho =",mp.nstr(rho,10),"  2log rho =",mp.nstr(2*mp.log(rho),10))
for m in [5,10,20,40,80,160]:
    R=sum(mp.mpf(u**(2*k))/k*mp.zeta(2*k,1+v) for k in range(m,m+2500))
    print(f"    m={m:>4}: (1/m)log R_m = {mp.nstr(mp.log(R)/m,10)}")

print()
print("="*76)
print("TABLE 3 data: theta and small-c prediction")
print("="*76)
rows=[]
for (x,y) in [(1,2),(0.5,3),(2,7),(1,10),(5,20),(20,100),(0.1,100),(0.1,199.9),(3,mp.mpf('3.0001')),(1,mp.mpf('1.0001')),(10,10000)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; s=mp.sqrt(x*y)
    t=tval(x,y); th=(t-s)/(v-s); vt=(v-t)/(v-s)
    pred=1-(v*mp.zeta(4,1+v))/(2*mp.zeta(3,1+v))
    rows.append((x,y,th,vt,pred))
    print(f"  x={mp.nstr(x,8):>10} y={mp.nstr(y,10):>12}  theta={mp.nstr(th,10):>13}  (v-t)/(v-s)={mp.nstr(vt,10):>13}  pred={mp.nstr(pred,10)}")

print()
print("="*76)
print("CONJECTURE scan: minimum theta over a wide grid (is it > 2/3 ?)")
print("="*76)
mn=mp.inf; arg=None
vals=[0.001,0.01,0.1,0.5,1,2,5,10,20,50,100,500,1000,10000]
for x in vals:
    for y in vals:
        if x==y: continue
        xm=mp.mpf(x); ym=mp.mpf(y); v=(xm+ym)/2; s=mp.sqrt(xm*ym)
        t=tval(xm,ym)
        if t is None: continue
        th=(t-s)/(v-s)
        if th<mn: mn=th; arg=(x,y)
print("  min theta =", mp.nstr(mn,12), " at (x,y) =", arg, "   2/3 =", mp.nstr(mp.mpf(2)/3,12))
print("  conjecture holds on grid:", mn > mp.mpf(2)/3)
