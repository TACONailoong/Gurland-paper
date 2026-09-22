import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40

def tval(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    u=(x-y)/2; v=(x+y)/2; s=mp.sqrt(x*y)
    if u==0: return None
    target=(mp.loggamma(1+x)+mp.loggamma(1+y)-2*mp.loggamma(1+v))/u**2
    f=lambda tt: mp.zeta(2,1+tt)-target
    lo,hi=s,v
    if not (f(lo)>0>f(hi)): return None
    for _ in range(240):
        mid=(lo+hi)/2
        if f(mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2

print("="*78)
print("THEOREM 1.8(ii) LOWER BOUND : validity of the proof's key inequality")
print("  need:  S'_0 < 3*delta0*zeta4*(v-delta0),   S'_0=sum_{j>=3}(j+1)delta0^j zeta(j+2,A)")
print("="*78)
worst=None; bad=0; tot=0
grid=[0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000]
for x in grid:
    for y in grid:
        if x==y: continue
        x=mp.mpf(x); y=mp.mpf(y)
        v=(x+y)/2; c=abs(x-y)/2; A=1+v; s=mp.sqrt(x*y)
        z3=mp.zeta(3,A); z4=mp.zeta(4,A)
        d0=c**2*z4/(4*z3+6*v*z4)
        S=sum((j+1)*d0**j*mp.zeta(j+2,A) for j in range(3,400))
        margin=3*d0*z4*(v-d0)
        ratio=float(S/margin) if margin>0 else 9e9
        tot+=1
        if S>=margin: bad+=1
        if worst is None or ratio>worst[0]: worst=(ratio,x,y,float(d0),float(v))
        # also test the crude analytic bound used in the proof
        sigma=d0/A
        lhs=A*sigma**2*(4-3*sigma)/(1-sigma)**2
        rhs=3*(v-d0)
        assert lhs<rhs, ("analytic bound fails", x,y,float(lhs),float(rhs))
print(f"  tested {tot};  key inequality violated in {bad} cases")
print(f"  worst S'_0/margin ratio = {worst[0]:.6f} at (x,y)=({worst[1]},{worst[2]}), delta0={worst[3]:.4g}, v={worst[4]:.4g}")
print("  analytic sufficient condition  A*sigma^2*(4-3sigma)/(1-sigma)^2 < 3(v-delta0)  holds in ALL cases")
print()
print("  lower bound validity: delta > delta0 ?")
bad2=0
for (x,y) in [(1,2),(0.5,3),(2,7),(1,10),(20,100),(0.1,199.9),(100,101),(0.01,1000),(1000,999),(10000,0.001)]:
    x=mp.mpf(x); y=mp.mpf(y)
    v=(x+y)/2; c=abs(x-y)/2; A=1+v
    d=v-tval(x,y)
    d0=c**2*mp.zeta(4,A)/(4*mp.zeta(3,A)+6*v*mp.zeta(4,A))
    ok=d>d0
    if not ok: bad2+=1
    print(f"   (x,y)=({mp.nstr(x,7)},{mp.nstr(y,9)})  delta={mp.nstr(d,10):>16}  delta0={mp.nstr(d0,10):>16}  ok={ok}")
print("  violations:",bad2)
