import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40
def Lam(r): return (1+r)*mp.log(1+r)+(1-r)*mp.log(1-r)
print("="*74)
print("CORRECTED COR 3.1 :  (x-y)^2/4 zeta(2,v) <= log G <= -log(1-rho^2) + v Lambda(rho)")
print("="*74)
ok=True
for (x,y) in [(1,2),(0.5,3),(2,7),(1,10),(0.1,199.9),(1,1.0001),(100,101),(0.01,1000)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; rho=abs(x-y)/(x+y)
    G=mp.gamma(x)*mp.gamma(y)/mp.gamma(v)**2; L=mp.log(G)
    lo=(x-y)**2/4*mp.zeta(2,v); hi=v*Lam(rho)-mp.log(1-rho**2)
    good=(lo<=L<=hi); ok&=good
    print(f"  (x,y)=({x},{y})  logG={mp.nstr(L,10):>14}  lo={mp.nstr(lo,10):>14}  hi={mp.nstr(hi,10):>14}  ok={good}")
print("  ALL:",ok)
print()
print("="*74)
print("CORRECTED COR 3.2 (beta bounds)")
print("="*74)
ok=True
for (x,y) in [(1,2),(1,1),(0.5,3),(2,7),(1,10),(10,20),(0.3,3)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; rho=abs(x-y)/(x+y)
    B=mp.gamma(x)*mp.gamma(y)/mp.gamma(x+y)
    lo=mp.sqrt(mp.pi)*v**mp.mpf('-0.5')/2**(2*v-1)*mp.e**((x-y)**2/4*mp.zeta(2,v))
    hi=mp.sqrt(mp.pi)*(v+1)**mp.mpf('0.5')/(2**(2*v-1)*v)*mp.e**(v*Lam(rho)-mp.log(1-rho**2))
    good=(lo<B<hi); ok&=good
    print(f"  (x,y)=({x},{y}) v={mp.nstr(v,5)}  B={mp.nstr(B,10):>15}  lo={mp.nstr(lo,10):>15}  hi={mp.nstr(hi,10):>15}  ok={good}")
print("  ALL:",ok)
print()
print("example (1,2): B=0.5, lo=%.6f hi=%.6f ; G(1,2)=%.6f, Gs(1,2)=%.6f" % (
   float(mp.sqrt(mp.pi)*mp.mpf(1.5)**mp.mpf('-0.5')/2**2*mp.e**mp.mpf(0.25)*0+0) if False else float(mp.sqrt(mp.pi)*mp.mpf('1.5')**mp.mpf('-0.5')/2**2*mp.e**(mp.mpf(1)/4*mp.zeta(2,mp.mpf('1.5')))),
   float(mp.sqrt(mp.pi)*mp.mpf('2.5')**mp.mpf('0.5')/(2**2*mp.mpf('1.5'))*mp.e**(mp.mpf('1.5')*Lam(mp.mpf(1)/3)-mp.log(1-mp.mpf(1)/9))),
   float(mp.gamma(1)*mp.gamma(2)/mp.gamma(mp.mpf('1.5'))**2),
   float(mp.gamma(2)*mp.gamma(3)/mp.gamma(mp.mpf('2.5'))**2)))
