import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40
def Gamm(x,y):
    v=(x+y)/2
    return mp.gamma(x)*mp.gamma(y)/mp.gamma(v)**2
def Lam(r): return (1+r)*mp.log(1+r)+(1-r)*mp.log(1-r)
print("="*72); print("COR 3.2 (beta bounds) verification"); print("="*72)
for (x,y) in [(1,2),(1,1),(0.5,3),(2,7),(1,10),(10,20)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2
    G=Gamm(x,y); B=mp.gamma(x)*mp.gamma(y)/mp.gamma(x+y)
    rho=abs(x-y)/(x+y)
    lo=mp.sqrt(mp.pi)*v**mp.mpf('-0.5')/2**(2*v-1)*mp.e**((x-y)**2/4*mp.zeta(2,v))
    hi=mp.sqrt(mp.pi)*(v+1)**mp.mpf('0.5')/(2**(2*v-1)*v)*mp.e**(v*Lam(rho))
    print(f"  (x,y)=({x},{y}) v={mp.nstr(v,4)} B={mp.nstr(B,10)}  lo={mp.nstr(lo,8)}  hi={mp.nstr(hi,8)}  lo<B<hi: {lo<B<hi}")
print()
print("="*72); print("COR 3.1 (Gurland bounds) verification"); print("="*72)
for (x,y) in [(1,2),(0.5,3),(2,7),(1,10),(0.1,199.9)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; rho=abs(x-y)/(x+y)
    G=Gamm(x,y); L=mp.log(G)
    lo=(x-y)**2/4*mp.zeta(2,v); hi=v*Lam(rho)
    print(f"  (x,y)=({x},{y}) logG={mp.nstr(L,10)}  lo={mp.nstr(lo,10)}  hi={mp.nstr(hi,10)}  ok={lo<=L<=hi}   logG>(x-y)^2/(2(x+y)): {L>(x-y)**2/(2*(x+y))}")
print()
print("="*72); print("MULTIVARIATE 3-point example: x1=x2=v(1+e), x3=v(1-2e)"); print("="*72)
for vv in ['1','3']:
    v=mp.mpf(vv); e=mp.mpf('0.1')
    x1=v*(1+e); x2=v*(1+e); x3=v*(1-2*e)
    lhs=mp.log(mp.gamma(x1))+mp.log(mp.gamma(x2))+mp.log(mp.gamma(x3))-3*mp.log(mp.gamma(v))
    ser=sum(((-1)**j)/j*mp.zeta(j,v)*(v**j*e**j*(2+(-2)**j)) for j in range(2,200))
    print(f"  v={vv} eps={e}:  exact={mp.nstr(lhs,14)}  series={mp.nstr(ser,14)}  diff={mp.nstr(lhs-ser,4)}")
    print(f"     first terms: 3v^2e^2 zeta(2,v)={mp.nstr(3*v**2*e**2*mp.zeta(2,v),10)}  2v^3e^3 zeta(3,v)={mp.nstr(2*v**3*e**3*mp.zeta(3,v),10)}  (9/2)v^4e^4 zeta(4,v)={mp.nstr(mp.mpf(9)/2*v**4*e**4*mp.zeta(4,v),10)}")
print()
print("="*72); print("(x,y)=(100,1): Q>1 but rho<1 ; tails and optimality data"); print("="*72)
x=mp.mpf(100); y=mp.mpf(1); u=(x-y)/2; v=(x+y)/2
print("  Q =",mp.nstr(abs(x-y)/(2*(1+mp.sqrt(x*y))),6),"  rho =",mp.nstr(abs(x-y)/(2+x+y),8))
for m in [20,100,1000]:
    R=sum(mp.mpf(u**(2*k))/k*mp.zeta(2*k,1+v) for k in range(m,m+2500))
    print(f"   R_{m} = {mp.nstr(R,6)}")
x=mp.mpf(2); y=mp.mpf(7); u=(x-y)/2; v=(x+y)/2; rho=abs(x-y)/(2+x+y)
print("  optimality at (2,7): 2log rho =", mp.nstr(2*mp.log(rho),10))
for m in [5,10,20,40,80,160]:
    R=sum(mp.mpf(u**(2*k))/k*mp.zeta(2*k,1+v) for k in range(m,m+2500))
    print(f"    m={m:>4}: (1/m)log R_m = {mp.nstr(mp.log(R)/m,10)}")
