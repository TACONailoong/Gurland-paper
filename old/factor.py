import mpmath as mp
mp.mp.dps = 40
def scan(v, c, K=400000):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3
    thr = (c*c+D*D)/(2*D)
    def term(u):
        return c*c/(u-D)**2 + mp.log(1-c*c/(u*u))   # log of the factor
    deficit=mp.mpf(0); surplus=mp.mpf(0); nbad=0
    m=0
    while True:
        u=A+m
        t=term(u)
        if t<0: deficit += -t; nbad+=1
        else: surplus += t
        # stop when u well beyond threshold and terms are tiny
        if u>thr+50 and t < mp.mpf(10)**(-30): break
        m+=1
        if m>K: break
    # tail beyond: all good, add integral estimate
    u0=A+m
    surplus += 2*c*c*D/(2*u0*u0)   # int_u0^inf 2 c^2 D u^-3 du
    return thr-A, nbad, deficit, surplus, m
print("     v        c        thr-A      #bad     deficit        surplus      S/D")
for v in ['5','20','100','1000','10000']:
    for f in ['0.1','0.5','0.9']:
        v0=mp.mpf(v); c0=v0*mp.mpf(f)
        if c0>=v0: continue
        t0,nb,de,su,m = scan(v0,c0)
        print(f"   {v:>6} {mp.nstr(c0,6):>9} {mp.nstr(t0,6):>10} {nb:>6} {mp.nstr(de,6):>14} {mp.nstr(su,6):>14}  {mp.nstr(su/de,6) if de>0 else 'inf'}")
