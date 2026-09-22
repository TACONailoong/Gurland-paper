import mpmath as mp
mp.mp.dps = 30
# rigorous sufficient condition, keeping s exact:
#   T(c) := s (A-1/2)^3 (A+1) (1-rho^2) / (A^3 (1+rho^2))  >  1
def T(v, c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); rho2=(c/A)**2
    return s*(A-mp.mpf(1)/2)**3*(A+1)*(1-rho2)/(A**3*(1+rho2))
print("--- T(c) > 1 ?  (rigorous sufficient condition for R'>0) ---")
for v in ['2','3','5','10','100','1000','10000']:
    v=mp.mpf(v); A=1+v; cmax=v*mp.mpf('0.99999')
    # find threshold where T = 1
    lo, hi = mp.mpf('1e-12'), cmax
    if T(v,hi) > 1:
        print(f"   v={mp.nstr(v,7):>9}: T>1 for ALL c<v"); continue
    for _ in range(200):
        mid=(lo+hi)/2
        if T(v,mid) > 1: lo=mid
        else: hi=mid
    c1=lo
    # compare with the TRUE threshold where R'>0 fails
    def Rprime(c):
        c=mp.mpf(c); s=mp.sqrt(v*v-c*c); Delta=c*c/(3*(v+s))
        Phi = sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,700))
        S   = sum(c**(2*k)*mp.zeta(2*k+2,A) for k in range(0,700))
        return 2*mp.zeta(3,A-Delta)*c/(3*s) - 2*(S-Phi)/c
    lo2, hi2 = mp.mpf('1e-12'), cmax
    if Rprime(hi2) > 0:
        print(f"   v={mp.nstr(v,7):>9}: rigorous T>1 up to c={mp.nstr(c1,6)}; R'>0 for ALL c<v"); continue
    for _ in range(120):
        mid=(lo2+hi2)/2
        if Rprime(mid) > 0: lo2=mid
        else: hi2=mid
    print(f"   v={mp.nstr(v,7):>9}: rigorous up to c={mp.nstr(c1,6)} (c/v={mp.nstr(c1/v,5)}); true R'>0 up to c={mp.nstr(lo2,6)} (c/v={mp.nstr(lo2/v,5)})")
