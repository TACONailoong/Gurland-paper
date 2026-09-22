import mpmath as mp
mp.mp.dps = 35
print("--- (1) zeta4/zeta3 <= 2A^3 / (3(A-1/2)^3 (A+1)) ? ---")
worst = (mp.mpf(-9), None)
for k in range(0, 60):
    A = mp.mpf(3)/2 + mp.mpf(k)/10
    lhs = mp.zeta(4,A)/mp.zeta(3,A); rhs = 2*A**3/(3*(A-mp.mpf(1)/2)**3*(A+1))
    r = lhs/rhs
    if r > worst[0]: worst = (r, A)
print("   max ratio lhs/rhs =", mp.nstr(worst[0],10), "at A =", mp.nstr(worst[1],6))
print()
print("--- (2) Psi(A) = (A-1)A^3/((A-1/2)^3 (A+1)) <= 1 - 1/(3A) ? ---")
worst2 = (mp.mpf(-9), None)
for k in range(0, 3000):
    A = mp.mpf(2) + mp.mpf(k)/100
    P = (A-1)*A**3/((A-mp.mpf(1)/2)**3*(A+1)); b = 1-1/(3*A)
    r = P/b
    if r > worst2[0]: worst2 = (r, A)
print("   max Psi/(1-1/(3A)) =", mp.nstr(worst2[0],10), "at A =", mp.nstr(worst2[1],6), "  (need <= 1)")
print()
print("--- (3) R'(c) > 0 and R(c) > 0 on the theorem's range c^2 <= A/4 ---")
def chk(v, c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); Delta=c*c/(3*(v+s))
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,900))
    S   = sum(c**(2*k)*mp.zeta(2*k+2,A) for k in range(0,900))
    Rp = 2*mp.zeta(3,A-Delta)*c/(3*s) - 2*(S-Phi)/c     # R'(c)
    return mp.zeta(2,A-Delta)-Phi, Rp
okR = True; okRp = True; n = 0
for v in ['2','3','5','10','50','200','1000','1e4','1e5']:
    v = mp.mpf(v); cmax = mp.sqrt((1+v)/4)
    for j in range(1, 12):
        c = cmax*mp.mpf(j)/12
        if c >= v: continue
        R, Rp = chk(v, c); n += 1
        okR &= (R > 0); okRp &= (Rp > 0)
print(f"   {n} pairs with c^2 <= A/4:  R>0 everywhere: {okR};  R'>0 everywhere: {okRp}")
print()
print("--- (4) how far does R'>0 extend?  (threshold c^2 = kappa*A) ---")
for v in ['2','10','100','1e4']:
    v = mp.mpf(v); A = 1+v
    lo, hi = mp.mpf(0), v*mp.mpf('0.999')
    for _ in range(60):
        mid = (lo+hi)/2
        R, Rp = chk(v, mid)
        if Rp > 0: lo = mid
        else: hi = mid
    print(f"   v={mp.nstr(v,7):>9}: R'>0 up to c^2 = {mp.nstr(lo*lo,6)}  (= {mp.nstr(lo*lo/A,6)} * A)")
