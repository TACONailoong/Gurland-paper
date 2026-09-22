import mpmath as mp
mp.mp.dps = 45
def theta_of(V, eps):
    V = mp.mpf(V); e = mp.mpf(eps)
    x = V*(1+e); y = V*(1-e)
    v = V; c = V*e; A = 1+v; s = mp.sqrt(x*y)
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 800))
    f = lambda d: mp.zeta(2, A-d) - Phi
    lo, hi = mp.mpf(0), v-s
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) < 0: lo = mid
        else: hi = mid
    d = (lo+hi)/2
    return 1 - (d)/(v-s)
print("--- fixed eps: theta at growing V (asymptote depends on eps!) ---")
print("      V        eps=1e-1        eps=1e-3        eps=1e-5")
for V in [mp.mpf(10)**3, mp.mpf(10)**4, mp.mpf(10)**5, mp.mpf(10)**6, mp.mpf(10)**8, mp.mpf(10)**10]:
    row = [theta_of(V, e) for e in ('0.1','0.001','0.00001')]
    print(f"  10^{int(mp.log10(V)):>2d}  " + "  ".join(mp.nstr(r,12).rjust(15) for r in row))
print()
print("--- joint limit eps = V^{-2}: theta -> 2/3 ---")
for V in [mp.mpf(10)**k for k in (2,3,4,5,6,7,8)]:
    e = V**-2
    th = theta_of(V, e)
    print(f"  V=10^{int(mp.log10(V)):>2d} (eps={mp.nstr(e,3)}): theta={mp.nstr(th,15)}  theta-2/3={mp.nstr(th-mp.mpf(2)/3,4)}  V*(theta-2/3)={mp.nstr(V*(th-mp.mpf(2)/3),6)}")
print()
print("--- joint limit eps = V^{-3/2} ---")
for V in [mp.mpf(10)**k for k in (3,4,5,6,7,8)]:
    e = V**mp.mpf('-1.5')
    th = theta_of(V, e)
    print(f"  V=10^{int(mp.log10(V)):>2d} (eps={mp.nstr(e,3)}): theta-2/3={mp.nstr(th-mp.mpf(2)/3,4)}")
