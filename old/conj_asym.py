import mpmath as mp
mp.mp.dps = 50
def theta_of(A, rho):
    c = rho*A; v = A-1; s = mp.sqrt(v*v-c*c); r = mp.mpf(rho)
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 600))
    f = lambda d: mp.zeta(2, A-d) - Phi
    lo, hi = mp.mpf(0), v-s
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) < 0: lo = mid
        else: hi = mid
    d = (lo+hi)/2
    return (v-d-s)/(v-s), d, v-s, Phi
print("--- theta-2/3 at fixed rho, growing A ---")
for rho in ['0.05','0.2','0.3333333333333333','0.5','0.7','0.9']:
    row = []
    for A in [mp.mpf(10)**3, mp.mpf(10)**4, mp.mpf(10)**5, mp.mpf(10)**6, mp.mpf(10)**7]:
        th, d, vs, Phi = theta_of(A, mp.mpf(rho))
        row.append(th - mp.mpf(2)/3)
    print(f"  rho={rho:8s}: " + "  ".join(mp.nstr(z,6).rjust(12) for z in row))
print()
print("--- A*(theta-2/3) : does it converge to a positive constant? ---")
for rho in ['0.05','0.2','0.3333333333333333','0.5','0.7','0.9']:
    row = []
    for A in [mp.mpf(10)**4, mp.mpf(10)**5, mp.mpf(10)**6, mp.mpf(10)**7]:
        th, d, vs, Phi = theta_of(A, mp.mpf(rho))
        row.append(A*(th - mp.mpf(2)/3))
    print(f"  rho={rho:8s}: " + "  ".join(mp.nstr(z,8).rjust(14) for z in row))
print()
print("--- expansion coefficients d1,d2,d3 (exact matching) at A=10 ---")
A = mp.mpf(10); z = lambda j: mp.zeta(j, A)
d1 = z(4)/(4*z(3))
d2 = (z(6)/3 - 3*z(4)*d1**2)/(2*z(3))
d3 = (z(8)/4 - 6*z(4)*d1*d2 - 4*z(5)*d1**3)/(2*z(3))
print("  d1 =", mp.nstr(d1,12), " d2 =", mp.nstr(d2,12), " d3 =", mp.nstr(d3,12))
for rho in ['0.1','0.3','0.5','0.7']:
    c = rho*A; v = A-1; s = mp.sqrt(v*v-c*c)
    th, d, vs, Phi = theta_of(A, mp.mpf(rho))
    print(f"   rho={rho}: delta={mp.nstr(d,12)}  d1c^2={mp.nstr(d1*c**2,12)}  d2c^4={mp.nstr(d2*c**4,12)}  (v-s)/3={mp.nstr(vs/3,12)}")
