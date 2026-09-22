import mpmath as mp
mp.mp.dps = 30
print("--- (1) zeta(3,z) >= 1/(2z^2) ---")
w = min(mp.zeta(3, z) * 2 * z**2 for z in [mp.mpf(10)**k for k in range(-2, 5)] + [mp.mpf('1.5'), 3, 7])
print("   worst ratio zeta(3,z)*2z^2 =", mp.nstr(w, 10))
print("--- (2) (v+s)(1+s) <= 2v(1+v) with s = sqrt(xy) ---")
ok = True
for (x, y) in [(1,1),(1,2),(0.1,199.9),(2,7),(100,1),(1e4,1),(3,3)]:
    x = mp.mpf(x); y = mp.mpf(y); v = (x+y)/2; s = mp.sqrt(x*y)
    ok &= (v+s)*(1+s) <= 2*v*(1+v) + mp.mpf(10)**-25
print("   holds:", ok)
print("--- (3) v(A+3)/(3A^2) <= 4/9 ---")
worst = max((v*(v+4)/(3*(v+1)**2), v) for v in [mp.mpf(k)/10 for k in range(1, 2000)] + [mp.mpf(2)])
print("   max =", mp.nstr(worst[0], 12), "at v =", worst[1])
print("--- (4) corr < gain whenever rho < sqrt5/3 ---")
thr = mp.sqrt(5)/3
print("   threshold rho* =", mp.nstr(thr, 12))
lhs = []
for v0 in [mp.mpf(s) for s in ('0.25','0.5','1','2','5','20','100','1000','1e5')]:
    for r0 in [mp.mpf(j)/100 for j in range(1, 100)]:
        c = v0*r0; A = 1+v0; s = v0*mp.sqrt(1-r0**2); rho = c/A
        if rho >= thr: continue
        gain = c**2*(mp.zeta(2,1+s)-mp.zeta(2,1+v0))
        corr = rho**4/(2*(1-rho**2))*(1+A/3)
        lhs.append((corr/gain, v0, r0, rho))
print("   worst corr/gain among rho<rho*:", mp.nstr(max(lhs)[0], 8), "at v =", max(lhs)[1], "shape", mp.nstr(max(lhs)[2],4))
# how far can the threshold be pushed?
worst_ok = mp.mpf(0)
for v0 in [mp.mpf(s) for s in ('0.25','0.5','1','2','5','20','100','1000','1e5','1e7')]:
    for j in range(1, 1000):
        r0 = mp.mpf(j)/1000
        c = v0*r0; A = 1+v0; s = v0*mp.sqrt(1-r0**2); rho = c/A
        gain = c**2*(mp.zeta(2,1+s)-mp.zeta(2,1+v0))
        corr = rho**4/(2*(1-rho**2))*(1+A/3)
        if corr < gain: worst_ok = max(worst_ok, rho)
print("   empirically, improvement holds up to rho =", mp.nstr(worst_ok, 8), "(so the rigorous rho* is not tight)")
