import mpmath as mp
mp.mp.dps = 40
def R_of(v, c):
    v = mp.mpf(v); c = mp.mpf(c)
    A = 1+v; s = mp.sqrt(v*v-c*c)
    Delta = c*c/(3*(v+s))
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 900))
    return mp.zeta(2, A-Delta) - Phi
print("--- is R(c) increasing in c for fixed v ?  (R(0)=0) ---")
for v in ['1.5', '2', '3', '5', '10', '100', '1000']:
    v = mp.mpf(v)
    cs = [v*mp.mpf(k)/10 for k in range(1, 10)]
    vals = [R_of(v, c) for c in cs]
    mono = all(vals[i] < vals[i+1] for i in range(len(vals)-1))
    print(f"  v={mp.nstr(v,6):>8}: R = " + " ".join(mp.nstr(z,3).rjust(10) for z in vals) + f"   monotone: {mono}")
print()
print("--- and J > 0 ?  R should be positive ---")
for v in ['1.5','2','10','1000']:
    v = mp.mpf(v)
    print(f"  v={v}: " + " ".join(mp.nstr(R_of(v, v*mp.mpf(k)/5),3).rjust(11) for k in range(1,5)))
