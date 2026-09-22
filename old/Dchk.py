import mpmath as mp
mp.mp.dps = 30
def D(r):
    r = mp.mpf(r); t = (1-mp.sqrt(1-r*r))/3
    S0 = 2*mp.atanh(r)/r + mp.log(1-r*r)/(r*r)
    return t/(1-t) - (S0-1)
print("--- D(rho) = t/(1-t) - (Sigma_0(rho)-1) > 0 ? ---")
worst=(mp.inf,None)
for j in range(1,2000):
    r = mp.mpf(j)/2000
    d = D(r)
    if d < worst[0]: worst=(d,r)
print("   min D =", mp.nstr(worst[0],8), "at rho =", mp.nstr(worst[1],5))
print("   D(0.01)=", mp.nstr(D('0.01'),8), "  D(0.2)=", mp.nstr(D('0.2'),8), "  D(0.5)=", mp.nstr(D('0.5'),8), "  D(0.9)=", mp.nstr(D('0.9'),8), "  D(0.99)=", mp.nstr(D('0.99'),8))
print("   series check: rho^4/360 + 11 rho^6/3024 at rho=0.2:", mp.nstr(mp.mpf('0.2')**4/360 + 11*mp.mpf('0.2')**6/3024, 8), " vs D:", mp.nstr(D('0.2'),8))
