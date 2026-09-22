import mpmath as mp
mp.mp.dps = 40
for (Ai, rho) in [('1e4','0.5'), ('1e6','0.5'), ('1e4','0.2'), ('1e4','0.9')]:
    A = mp.mpf(Ai); r = mp.mpf(rho); v = A-1; c = r*A
    s = mp.sqrt(v*v-c*c); Delta = (v-s)/3
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,400))
    f = lambda dd: mp.zeta(2,A-dd)-Phi
    lo,hi = mp.mpf(0), v-s
    for _ in range(300):
        m=(lo+hi)/2
        if f(m)<0: lo=m
        else: hi=m
    d=(lo+hi)/2
    R  = mp.zeta(2,A-Delta)-Phi
    LHS = 2*mp.quad(lambda z: mp.zeta(3,z), [A-Delta, A])
    print(f"A={Ai} rho={rho}: theta={mp.nstr(1-d/(v-s),12)}  delta/(v-s)={mp.nstr(d/(v-s),10)}")
    print(f"    R = {mp.nstr(R,10)}   LHS = 2*int_(A-Delta)^A zeta3 = {mp.nstr(LHS,10)}   Phi-zeta2 = {mp.nstr(Phi-mp.zeta(2,A),10)}   LHS-R = {mp.nstr(LHS-R,8)}")
