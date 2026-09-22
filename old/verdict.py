import mpmath as mp
mp.mp.dps = 50
def probe(v, c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); Delta=(v-s)/3; r=c/A
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,600))
    R  = mp.zeta(2,A-Delta)-Phi
    # true theta
    f = lambda d: mp.zeta(2,A-d)-Phi
    lo,hi = mp.mpf(0), v-s
    for _ in range(300):
        m=(lo+hi)/2
        if f(m)<0: lo=m
        else: hi=m
    d=(lo+hi)/2
    th = (v-d-s)/(v-s)
    return r, R, th-mp.mpf(2)/3
print(" v=9999 (A=1e4): scan c, report R and theta-2/3")
for c in ['100','200','300','500','1000','2000','4000','8000']:
    r,R,dth = probe(9999, c)
    print(f"   c={c:>5} rho={mp.nstr(r,6):>10}  R={mp.nstr(R,8):>14}  theta-2/3={mp.nstr(dth,8):>14}   sign R: {'+' if R>0 else '-'}")
