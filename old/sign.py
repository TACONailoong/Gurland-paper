import mpmath as mp
mp.mp.dps = 40
print("--- which way?  delta*c^2*zeta3(A-delta)   vs   s(v-s)E  ---")
def L(v,c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3
    Phi=sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,700))
    f=lambda d: mp.zeta(2,A-d)-Phi
    lo,hi=mp.mpf(0),v-s
    for _ in range(300):
        m=(lo+hi)/2
        if f(m)<0: lo=m
        else: hi=m
    d=(lo+hi)/2
    E=sum(mp.mpf(m)/(m+1)*c**(2*m)*mp.zeta(2*m+2,A) for m in range(1,600))
    return d*c*c*mp.zeta(3,A-d), s*(v-s)*E
for v in ['2','5','20','100','10000']:
    for f in ['0.1','0.5','0.9']:
        v0=mp.mpf(v); c0=v0*mp.mpf(f)
        a,b=L(v0,c0)
        print(f"  v={v:>6} c/v={f}:  delta c^2 z3 = {mp.nstr(a,6):>14}   s(v-s)E = {mp.nstr(b,6):>14}   {'>' if a>b else '<'}")
