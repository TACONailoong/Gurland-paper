import mpmath as mp
mp.mp.dps = 40
def ratio(v,c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c)
    Phi=sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,700))
    f=lambda d: mp.zeta(2,A-d)-Phi
    lo,hi=mp.mpf(0),v-s
    for _ in range(300):
        m=(lo+hi)/2
        if f(m)<0: lo=m
        else: hi=m
    return (lo+hi)/2/(v-s)
print("--- is  delta/(v-s)  monotone DEcreasing in c (fixed v)? ---")
for v in ['1.5','2','5','20','100','1000','10000']:
    v=mp.mpf(v); vals=[]
    for j in range(1,13):
        c=v*mp.mpf(j)/13
        vals.append(ratio(v,c))
    mono = all(vals[i]>vals[i+1] for i in range(len(vals)-1))
    print(f"  v={mp.nstr(v,7):>9}: " + " ".join(mp.nstr(z,6) for z in vals) + f"   decreasing: {mono}")
print()
print("--- limit as c->0  should be 1/3 ---")
for v in ['2','5','100','10000']:
    v=mp.mpf(v)
    r=ratio(v, v*mp.mpf('1e-6'))
    print(f"  v={v}: delta/(v-s) at c=1e-6 v -> {mp.nstr(r,12)}  (1/3 = {mp.nstr(mp.mpf(1)/3,12)})")
