import mpmath as mp
mp.mp.dps = 30
def cond(v,c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3; r2=(c/A)**2
    z3=mp.zeta(3,A); z4=mp.zeta(4,A)
    lhs=2*D*z3
    rhs=D*D*z4 + c*c*z4/(2*(1-r2))
    return lhs/rhs
print("  ratio  (2D z3) / (D^2 z4 + c^2 z4/(2(1-rho^2)))   must exceed 1")
for v in ['2','3','5','10','30','100','1000','10000']:
    row=[]
    for f in ['0.05','0.2','0.4','0.6','0.8','0.95']:
        v0=mp.mpf(v); c0=v0*mp.mpf(f)
        if c0>=v0: continue
        r=cond(v0,c0)
        row.append(f"{f}:{mp.nstr(r,4)}")
    print(f"   v={v:>6}: " + "  ".join(row))
print()
print("  --- as c^2 (i.e. rho^2) grows, ratio falls; find the crossover c/v ---")
for v in ['5','20','100','1000','10000']:
    v0=mp.mpf(v); lo,hi=mp.mpf('0.001'),mp.mpf('0.999')
    def ok(c): return cond(v0,c)>1
    if not ok(lo*v0): print(f"   v={v}: fails already at c=0.001v"); continue
    if ok(hi*v0): print(f"   v={v}: holds for ALL c<v"); continue
    for _ in range(80):
        mid=(lo+hi)/2
        if ok(mid*v0): lo=mid
        else: hi=mid
    print(f"   v={v:>6}: holds for c/v < {mp.nstr(lo,6)}  (c^2/A = {mp.nstr((lo*v0)**2/(1+v0),6)})")
