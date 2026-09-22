import mpmath as mp
mp.mp.dps = 35
print("--- (a) g~(y) <= y/(2(1-y)^2) on (0,1) ? ---")
worst=(mp.mpf(-9),None)
for j in range(1,20000):
    y=mp.mpf(j)/20000
    g = 1/(1-y)+mp.log(1-y)/y
    b = y/(2*(1-y)**2)
    r = g/b
    if r>worst[0]: worst=(r,y)
print("   max g~/bound =", mp.nstr(worst[0],10), "at y =", mp.nstr(worst[1],5), " (need <= 1)")
print()
print("--- (b) chain:  E <= (c^2/2) Sigma_2  ?  and  zeta3(A-D) > (3s/2) Sigma_2 ? ---")
def chain(v,c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3
    E = sum(mp.mpf(m)/(m+1)*c**(2*m)*mp.zeta(2*m+2,A) for m in range(1,600))
    n0=int(A)+5
    S2 = sum(1/((A+n)**2-c*c)**2 for n in range(0,n0)) + 1/(3*(A+n0)**3)
    return E/(c*c/2*S2), mp.zeta(3,A-D)/(mp.mpf(3)/2*s*S2)
print("     v        c        E/((c^2/2)S2)      zeta3/((3s/2)S2)")
for v in ['2','5','20','100','1000','10000']:
    for f in ['0.05','0.3','0.7','0.95']:
        v0=mp.mpf(v); c0=v0*mp.mpf(f)
        if c0>=v0: continue
        r1,r2 = chain(v0,c0)
        flag = "OK" if (r1<=1 and r2>1) else "**FAIL**"
        print(f"   {v:>6} {mp.nstr(c0,6):>9}   {mp.nstr(r1,8):>14}   {mp.nstr(r2,8):>14}   {flag}")
