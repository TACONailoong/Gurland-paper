import mpmath as mp
mp.mp.dps = 40
def W(v, c):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; rho2=(c/A)**2; s=mp.sqrt(v*v-c*c)
    R1 = 2*A**3/(3*(A-mp.mpf(1)/2)**3*(A+1))
    return 1/(3*s*R1) - mp.mpf(1)/2 - rho2/(1-rho2)
print("--- W(c) = 1/(3 s R1) - 1/2 - rho^2/(1-rho^2)  must be > 0 ---")
for v in ['1.0001','1.01','1.1','1.5','2','3','5','10','50','100','1000','1e4','1e6','1e9']:
    v=mp.mpf(v); A=1+v; worst=(mp.inf,None)
    N=4000
    for j in range(1,N):
        c = v*mp.mpf(j)/N
        w = W(v,c)
        if w < worst[0]: worst=(w,c)
    print(f"   v={mp.nstr(v,8):>10}:  min W = {mp.nstr(worst[0],10):>14}  at c/v = {mp.nstr(worst[1]/v,6)}")
