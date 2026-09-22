import mpmath as mp
mp.mp.dps = 40
def test(v, c, N=400000):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); Delta=(v-s)/3
    # Sigma = sum_{n>=0} 1/((A+n)^2 ((A+n)^2-c^2)), tail by integral
    n0 = int(A)+2
    S = sum(1/((A+n)**2*((A+n)**2-c*c)) for n in range(0, n0))
    # tail: 1/((A+n)^2-c^2) <= 1/((A+n)^2-A^2) = 1/((n-1)(2A+n-1)) for n>=1
    S += mp.quad(lambda t: 1/((A+t)**2*((A+t)**2-c*c)), [n0, mp.inf])
    lhs = mp.zeta(3, A-Delta)
    rhs = mp.mpf(3)/2*s*S
    return lhs, rhs, lhs/rhs
print("--- sharp chain:  zeta(3,A-Delta)  vs  (3s/2) * sum_n 1/(u^2(u^2-c^2)) ---")
worst=(mp.inf,None,None)
for vi in ['1.0001','1.01','1.1','1.5','2','3','5','10','50','100','1000','1e4','1e6']:
    v=mp.mpf(vi); mn=(mp.inf,None)
    for j in range(1,200):
        c=v*mp.mpf(j)/200
        lhs,rhs,r = test(v,c)
        if r<mn[0]: mn=(r,c)
    print(f"   v={mp.nstr(v,9):>10}:  min ratio LHS/RHS = {mp.nstr(mn[0],10):>12}  at c/v = {mp.nstr(mn[1]/v,5)}")
    if mn[0]<worst[0]: worst=(mn[0],v,mn[1])
print()
print("overall worst ratio:", mp.nstr(worst[0],10), "at v =", mp.nstr(worst[1],8), ", c/v =", mp.nstr(worst[2]/worst[1],5))
