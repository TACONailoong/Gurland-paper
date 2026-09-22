import mpmath as mp
mp.mp.dps = 45
def solve(A, c):
    v = A-1; s = mp.sqrt(v*v-c*c); rho = c/A
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 700))
    f = lambda d: mp.zeta(2, A-d) - Phi
    lo, hi = mp.mpf(0), v-s
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) < 0: lo = mid
        else: hi = mid
    d = (lo+hi)/2
    e1 = c**2*mp.zeta(4,A)/(4*mp.zeta(3,A))
    return dict(rho=rho, delta=d, e1=e1, ratio=d/e1, theta=(v-d-s)/(v-s),
                ds=(d/(v-s)), pred=(v+s)**2*mp.zeta(4,A)/(4*mp.zeta(3,A)))
print("--- sharpened bracket: is  delta/e1 - 1 = O(rho^2) ? ---")
print("   (x,y)            rho        delta/e1-1     /rho^2     1/(1-rho^2)")
for (x,y) in [(1000,999),(1e4,1e4-1),(100,99),(20,100),(1,2),(2,7),(0.1,199.9),(1e6,2e6)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; c=abs(x-y)/2; A=1+v
    D = solve(A,c)
    r2 = D['rho']**2
    print(f"  ({mp.nstr(x,8)},{mp.nstr(y,8)})".ljust(24), mp.nstr(D['rho'],5).ljust(10),
          mp.nstr(D['ratio']-1,5).ljust(13), mp.nstr((D['ratio']-1)/r2,5).ljust(11), mp.nstr(1/(1-r2),5))
print()
print("--- delta/(v-s) vs prediction (v+s)^2 zeta4/(4 zeta3) ; and theta ---")
for (x,y) in [(1e3,1e3-1),(1e5,1e5-1),(1e7,1e7-1),(1e3,2e3),(1e5,2e5),(1e7,2e7)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; c=abs(x-y)/2; A=1+v
    D = solve(A,c)
    print(f"  ({mp.nstr(x,6)},{mp.nstr(y,6)}) rho={mp.nstr(D['rho'],5)}: theta={mp.nstr(D['theta'],12)}"
          f"  delta/(v-s)={mp.nstr(D['ds'],12)}  pred={mp.nstr(D['pred'],12)}  pred-ratio={mp.nstr(D['pred']/D['ds']-1,4)}")
