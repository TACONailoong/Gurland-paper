import mpmath as mp
mp.mp.dps = 40
def data(x, y):
    x = mp.mpf(x); y = mp.mpf(y)
    u = (x-y)/2; v = (x+y)/2; c = abs(u); s = mp.sqrt(x*y); A = 1+v
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 400))   # log G* / c^2
    # solve zeta(2, A - delta) = Phi for delta
    f = lambda d: mp.zeta(2, A-d) - Phi
    lo, hi = mp.mpf(0), v-s
    for _ in range(200):
        mid = (lo+hi)/2
        if f(mid) < 0: lo = mid
        else: hi = mid
    d = (lo+hi)/2
    theta = (v-d-s)/(v-s) if v != s else None
    Delta = c**2/(3*(v+s))      # the conjectured comparison point
    J = mp.zeta(2, A-Delta) - Phi
    return dict(v=v, c=c, s=s, A=A, rho=c/A, delta=d, theta=theta, Delta=Delta, J=J, rho_g=c/v)
print(" (x,y)         rho     theta       theta-2/3      J(Delta)      J/Phi")
worst = (9, None, None)
for (x, y) in [(1,2),(2,7),(20,100),(100,101),(0.1,199.9),(100,1),(1000,999),(1e4,1),
               (5,6),(1.9,2.1),(2,2.0001),(10,11),(0.5,3),(3,3.001),(50,60),(2,1000),
               (1e6,2e6),(1e8,1e8+1),(1.0001,1.0002)]:
    D = data(x,y)
    th = D['theta'] - mp.mpf(2)/3
    if th < worst[0]: worst = (th, (x,y), D['rho'])
    print(f" ({x},{y})".ljust(15), mp.nstr(D['rho'],5).ljust(8), mp.nstr(D['theta'],10).ljust(12),
          mp.nstr(th,4).ljust(12), mp.nstr(D['J'],6).ljust(13), mp.nstr(D['J']/(D['J']+ (mp.zeta(2,D['A'])-D['J'])),4))
print()
print("worst theta-2/3 over this sample:", mp.nstr(worst[0],6), "at", worst[1], "rho=", mp.nstr(worst[2],5))
print()
print("--- termwise comparison  K(u)  vs  1/(u-Delta)^2  ---")
for (x,y) in [(20,100),(1000,999),(1e6,2e6)]:
    D = data(x,y); c, A, Delta = D['c'], D['A'], D['Delta']
    bad = [n for n in range(0, 60) if mp.log(1/(1-c**2/(A+n)**2))/c**2 > 1/(A+n-Delta)**2]
    print(f"  (x,y)=({x},{y}) rho={mp.nstr(D['rho'],4)}: first {len(bad)} terms violate; u* threshold = {mp.nstr(Delta + mp.sqrt(Delta**2+c**2)/1,6)}")
