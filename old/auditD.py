import mpmath as mp
mp.mp.dps = 30
def N(e,q):
    return (40*e**3+40*e**2*q-12*e**2+83*e*q**2-32*e*q-24*e
            +2*q**3+q**2-8*q-4)
bad=[]
for j in range(1,2000):
    e=mp.mpf(j)/2000; qhi=1-e
    # convexity in q and endpoint signs
    d2 = 166*e+12*0+2          # d2N/dq2 with q=0 (linear in q)
    n0=N(e,0); nh=N(e,qhi)
    for k in range(0,11):
        q=qhi*mp.mpf(k)/10
        if N(e,q)>max(n0,nh)+mp.mpf('1e-20'): bad.append((e,q))
print("convexity/endpoint argument violations:", bad[:3], " count:", len(bad))
print("N(e,0)   sign check:", mp.nstr(N(mp.mpf('0.3'),0),8), mp.nstr(N(mp.mpf('0.7'),0),8))
print("N(e,1-e) sign check:", mp.nstr(N(mp.mpf('0.3'),mp.mpf('0.7')),8), mp.nstr(N(mp.mpf('0.7'),mp.mpf('0.3')),8))
print("their factorisations:")
e,q=mp.symbols('e q'); import sympy as sp
print("  N(e,0) =", sp.factor(sp.expand(40*e**3-12*e**2-24*e-4)))
print("  N(e,1-e) =", sp.factor(sp.expand(N(e,1-e) if False else (40*e**3+40*e**2*(1-e)-12*e**2+83*e*(1-e)**2-32*e*(1-e)-24*e+2*(1-e)**3+(1-e)**2-8*(1-e)-4))))
