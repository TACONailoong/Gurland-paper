import mpmath as mp
mp.mp.dps = 35
print("--- identity:  Phi = (1/c^2) int_0^c [psi(A+t)-psi(A-t)] dt ---")
for (x,y) in [(1,2),(20,100),(1000,999),(0.1,199.9)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; c=abs(x-y)/2; A=1+v
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 300))
    I = mp.quad(lambda t: mp.digamma(A+t)-mp.digamma(A-t), [0, c])/c**2
    print(f"  ({x},{y}): Phi={mp.nstr(Phi,18)}  integral={mp.nstr(I,18)}  diff={mp.nstr(I-Phi,3)}")
print()
print("--- triangle average:  Phi = int_-c^c psi'(A+z) (c-|z|)/c^2 dz ---")
for (x,y) in [(1,2),(20,100)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; c=abs(x-y)/2; A=1+v
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k, A) for k in range(1, 300))
    I = mp.quad(lambda z: mp.zeta(2, A+z)*(c-abs(z)), [-c, c])/c**2
    print(f"  ({x},{y}): Phi={mp.nstr(Phi,18)}  tri={mp.nstr(I,18)}  diff={mp.nstr(I-Phi,3)}")
