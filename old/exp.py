import mpmath as mp
mp.mp.dps = 40
def pieces(A, rho):
    A=mp.mpf(A); r=mp.mpf(rho); v=A-1; c=r*A
    s=mp.sqrt(v*v-c*c); Delta=(v-s)/3
    Phi = sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,500))
    R = mp.zeta(2,A-Delta)-Phi
    t = Delta/A
    D = t/(1-t) - (2*mp.atanh(r)/r + mp.log(1-r*r)/(r*r) - 1)
    S1 = -mp.log(1-r*r)/(r*r) - 1
    E1 = t*(2-t)/(2*(1-t)**2) - S1/2
    return R, D, E1
print("  A      rho     A*R            D              (A*R-D)*A       E1")
for Ai in ['1e3','1e4','1e5']:
    for r in ['0.2','0.5','0.9']:
        R,D,E1 = pieces(Ai,r); A=mp.mpf(Ai)
        print(f" {Ai:>6} {r:>5}  {mp.nstr(A*R,10):>14} {mp.nstr(D,10):>14} {mp.nstr((A*R-D)*A,10):>14} {mp.nstr(E1,10):>14}")
print()
print("--- small rho: is A*R - D  ~ E1/A ? ---")
for Ai in ['1e4','1e6']:
    for r in ['0.005','0.02','0.05']:
        R,D,E1 = pieces(Ai,r); A=mp.mpf(Ai)
        print(f" A={Ai} rho={r}: A*R-D = {mp.nstr(A*R-D,6)}   E1/A = {mp.nstr(E1/A,6)}   ratio = {mp.nstr((A*R-D)/(E1/A),6)}")
