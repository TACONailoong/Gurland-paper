import mpmath as mp
mp.mp.dps = 50
print("--- decisive direct check: R = zeta(2,A-Delta) - Phi > 0 ---")
for (v,c) in [(9999,100),(9999,200),(9999,2000),(9999,8000),(99,80),(1e5,5e4),(1e6,9e5)]:
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3
    Phi=sum(c**(2*k-2)/k*mp.zeta(2*k,A) for k in range(1,700))
    R=mp.zeta(2,A-D)-Phi
    tag = "POSITIVE" if R>0 else "NEGATIVE"
    print(f"   v={mp.nstr(v,7)} c={mp.nstr(c,7)} rho={mp.nstr(c/A,5)}: R={mp.nstr(R,6)}  {tag}")
