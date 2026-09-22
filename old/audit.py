import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 50
print("AUDIT 1 : h_k(A) = A^{2k-1}/k * zeta(2k,A)  <=  1/(kA) + 1/(k(2k-1))")
ok=True
for A in [mp.mpf('0.5'), mp.mpf(1), mp.mpf(3), mp.mpf(50), mp.mpf(1000)]:
    for k in [1,2,3,10,50,200]:
        h = A**(2*k-1)/k*mp.zeta(2*k,A)
        b = 1/(k*A)+1/(k*(2*k-1))
        ok &= (h<=b)
        if k in (1,3,50): print(f"   A={mp.nstr(A,6):>8} k={k:>4}: h_k={mp.nstr(h,10):>16} bound={mp.nstr(b,10):>16} {'ok' if h<=b else 'FAIL'}")
print("   lemma 1 holds everywhere:", ok)
print()
print("AUDIT 2 : h_k(A) -> 1/(k(2k-1)) as A -> infinity")
for k in [1,2,5]:
    for A in [10,100,1000,10000]:
        h = A**(2*k-1)/k*mp.zeta(2*k,A)
        tgt = mp.mpf(1)/(k*(2*k-1))
        print(f"   k={k} A={A:>6}: h_k={mp.nstr(h,12):>16} target={mp.nstr(tgt,12):>16} diff={mp.nstr(h-tgt,4)}")
print()
print("AUDIT 3 : sum_k rho^{2k-2}/(k(2k-1)) = Lambda(rho)/rho^2 ")
for r in ['0.1','0.5','0.9','0.99']:
    r=mp.mpf(r)
    s=mp.nsum(lambda k: r**(2*k-2)/(k*(2*k-1)), [1, mp.inf])
    Lam=(1+r)*mp.log(1+r)+(1-r)*mp.log(1-r)
    print(f"   rho={r}: series={mp.nstr(s,14):>18} Lambda/rho^2={mp.nstr(Lam/r**2,14):>18}")
print()
print("AUDIT 4 : Taylor expansion of zeta(2,A-delta) has positive coefficients")
A=mp.mpf(3); d=mp.mpf('0.7')
ser=sum((j+1)*d**j*mp.zeta(j+2,A) for j in range(0,200))
print("   zeta(2,A-delta) =", mp.nstr(mp.zeta(2, A-d), 20), "  Taylor series =", mp.nstr(ser, 20))
print()
print("AUDIT 5 : v-s = c^2/(v+s)  and  s<v")
for (x,y) in [(1,2),(0.1,199.9),(20,100)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; s=mp.sqrt(x*y); c=abs(x-y)/2
    print(f"   ({x},{y}): v-s={mp.nstr(v-s,12)}  c^2/(v+s)={mp.nstr(c**2/(v+s),12)}  s<v: {s<v}")
