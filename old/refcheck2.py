import mpmath as mp
from mpmath import binomial
mp.mp.dps=25
cm=lambda m: mp.mpf(binomial(2*m,m))/(4**m*(2*m-1))
print("(2) c_m - c_{m-1}/4  vs  3(m-2)c_m/(2(2m-3)):")
for m in [4,5,10,50]:
    lhs=cm(m)-cm(m-1)/4; rhs=3*(m-2)*cm(m)/(2*(2*m-3))
    print(f"    m={m}: {mp.nstr(lhs,12)} vs {mp.nstr(rhs,12)}  equal: {abs(lhs-rhs)<mp.mpf(10)**-20}")
print("(3) -log(1-y) -(y+y^2/2+y^3/(3(1-y))) at y=0.5:", mp.nstr(-mp.log(mp.mpf('0.5'))-(mp.mpf('0.5')+mp.mpf('0.125')+mp.mpf('0.125')/(3*mp.mpf('0.5'))),8), "(<0 means the '<=' form is right)")
print("(5) Psi(0.1) =", mp.nstr((mp.mpf('0.1')-1)*mp.mpf('0.1')**3/((mp.mpf('0.1')-mp.mpf('0.5'))**3*mp.mpf('1.1')),6), " vs 1-1/(3*0.1) =", mp.nstr(1-1/(3*mp.mpf('0.1')),6))
print("    Psi(1) =", mp.nstr(0,3), " 1-1/3 =", mp.nstr(1-mp.mpf(1)/3,6))
