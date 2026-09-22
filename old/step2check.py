
from fractions import Fraction as F
from math import comb

c = {1: F(1,2)}
for j in range(2, 60):
    c[j] = c[j-1] * F(2*j-3, 2*j)
a = {0: F(0), 1: F(1,6)}
for m in range(2, 60):
    a[m] = c[m] - a[m-1]/3
b = lambda m: F(1, (m+1)*(2*m+1))
r = lambda m: a[m]/c[m]

print("c1..c5:", [c[i] for i in range(1,6)])
print("a2-b2 =", a[2]-b(2), " expect 1/360 =", F(1,360))
print("a3-b3 =", a[3]-b(3), " expect 11/3024 =", F(11,3024))
print("r2,r3 =", r(2), r(3), " expect 5/9, 17/27")
ok_bounds = all(F(3*(m-1),4*m) <= r(m) < F(3,4) for m in range(2,60))
print("bounds 3(m-1)/(4m) <= r_m < 3/4 for m>=2:", ok_bounds)
print("binom c_m matches:", all(c[m] == F(comb(2*m,m), 4**m * (2*m-1)) for m in range(1,25)))
Q = lambda m: F(3*(m-2), 2*(2*m-3)) * c[m] / b(m)
print("Q4 =", Q(4), " expect 135/128 =", F(135,128))
ratio = lambda m: (F(m-1)*F(m+2)*F(2*m-3)*F(2*m+3)) / (2*F(m-2)*F(m+1)**2*F(2*m+1))
print("Q_{m+1}/Q_m identity ok:", all(Q(m+1)/Q(m) == ratio(m) for m in range(4,50)))
print("Q increasing for m>=4:", all(Q(m) < Q(m+1) for m in range(4,50)))
print("poly 2m^3-5m^2+5m+22 > 0 for m>=4:", all(2*m**3-5*m**2+5*m+22 > 0 for m in range(4,50)))
print("a_m > b_m for 2<=m<60:", all(a[m] > b(m) for m in range(2,60)))
print("a_m > 3(m-2)/(2(2m-3)) c_m for m>=4:", all(a[m] > F(3*(m-2),2*(2*m-3))*c[m] for m in range(4,60)))
# final inequality check: 6S'(y) < 6h'(y)
import mpmath as mp
mp.mp.dps = 30
S = lambda y: mp.nsum(lambda m: y**m/((m+1)*(2*m+1)), [0, mp.inf])
Sp = lambda y: mp.nsum(lambda m: m*y**(m-1)/((m+1)*(2*m+1)), [1, mp.inf])
hp = lambda y: 3/(2*mp.sqrt(1-y)*(2+mp.sqrt(1-y))**2)
print("max of 6S'/6h' over y in (0,1):",
      max(float(6*Sp(mp.mpf(k)/100)/ (6*hp(mp.mpf(k)/100))) for k in list(range(1,100))+[990,995,999]))
