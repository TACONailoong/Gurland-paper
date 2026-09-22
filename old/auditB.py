import mpmath as mp
mp.mp.dps = 40
from mpmath import binomial
c=lambda m: mp.mpf(binomial(2*m,m))/(4**m*(2*m-1))
bb=lambda m: mp.mpf(1)/((m+1)*(2*m+1))
a=[None,mp.mpf(1)/6]
for m in range(2,40): a.append(c(m)-a[m-1]/3)
print("--- (VII): a1=b1 and a_m>b_m for m>=2 ---")
print("   a1-b1 =", mp.nstr(a[1]-bb(1),5), "  a2-b2 =", mp.nstr(a[2]-bb(2),6),
      "  a3-b3 =", mp.nstr(a[3]-bb(3),6))
print("   violations:", [m for m in range(2,40) if not a[m]>bb(m)])
print()
print("--- (IX): 3(m-1)/(4m) <= r_m < 3/4 ---")
r=[None,mp.mpf(1)/6/c(1)]
for m in range(2,60): r.append(a[m]/c(m))
print("   violations:", [m for m in range(2,60) if not (mp.mpf(3)*(m-1)/(4*m) <= r[m] < mp.mpf(3)/4)])
print("   r_2..r_6 =", [mp.nstr(r[m],7) for m in range(2,7)])
print()
def Q(m): return (mp.mpf(3)*(m-2)*c(m)/(2*(2*m-3)))/bb(m)
print("--- Q_m>1 for m>=4 ? ---")
print("   Q_4 =", mp.nstr(Q(4),10), "  violations:", [m for m in range(4,80) if not Q(m)>1])
print("   Q_{m+1}/Q_m m=4..8:", [mp.nstr(Q(m+1)/Q(m),8) for m in range(4,9)])
print("   their formula at m=5 (should be Q_6/Q_5):",
      mp.nstr((mp.mpf(4)*7*7*13)/(2*3*36*11),8), " actual:", mp.nstr(Q(6)/Q(5),8))
