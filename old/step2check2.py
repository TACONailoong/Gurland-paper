
from fractions import Fraction as F
c = {1: F(1,2)}
for j in range(2, 70): c[j] = c[j-1] * F(2*j-3, 2*j)
a = {0: F(0), 1: F(1,6)}
for m in range(2, 70): a[m] = c[m] - a[m-1]/3
b = lambda m: F(1, (m+1)*(2*m+1))
d = lambda m: F(3*(m-2), 2*(2*m-3)) * c[m]
Q = lambda m: d(m)/b(m)
ratio = lambda m: F((m-1)*(m+2)*(2*m-3)*(2*m+3), 2*(m-2)*(m+1)**2*(2*m+1))
print("Q(4) =", Q(4), "= 135/128?", Q(4) == F(135,128))
print("Q(m+1)=Q(m)*ratio(m) for m>=4:", all(Q(m+1) == Q(m)*ratio(m) for m in range(4,60)))
# n-form (m = n+4)
ratn = lambda n: F((n+3)*(n+6)*(2*n+5)*(2*n+11), 2*(n+2)*(n+5)**2*(2*n+9))
print("n-form == m-form:", all(ratn(n) == ratio(n+4) for n in range(0,40)))
print("ratn(n)-1 numerator identity:", all(ratn(n) - 1 == F(2*n**3+19*n**2+61*n+90, 2*(n+2)*(n+5)**2*(2*n+9)) for n in range(0,40)))
print("Q increasing:", all(Q(m) < Q(m+1) for m in range(4,60)))
print("d(m) < a(m) for m>=4:", all(d(m) < a(m) for m in range(4,60)))
print("b(m) < d(m) for m>=4:", all(b(m) < d(m) for m in range(4,60)))
print("b(2)<a(2), b(3)<a(3):", b(2) < a[2], b(3) < a[3], " a2-b2 =", a[2]-b(2), " a3-b3 =", a[3]-b(3))
