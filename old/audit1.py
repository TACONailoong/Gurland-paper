import mpmath as mp
mp.mp.dps = 40
ok = lambda b: "OK " if b else "**FAIL**"

# ---- Lemma (VI): F(y) < 9/(sqrt(1-y)(2+sqrt(1-y))^2) ----
def F(y):  return 6*sum(mp.mpf(m)/((m+1)*(2*m+1))*y**(m-1) for m in range(1,4000))
def G(y):
    q=mp.sqrt(1-y); return 9/(q*(2+q)**2)
print("--- Lemma (VI):  F(y) < G(y) on (0,1) ---")
worst=(mp.mpf(-9),None)
for j in range(1,4000):
    y=mp.mpf(j)/4000
    r=F(y)/G(y)
    if r>worst[0]: worst=(r,y)
print("   max F/G =", mp.nstr(worst[0],12), "at y =", mp.nstr(worst[1],5), ok(worst[0]<1))

# ---- coefficient check (VII) ----
print()
print("--- (VII) coefficients: a1=b1, a_m>b_m (m>=2) ---")
from mpmath import binomial
c=lambda m: mp.mpf(binomial(2*m,m))/(4**m*(2*m-1))
a=[None,mp.mpf(1)/6]
for m in range(2,30): a.append(c(m)-a[m-1]/3)
b=lambda m: mp.mpf(1)/((m+1)*(2*m+1))
bad=[m for m in range(2,30) if not a[m]>b[m]]
print("   a1-b1 =", mp.nstr(a[1]-b(1),5), " a2-b2 =", mp.nstr(a[2]-b(2),6),
      " a3-b3 =", mp.nstr(a[3]-b(3),6), " violations m>=2:", bad)

# ---- (IX) induction claim ----
print()
print("--- (IX): 3(m-1)/(4m) <= r_m < 3/4 ---")
r=[None,mp.mpf(1)/6/c(1)]
for m in range(2,60): r.append(a[m]/c(m))
bad2=[m for m in range(2,60) if not (mp.mpf(3)*(m-1)/(4*m) <= r[m] < mp.mpf(3)/4)]
print("   violations:", bad2, "  r_2..r_5 =", [mp.nstr(r[m],6) for m in range(2,6)])

# ---- (X): G(y) < L_e(y) ----
print()
print("--- (X):  G(y) < B^3(z+1)/(s z^3)  at many (v,c) ---")
def right(v,c):
    A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3; z=A-D
    B=A-mp.mpf(1)/2
    return B**3*(z+1)/(s*z**3)
worst2=(mp.mpf(9),None,None)
for vi in ['0.3','1','2','5','20','100','1000','10000']:
    v=mp.mpf(vi)
    for j in range(1,60):
        c=v*mp.mpf(j)/60
        A=1+v; B=A-mp.mpf(1)/2; y=(c/B)**2
        r=G(y)/right(v,c)
        if r>worst2[0]: worst2=(r,v,c)
print("   max G/right =", mp.nstr(worst2[0],12), "at v =", mp.nstr(worst2[1],6), ", c =", mp.nstr(worst2[2],6), ok(worst2[0]<1))
