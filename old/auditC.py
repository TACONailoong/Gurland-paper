import mpmath as mp
mp.mp.dps = 40
from mpmath import binomial
c=lambda m: mp.mpf(binomial(2*m,m))/(4**m*(2*m-1))
bb=lambda m: mp.mpf(1)/((m+1)*(2*m+1))
N=80
a=[None,mp.mpf(1)/6]
for m in range(2,N+1): a.append(c(m)-a[m-1]/3)
r=[None,mp.mpf(1)/6/c(1)]
for m in range(2,N+1): r.append(a[m]/c(m))
print("--- (IX) violations:", [m for m in range(2,N+1) if not (mp.mpf(3)*(m-1)/(4*m) <= r[m] < mp.mpf(3)/4)][:5],
      "  r_2..r_6 =", [mp.nstr(r[m],7) for m in range(2,7)])
def Q(m): return (mp.mpf(3)*(m-2)*c(m)/(2*(2*m-3)))/bb(m)
print("--- Q_m>1 violations (m>=4):", [m for m in range(4,N) if not Q(m)>1][:5],
      "  Q_4=", mp.nstr(Q(4),10))
print("--- Q_{m+1}/Q_m m=4..7:", [mp.nstr(Q(m+1)/Q(m),8) for m in range(4,8)],
      " their formula m=5:", mp.nstr((mp.mpf(4)*7*7*13)/(2*3*36*11),8))

print()
print("=== X): G(y) < B^3(z+1)/(s z^3) ; and the FULL chain ===")
F=lambda y: 6*sum(mp.mpf(m)/((m+1)*(2*m+1))*y**(m-1) for m in range(1,3000))
G=lambda y: 9/(mp.sqrt(1-y)*(2+mp.sqrt(1-y))**2)
def E(v,c):
    A=1+v
    return sum(mp.mpf(m)/(m+1)*c**(2*m)*mp.zeta(2*m+2,A) for m in range(1,600))
def data(v,c):
    A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3; z=A-D; B=A-mp.mpf(1)/2; y=(c/B)**2
    return A,s,D,z,B,y
worstX=(mp.mpf(-9),None,None); worstMid=(mp.mpf(-9),None,None); worstFin=(mp.mpf(-9),None,None)
for vi in ['0.3','0.6','1','1.5','2','5','20','100','1000','10000']:
    v=mp.mpf(vi)
    for j in range(1,80):
        cc=v*mp.mpf(j)/80
        A,s,D,z,B,y=data(v,cc)
        rX=G(y)/(B**3*(z+1)/(s*z**3))
        rM=E(v,cc)/(cc**2/(6*B**3)*F(y))
        rZ=E(v,cc)/(cc**2*zeta3) if False else E(v,cc)/(cc**2/(3*s)*mp.zeta(3,z))
        if rX>worstX[0]: worstX=(rX,v,cc)
        if rM>worstMid[0]: worstMid=(rM,v,cc)
        if rZ>worstFin[0]: worstFin=(rZ,v,cc)
print("  (X)  max G/right    =", mp.nstr(worstX[0],10), "at v=",mp.nstr(worstX[1],5),"c=",mp.nstr(worstX[2],5),  "OK" if worstX[0]<1 else "**FAIL**")
print("  mid  max E/(c^2F/(6B^3)) =", mp.nstr(worstMid[0],10), "at v=",mp.nstr(worstMid[1],5),"c=",mp.nstr(worstMid[2],5), "OK" if worstMid[0]<1 else "**FAIL**")
print("  (I)  max E/(c^2 zeta3/(3s)) =", mp.nstr(worstFin[0],10), "at v=",mp.nstr(worstFin[1],5),"c=",mp.nstr(worstFin[2],5), "OK" if worstFin[0]<1 else "**FAIL**")
