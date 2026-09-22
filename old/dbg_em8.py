import mpmath as mp
mp.mp.dps = 45
def intK(T,c):
    return (-T*mp.log(T**2/(T**2-c**2)) + c*mp.log((T+c)/(T-c)))/c**2
x,y = mp.mpf(1000), mp.mpf(999); v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2
K = lambda t: mp.log(1/(1-c**2/t**2))/c**2
L = lambda t: 1/(t-w)**2
G = lambda t: K(t)-L(t)
S200 = sum(G(m+v) for m in range(1,201))
print("S200 =", mp.nstr(S200,22))
T0 = 200+v
I = intK(T0,c) - 1/(T0-w)
print("int_T0 G =", mp.nstr(I,22), "  G(T0)/2 =", mp.nstr(G(T0)/2,10))
for r in (1,2,3,4,5):
    Kd = mp.diff(K, T0, 2*r-1) if r<=3 else None
    pass
# exact tail pieces
for j in [1,3,5,7]:
    print(f"   K^({j})(T0) = {mp.nstr(mp.diff(K,T0,j),8)}   L^({j})(T0)*? = {mp.nstr((-1)**j*mp.factorial(j+1)*(T0-w)**(-(j+2)),8)}")
# brute-force tail: sum over huge range using asymptotic? instead compute partial tails
for N in [200,1000,5000,20000]:
    Tn = N+v
    print(f"   N={N:6d}: partial sum+int = {mp.nstr(sum(G(m+v) for m in range(1,N+1)) + (intK(Tn,c)-1/(Tn-w)),22)}")
