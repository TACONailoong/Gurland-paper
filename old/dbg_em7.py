import mpmath as mp
mp.mp.dps = 45
def intK(T,c):
    return (-T*mp.log(T**2/(T**2-c**2)) + c*mp.log((T+c)/(T-c)))/c**2
x,y = mp.mpf(1000), mp.mpf(999); v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2
K = lambda t: mp.log(1/(1-c**2/t**2))/c**2
L = lambda t: 1/(t-w)**2
G = lambda t: K(t)-L(t)
Kd = {0:K}
for j in [1,3,5]:
    Kd[j] = mp.diff(K, A, j)
Ld = lambda j,t: (-1)**j*mp.factorial(j+1)*(t-w)**(-(j+2))
B = {1: mp.mpf(1)/6, 2: mp.mpf(-1)/30, 3: mp.mpf(1)/42}
I = intK(A,c) - 1/(A-w)
print("I (closed)      =", mp.nstr(I,20))
for r in (1,2,3):
    E = I + G(A)/2
    for k in range(1,r+1):
        E += -B[k]/mp.factorial(2*k)*(Kd[2*k-1]-Ld(2*k-1,A))
    print(f"  E_{r} = {mp.nstr(E,22)}")
print("  G(A) =", mp.nstr(G(A),10), "  G'(A)=", mp.nstr(Kd[1]-Ld(1,A),10), " G'''(A)=", mp.nstr(Kd[3]-Ld(3,A),10))
for N in [200,2000,20000,200000]:
    T0 = N+v
    direct = sum(G(m+v) for m in range(1,N+1)) + (intK(T0,c)-1/(T0-w))
    print(f"  N={N:7d}  sum+int = {mp.nstr(direct,22)}   +G(T0)/2 = {mp.nstr(direct+G(T0)/2,22)}")
