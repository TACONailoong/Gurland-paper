import mpmath as mp
mp.mp.dps = 40
def explicitK(T,c):
    return (-T*mp.log(T**2/(T**2-c**2)) + c*mp.log((T+c)/(T-c)))/c**2
B2K = {1: mp.mpf(1)/6, 2: mp.mpf(-1)/30, 3: mp.mpf(1)/42}
def K(t,c): return mp.log(1/(1-c**2/t**2))/c**2
def Kder(j,t,c):
    if j==0: return K(t,c)
    return mp.diff(lambda z: K(z,c), t, j)
for (xx,yy) in [(1000,999),(100,1)]:
    x,y=mp.mpf(xx),mp.mpf(yy); v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2
    L=lambda t: 1/(t-w)**2
    G=lambda t: K(t,c)-L(t)
    # certificate
    I = explicitK(A,c) - 1/(A-w)
    E={}
    for r in (1,2,3):
        e = I + G(A)/2
        for k in range(1,r+1):
            e += -B2K[k]/mp.factorial(2*k)*(Kder(2*k-1,A,c) - (-1)**(2*k-1)*mp.factorial(2*k)*(A-w)**(-(2*k+1)))
        E[r]=e
    # reference route 2
    print(f"### ({xx},{yy})  E_3 = {mp.nstr(E[3],22)}")
    for N in [200, 2000, 20000]:
        T0=N+v
        Phi = sum(K(m+v,c) for m in range(1,N+1)) + explicitK(T0,c) - K(T0,c)/2
        for k in range(1,4):
            Phi += -B2K[k]/mp.factorial(2*k)*Kder(2*k-1,T0,c)
        true = Phi - mp.zeta(2, 1+(v+s)/2)
        print(f"    N={N:6d}: Phi={mp.nstr(Phi,22)}  true={mp.nstr(true,22)}  true-E_3={mp.nstr(true-E[3],6)}")
    print(f"    zeta(2,1+(v+s)/2) = {mp.nstr(mp.zeta(2,1+(v+s)/2),22)}")
