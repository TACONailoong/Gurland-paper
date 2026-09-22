import mpmath as mp
mp.mp.dps = 45
def intK(T,c):
    return (-T*mp.log(T**2/(T**2-c**2)) + c*mp.log((T+c)/(T-c)))/c**2
for (xx,yy) in [(1000,999),(100,1),(0.01,10000)]:
    x,y=mp.mpf(xx),mp.mpf(yy); v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v
    w=(v-s)/2
    K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
    L=lambda t: 1/(t-w)**2
    G=lambda t: K(t)-L(t)
    print(f"### (x,y)=({xx},{yy})  v={v} c={c} w={mp.nstr(w,10)}")
    for N in [400,4000,40000]:
        S = sum(G(m+v) for m in range(1,N+1)) + (intK(N+v,c) - 1/(N+v-w))
        print(f"   direct sum N={N:6d}: {mp.nstr(S,20)}")
    # independent: Phi = sum K  (direct + tail), and exact zeta for sum L
    for N in [4000,40000]:
        Phi = sum(K(m+v) for m in range(1,N+1)) + intK(N+v,c)
        zl  = mp.zeta(2, 1+v-w)
        print(f"   route2  N={N:6d}: Phi={mp.nstr(Phi,20)} zetaL={mp.nstr(zl,20)} diff={mp.nstr(Phi-zl,20)}")
