import mpmath as mp, sympy as sp
mp.mp.dps = 35
def intK(T, c):   # closed form of int_T^inf K
    return (-T*mp.log(T**2/(T**2-c**2)) + c*mp.log((T+c)/(T-c)))/c**2
for (xx,yy) in [(100,1),(1000,999),(0.01,10000),(1,2),(2,7)]:
    x,y=mp.mpf(xx),mp.mpf(yy); v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2
    K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
    L=lambda t: 1/(t-w)**2
    G=lambda t: K(t)-L(t)
    q  = mp.quad(G,[A,mp.inf])
    ex = intK(A,c) - 1/(A-w)
    print(f"(x,y)=({xx},{yy}):  quad={mp.nstr(q,14)}  closed={mp.nstr(ex,14)}  diff={mp.nstr(q-ex,4)}  rel={mp.nstr(abs((q-ex)/ex),4)}")
    # tail comparison
    T=A+400
    qt = mp.quad(G,[T,mp.inf]); ext = intK(T,c) - 1/(T-w)
    print(f"        tail: quad={mp.nstr(qt,12)}  closed={mp.nstr(ext,12)}  diff={mp.nstr(qt-ext,4)}")
