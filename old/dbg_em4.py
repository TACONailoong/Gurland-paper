import mpmath as mp, sympy as sp
mp.mp.dps = 40
t, cc = sp.symbols('t c', positive=True)
Kexpr = -sp.log(1 - cc**2/t**2)/cc**2
Lexpr = 1/(t-cc)**2   # placeholder, real L uses w not c
for (xx,yy) in [(100,1),(1000,999)]:
    x,y = mp.mpf(xx), mp.mpf(yy); v=(x+y)/2; c=abs(x-y)/2; A=1+v
    w = (v-mp.sqrt(x*y))/2
    Kf = sp.lambdify((t,cc), Kexpr, modules='mpmath')
    Gf = lambda z: Kf(z,c) - 1/(z-w)**2
    # exact derivatives of G at A
    terms=[]
    for k in range(1,13):
        dK = sp.lambdify((t,cc), sp.diff(Kexpr,t,2*k-1), modules='mpmath')
        val = dK(A,c) - (-1)**(2*k-1)*mp.factorial(2*k)*(A-w)**(-(2*k+1))
        B = mp.bernoulli(2*k)
        terms.append(-B/mp.factorial(2*k)*val)
    print(f"--- (x,y)=({xx},{yy})  A={mp.nstr(A,8)} c={mp.nstr(c,8)} w={mp.nstr(w,8)}")
    for k,T in enumerate(terms,1):
        print(f"   T_{k:2d} = {mp.nstr(T,6)}")
    # partial sums using the closed-form integral
    intK = lambda T_: (-T_*mp.log(T_**2/(T_**2-c**2)) + c*mp.log((T_+c)/(T_-c)))/c**2
    S0 = intK(A) - 1/(A-w) + Gf(A)/2
    acc = S0
    for k,T in enumerate(terms,1):
        acc += T
        print(f"   E_{k} = {mp.nstr(acc,12)}   (cumulative)")
