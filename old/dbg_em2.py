import mpmath as mp, sympy as sp
mp.mp.dps = 35
_t, _cc = sp.symbols('t c', positive=True)
_Kexpr = -sp.log(1 - _cc**2/_t**2)/_cc**2
_ORD=[0,1,3,4,5,6,8,7,10]
_Kfun={j: sp.lambdify((_t,_cc), sp.diff(_Kexpr,_t,j), modules='mpmath') for j in _ORD}
B2K={1: mp.mpf(1)/6, 2: mp.mpf(-1)/30, 3: mp.mpf(1)/42}
def Kder(j,t,c): return _Kfun[j](mp.mpf(t), mp.mpf(c))
def Lder(j,t,w): return (-1)**j*mp.factorial(j+1)*mp.mpf(t-w)**(-(j+2))
def Gder(j,t,c,w): return Kder(j,t,c)-Lder(j,t,w)

for (xx,yy) in [(100,1),(1000,999),(1,2),(0.01,10000)]:
    x,y=mp.mpf(xx),mp.mpf(yy); v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2
    G=lambda t: Gder(0,t,c,w)
    I=mp.quad(G,[A,mp.inf]); g0=G(A)
    N=400
    true=sum(G(m+v) for m in range(1,N+1))+mp.quad(G,[N+v,mp.inf])
    print("="*70)
    print(f"(x,y)=({xx},{yy})  A={mp.nstr(A,10)} c={mp.nstr(c,8)} w={mp.nstr(w,8)}")
    print(f"   I={mp.nstr(I,12)}  G(A)/2={mp.nstr(g0/2,12)}  true={mp.nstr(true,12)}")
    for r in (1,2,3):
        E=I+g0/2
        for k in range(1,r+1):
            E += -B2K[k]/mp.factorial(2*k)*Gder(2*k-1,A,c,w)
        mm=2*r+2
        integ=mp.quad(lambda t: abs(Gder(mm,t,c,w)),[A,mp.inf])
        Rem=2*mp.zeta(mm)/(2*mp.pi)**mm*integ
        # size of the first omitted EM term
        nxt=abs(-B2K.get(r+1, mp.mpf(0))/mp.factorial(2*(r+1))*Gder(2*(r+1)-1,A,c,w)) if r<3 else mp.mpf(0)
        print(f"   r={r}: E={mp.nstr(E,12)}  |E-true|={mp.nstr(abs(E-true),3)}  Rem={mp.nstr(Rem,3)}"
              f"  int|G^({mm})|={mp.nstr(integ,3)}  valid={abs(E-true)<=Rem}")
