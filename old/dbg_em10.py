import mpmath as mp, sympy as _sp
mp.mp.dps = 35
B2K = {1: mp.mpf(1)/6, 2: mp.mpf(-1)/30, 3: mp.mpf(1)/42}
_t, _cc = _sp.symbols('t c', positive=True)
_Kexpr = -_sp.log(1 - _cc**2/_t**2)/_cc**2
_ORD = [0,1,3,4,5,6,8]
_Kfun = {j: _sp.lambdify((_t,_cc), _sp.diff(_Kexpr,_t,j), modules='mpmath') for j in _ORD}
_Kint = _sp.lambdify((_t,_cc), _sp.integrate(_Kexpr,(_t,_t,_sp.oo)), modules='mpmath')
def Kder(j,t,c): return _Kfun[j](mp.mpf(t), mp.mpf(c))
def Lder(j,t,w): return (-1)**j*mp.factorial(j+1)*mp.mpf(t-w)**(-(j+2))
def Gder(j,t,c,w): return Kder(j,t,c)-Lder(j,t,w)

def em_data(x, y, orders=(1,2,3), with_true=False):
    v = (x + y) / 2; c = abs(x - y) / 2; s = mp.sqrt(x * y); A = 1 + v
    w = (v - s) / 2
    G = lambda t: Gder(0, t, c, w)
    I = mp.quad(G, [A, mp.inf]); g0 = G(A)
    out = {}
    for r in orders:
        E = I + g0 / 2
        for k in range(1, r + 1):
            E += -B2K[k] / mp.factorial(2 * k) * Gder(2 * k - 1, A, c, w)
        m = 2 * r + 2
        Rem = (2 * mp.zeta(m) / (2 * mp.pi) ** m
               * mp.quad(lambda t: abs(Gder(m, t, c, w)), [A, mp.inf]))
        out[r] = (E, Rem)
    true = None
    if with_true:
        N = max(200, 5 * int(c) + 1)
        T0 = N + v
        Phi = sum(Kder(0, m + v, c) for m in range(1, N + 1))
        Phi += _Kint(T0, c) - Kder(0, T0, c) / 2
        for k in range(1, 4):
            Phi += -B2K[k] / mp.factorial(2 * k) * Kder(2 * k - 1, T0, c)
        true = Phi - mp.zeta(2, 1 + (v + s) / 2)
    return out, true, I, g0, A, c, w, v

for (x,y) in [(1000,999)]:
    out, true, I, g0, A, c, w, v = em_data(x,y,with_true=True)
    print(f"types: v={type(v).__name__} c={type(c).__name__} A={type(A).__name__}")
    print("  A =", mp.nstr(mp.mpf(A),20), " c =", mp.nstr(mp.mpf(c),20), " w =", mp.nstr(w,20))
    print("  I(quad) =", mp.nstr(I,22), "  I(closed) =", mp.nstr(_Kint(A,c)-1/(A-w),22))
    print("  g0/2 =", mp.nstr(g0/2,12))
    print("  E_3 =", mp.nstr(out[3][0],22), " Rem =", mp.nstr(out[3][1],4))
    print("  true =", mp.nstr(true,22), "   true-E_3 =", mp.nstr(true-out[3][0],4))
