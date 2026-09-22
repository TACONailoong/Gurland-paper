import mpmath as mp, sympy as sp
mp.mp.dps = 35
t, c = sp.symbols('t c', positive=True)
Kexpr = -sp.log(1 - c**2/t**2)/c**2
Iexpr = sp.integrate(Kexpr, (t, t, sp.oo))
print("sympy integral =", sp.simplify(Iexpr))
print("   as string:", sp.sstr(sp.simplify(Iexpr)))
Ki = sp.lambdify((t,c), Iexpr, modules='mpmath')
def explicit(T, cc):
    return (-T*mp.log(T**2/(T**2-cc**2)) + cc*mp.log((T+cc)/(T-cc)))/cc**2
K = lambda z,cc: mp.log(1/(1-cc**2/z**2))/cc**2
for (T,cc) in [(1000.5,0.5),(1400.5,0.5),(51.5,49.5),(3.5,0.5),(5200.0,4999.995)]:
    a_m = Ki(mp.mpf(T), mp.mpf(cc)); a_f = Ki(T, cc); e = explicit(mp.mpf(T), mp.mpf(cc))
    q = mp.quad(lambda z: K(z,mp.mpf(cc)), [mp.mpf(T), mp.inf])
    print(f"T={T} c={cc}: mpf-args={mp.nstr(a_m,18)} float-args={mp.nstr(a_f,18)} explicit={mp.nstr(e,18)} quad={mp.nstr(q,18)}")
