import sys
sys.path.insert(0, r"D:\AI agent workplace\deepseek agent\projects\gurland-paper")
import mpmath as mp
mp.mp.dps = 35
import sympy as sp
t, c = sp.symbols('t c', positive=True)
Kexpr = -sp.log(1 - c**2/t**2)/c**2
Ki = sp.lambdify((t, c), sp.integrate(Kexpr, (t, t, sp.oo)), modules='mpmath')
def explicit(T, cc):
    return (-T*mp.log(T**2/(T**2-cc**2)) + cc*mp.log((T+cc)/(T-cc)))/cc**2
for (T, cc) in [(51.5, 49.5), (1400.5, 0.5), (5200.0, 4999.995), (3.5, 0.5)]:
    a = Ki(mp.mpf(T), mp.mpf(cc)); b = explicit(mp.mpf(T), mp.mpf(cc))
    print(f"T={T} c={cc}: sympy={mp.nstr(a,18)} explicit={mp.nstr(b,18)} diff={mp.nstr(a-b,3)}")
