import mpmath as mp, sympy as sp
mp.mp.dps = 40
x, y = mp.mpf(100), mp.mpf(1)
v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2
print("A,c,w =", A, c, w)

t_, c_ = sp.symbols('t c', positive=True)
Kexpr = -sp.log(1 - c_**2/t_**2)/c_**2
f8 = sp.lambdify((t_, c_), sp.diff(Kexpr, t_, 8), modules='mpmath')
f0 = sp.lambdify((t_, c_), Kexpr, modules='mpmath')
def Ld(j,t): return (-1)**j * mp.factorial(j+1) * mp.mpf(t-w)**(-(j+2))
def Kd(j,t):
    if j==0: return f0(mp.mpf(t), c)
    return f8(mp.mpf(t), c) if j==8 else None
print("K(A)      =", mp.nstr(Kd(0,A),12))
print("K''''''''(A) =", mp.nstr(f8(mp.mpf(A), c),12))
print("L''''''''(A) =", mp.nstr(Ld(8,A),12))
G8 = lambda t: f8(mp.mpf(t), c) - Ld(8, t)
for dt in [0, mp.mpf('0.5'), 1, 2, 5, 10, 50, 200]:
    print(f"   G8(A+{dt}) = {mp.nstr(G8(A+dt),8)}")
print()
print("quad over [A, inf)      :", mp.nstr(mp.quad(lambda t: abs(G8(t)), [A, mp.inf]),12))
print("quad over [A, 2A, 10A] :", mp.nstr(mp.quad(lambda t: abs(G8(t)), [A, 2*A, 10*A, mp.inf]),12))
print("quad over [A, 10A]     :", mp.nstr(mp.quad(lambda t: abs(G8(t)), [A, 10*A]),12))
print("substitution t=A/u     :", mp.nstr(mp.quad(lambda u: abs(G8(A/u))*A/u**2, [0,1]),12))
print("quad with degree=8     :", mp.nstr(mp.quad(lambda t: abs(G8(t)), [A, 2*A, 20*A, mp.inf], maxdegree=8),12))
print()
print("crude: |G8(A)| * A/8   =", mp.nstr(abs(G8(A))*A/8,12))
