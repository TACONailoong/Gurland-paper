
import sympy as sp
e,a,b,t = sp.symbols('e a b t')
N = 40*e**3 + 40*e**2*sp.Symbol('q') - 12*e**2 + 83*e*sp.Symbol('q')**2 - 32*e*sp.Symbol('q') - 24*e + 2*sp.Symbol('q')**3 + sp.Symbol('q')**2 - 8*sp.Symbol('q') - 4
q = sp.Symbol('q')
Nf = sp.lambdify((e,q), N)
def Nv(ee,qq): return Nf(ee,qq)
lhs = (1-t)*Nv(e,a) + t*Nv(e,b) - Nv(e, (1-t)*a + t*b)
rhs = t*(1-t)*(b-a)**2*((83*e+1) + 2*(3*a + (b-a)*(1+t)))
print("identity holds:", sp.simplify(sp.expand(lhs - rhs)) == 0)
