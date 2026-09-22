import sympy as sp
A = sp.symbols('A', positive=True)
Psi = (A-1)*A**3/((A-sp.Rational(1,2))**3*(A+1))
diff = sp.simplify((1-sp.Rational(1,3)/A) - Psi)
print("(1-1/(3A)) - Psi  =", sp.factor(sp.together(diff)))
num = sp.factor(sp.together(diff).as_numer_denom()[0])
print("numerator       =", sp.expand(num))
print("roots           =", sp.nroots(num))
# ratio inequality for zeta4/zeta3
print()
z3 = sp.Rational(1,2)/A**2 + sp.Rational(1,2)/A**3        # lower bound for zeta(3,A)
z4 = 1/(3*(A-sp.Rational(1,2))**3)                        # upper bound for zeta(4,A)
d = sp.simplify(sp.Rational(2)*A**3/(3*(A-sp.Rational(1,2))**3*(A+1)) - z4/z3)
print("(2A^3/(3(A-1/2)^3(A+1))) - z4bound/z3bound =", sp.factor(sp.together(d)))
# final condition, as a polynomial
print()
rho2 = sp.symbols('rho2', positive=True)
cond = sp.simplify(1 - Psi*(1+rho2)/(1-rho2))
print("1 - Psi(1+rho2)/(1-rho2) with rho2 = 1/(3A):",
      sp.simplify(cond.subs(rho2, sp.Rational(1,3)/A)))
