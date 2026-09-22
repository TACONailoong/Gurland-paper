import sympy as sp, mpmath as mp
m=sp.symbols('m', positive=True, integer=True)
A,t,y=sp.symbols('A t y', positive=True)
cm = sp.binomial(2*m,m)/(4**m*(2*m-1))
print("(2) c_m - c_{m-1}/4  vs  3(m-2)c_m/(2(2m-3)):",
      sp.simplify(sp.expand(1-sp.Rational(1,4)*sp.Rational(2*m,2*m-3) - sp.Rational(3,2)*(m-2)/(2*m-3))))
print("(3) -log(1-y) - (y + y^2/2 + y^3/(3(1-y))) at y=0.5:",
      mp.nstr(-mp.log(0.5)-(0.5+0.125+0.125/(3*0.5)),8), " (should be < 0)")
print("(5) 4A^4-14A^3+21A^2-8A+1 with A=t+1:",
      sp.expand(4*(t+1)**4-14*(t+1)**3+21*(t+1)**2-8*(t+1)+1))
print("    Psi(0.1) vs 1-1/(3*0.1):",
      mp.nstr((mp.mpf('0.1')-1)*mp.mpf('0.1')**3/((mp.mpf('0.1')-mp.mpf('0.5'))**3*mp.mpf('1.1')),6),
      "vs", mp.nstr(1-1/(3*mp.mpf('0.1')),6))
