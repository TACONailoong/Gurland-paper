import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 45

def parts(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    u=(x-y)/2; v=(x+y)/2; c=abs(u); s=mp.sqrt(x*y); A=1+v
    return u,v,c,s,A

def lnGs(x,y):
    u,v,c,s,A = parts(x,y)
    return mp.loggamma(A+c)+mp.loggamma(A-c)-2*mp.loggamma(A)

# ---- 1) confirm the algebra slip: T' bound must carry A^{-2}
print("="*74)
print("1) T' bound check:  T' <= zeta(4,A) c^4 /(3 A^2 (1-rho^2))   [correct]  vs")
print("                   T' <= zeta(4,A) c^4 /(3 (1-rho^2))       [as written]")
print("="*74)
for (x,y) in [(1,10),(0.1,199.9),(0.5,3),(20,100)]:
    u,v,c,s,A = parts(x,y)
    Phi = lnGs(x,y)/c**2
    Tp  = Phi - mp.zeta(2,A) - c**2/2*mp.zeta(4,A)
    rho = c/A
    b_ok  = mp.zeta(4,A)*c**4/(3*A**2*(1-rho**2))
    b_bad = mp.zeta(4,A)*c**4/(3*(1-rho**2))
    print(f"  (x,y)=({x},{y}) rho={mp.nstr(rho,6)}")
    print(f"      T'={mp.nstr(Tp,8):>14}   correct bound={mp.nstr(b_ok,8):>14} ok={Tp<=b_ok}"
          f"   as-written bound={mp.nstr(b_bad,8):>14} ok={Tp<=b_bad}")

# ---- 2) Euler-Maclaurin viability for the unconditional midpoint lemma
print()
print("="*74)
print("2) Euler-Maclaurin test for  sum_{n>=1} G(n+v) < 0,  G=K-L")
print("   K(x)=c^-2 ln(1/(1-c^2/x^2)),  L(x)=1/(x-w)^2,  w=(v+s)/2")
print("="*74)
def K(x,c): return mp.log(1/(1-c**2/x**2))/c**2
def L(x,w): return 1/(x-w)**2
def EM_terms(x,y,order=6):
    u,v,c,s,A = parts(x,y); w=(v+s)/2
    G  = lambda t: K(t,c)-L(t,w)
    I  = mp.quad(lambda t: K(t,c)-L(t,w), [A, mp.inf])
    tot = I + G(A)/2
    terms=[I, G(A)/2]
    for j in [1,3,5]:
        d = mp.diff(G, A, j)
        coeff = {1: mp.mpf(-1)/12, 3: mp.mpf(1)/720, 5: mp.mpf(-1)/30240}[j]
        tot += coeff*d
        terms.append(coeff*d)
    return tot, terms
worst=[]; okcnt=0; tot=0
for xv in [0.001,0.01,0.1,0.5,1,2,5,20,100,1000]:
    for yv in [0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000]:
        if xv==yv: continue
        try:
            t,terms = EM_terms(xv,yv)
        except Exception as e:
            print("  err",xv,yv,repr(e)[:60]); continue
        true = sum(K(n+ (mp.mpf(xv)+mp.mpf(yv))/2, abs(mp.mpf(xv)-mp.mpf(yv))/2) - L(n+(mp.mpf(xv)+mp.mpf(yv))/2, (mp.mpf(xv)+mp.mpf(yv))/2+mp.sqrt(mp.mpf(xv)*mp.mpf(yv)))/2*0 - L(n+(mp.mpf(xv)+mp.mpf(yv))/2, ((mp.mpf(xv)+mp.mpf(yv))/2+mp.sqrt(mp.mpf(xv)*mp.mpf(yv)))/2) for n in range(1,6000))
        tot+=1
        if t<0: okcnt+=1
        worst.append((float(t), xv, yv, float(true)))
worst.sort(reverse=True)
print("  EM partial sum < 0 in", okcnt, "/", tot, "cases")
print("  largest (worst) EM partial sums:")
for w in worst[:8]:
    print(f"     EM={w[0]:+.6e}  true sum={w[3]:+.6e}   (x,y)=({w[1]},{w[2]})")
