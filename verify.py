#!/usr/bin/env python3
"""
Numerical verification suite for
  "A sharp localization constant for Gurland's ratio of the gamma function".

Usage:  python verify.py            (all sections)
        python verify.py 3 5        (only sections 3 and 5)

Every numbered section checks the corresponding statement of the paper and reports
PASS/FAIL; the exit code is 0 iff every executed check passed.  Requires mpmath.
"""
import sys, time
import mpmath as mp

mp.mp.dps = 35
RESULTS = []

def check(tag, ok, detail=""):
    RESULTS.append((tag, bool(ok), detail))
    print(f"   [{'ok' if ok else 'FAIL'}] {tag}" + (f"   {detail}" if detail else ""))
    return bool(ok)

def Lam(r):                       # Lambda(rho), eq. (16) = eq:Lambda
    return (1 + r) * mp.log(1 + r) + (1 - r) * mp.log(1 - r)

def logGs(x, y):                  # log G_star
    v = (x + y) / 2
    return mp.loggamma(1 + x) + mp.loggamma(1 + y) - 2 * mp.loggamma(1 + v)

def logG(x, y):                   # log G
    v = (x + y) / 2
    return mp.loggamma(x) + mp.loggamma(y) - 2 * mp.loggamma(v)

def series(x, y, shift, kmax=3000):
    u = (x - y) / 2; v = (x + y) / 2
    s = mp.mpf(0)
    for k in range(1, kmax + 1):
        t = mp.mpf(u ** (2 * k)) / k * mp.zeta(2 * k, shift + v)
        s += t
        if abs(t) < abs(s) * mp.mpf(10) ** (-mp.mp.dps + 4) and k > 3:
            break
    return s

def tval(x, y):                   # the parameter t, eq. (10) = eq:Wt
    x = mp.mpf(x); y = mp.mpf(y)
    u = (x - y) / 2; v = (x + y) / 2; s = mp.sqrt(x * y)
    if u == 0:
        return None
    target = logGs(x, y) / u ** 2
    f = lambda tt: mp.zeta(2, 1 + tt) - target
    lo, hi = s, v
    if not (f(lo) > 0 > f(hi)):
        return None
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) > 0: lo = mid
        else: hi = mid
    return (lo + hi) / 2

GRID = [0.001, 0.01, 0.1, 0.5, 1, 2, 5, 20, 100, 1000, 10000]

# ---------------------------------------------------------------- section 1
def sec1():
    print("\n[1] Theorem 1.1 -- exact expansion, both normalisations")
    for (x, y) in [(1, 2), (1, 10), (0.5, 3), (0.2, 0.9), (3, mp.mpf('3.0001'))]:
        x = mp.mpf(x); y = mp.mpf(y)
        d1 = abs(logGs(x, y) - series(x, y, 1))
        d2 = abs(logG(x, y) - series(x, y, 0))
        check(f"expansion at ({mp.nstr(x,7)},{mp.nstr(y,9)})",
              d1 < mp.mpf(10) ** (-(mp.mp.dps - 6)) and d2 < mp.mpf(10) ** (-(mp.mp.dps - 6)),
              f"|dG*|={mp.nstr(d1,3)}, |dG|={mp.nstr(d2,3)}")
    check("optimality of the region: |x-y|<x+y fails at (1,-1/2)",
          not (abs(mp.mpf(1) + mp.mpf('0.5')) < mp.mpf(1) - mp.mpf('0.5')))
    check("optimality of the region: |x-y|<2+x+y fails at (3,-2)",
          not (abs(mp.mpf(3) + 2) < 2 + 3 - 2))

# ---------------------------------------------------------------- section 2
def sec2():
    print("\n[2] Corollary 1.2 -- the hypothesis Q<1 is redundant")
    x, y = mp.mpf(100), mp.mpf(1)
    Q = abs(x - y) / (2 * (1 + mp.sqrt(x * y)))
    u = (x - y) / 2; v = (x + y) / 2
    tails = {}
    for m in (20, 100, 1000):
        tails[m] = sum(mp.mpf(u ** (2 * k)) / k * mp.zeta(2 * k, 1 + v)
                       for k in range(m, m + 2500))
    check("Q(100,1)>1 yet the series still converges", Q > 1 and tails[1000] < mp.mpf(10) ** -30,
          f"Q={mp.nstr(Q,4)}, R_20={mp.nstr(tails[20],4)}, R_1000={mp.nstr(tails[1000],4)}")

# ---------------------------------------------------------------- section 3
def sec3():
    print("\n[3] Theorem 1.3 -- two-sided remainder estimate; the rate rho is optimal")
    for (x, y) in [(1, 10), (0.5, 3), (2, 7), (0.3, 3), (0.1, 199.9)]:
        x = mp.mpf(x); y = mp.mpf(y)
        u = (x - y) / 2; v = (x + y) / 2
        rho = abs(x - y) / (2 + x + y)
        for m in (2, 5, 12):
            R = logGs(x, y) - sum(mp.mpf(u ** (2 * k)) / k * mp.zeta(2 * k, 1 + v)
                                  for k in range(1, m))
            lo = mp.mpf(u ** (2 * m)) / m * mp.zeta(2 * m, 1 + v)
            hi = rho ** (2 * m) / (m * (1 - rho ** 2)) * (1 + (1 + v) / (2 * m - 1))
            check(f"lo<=R<=hi at ({mp.nstr(x,6)},{mp.nstr(y,8)}), m={m}", lo <= R <= hi)
    x, y = mp.mpf(2), mp.mpf(7); u = (x - y) / 2; v = (x + y) / 2
    rho = abs(x - y) / (2 + x + y)
    vals = []
    for m in (10, 40, 160):
        R = sum(mp.mpf(u ** (2 * k)) / k * mp.zeta(2 * k, 1 + v) for k in range(m, m + 2500))
        vals.append(mp.log(R) / m)
    check("(1/m)log R_m increases towards 2 log rho",
          vals[0] < vals[1] < vals[2] < 2 * mp.log(rho),
          f"2log rho={mp.nstr(2*mp.log(rho),8)}, values={[mp.nstr(z,8) for z in vals]}")
    # Corollary 1.5 -- at m = 2 the upper bound of (8) = eq:Wbil improves on the upper
    # bound of Wisniewska's bilateral inequality.  The two ingredients of the
    # short proof are checked first: zeta(3,z) >= 1/(2z^2) and
    # v(A+3)/(3A^2) = v(v+4)/(3(v+1)^2) <= 4/9, which yield rho < sqrt(5)/3.
    okz = all(mp.zeta(3, z) * 2 * z ** 2 >= 1 for z in
              [mp.mpf(10) ** k for k in range(-3, 5)] + [mp.mpf('1.5'), 3, 7])
    okc = all(v * (v + 4) / (3 * (v + 1) ** 2) <= mp.mpf(4) / 9
              for v in [mp.mpf(k) / 20 for k in range(1, 4000)] + [mp.mpf(2)])
    check("Cor 1.5: zeta(3,z) >= 1/(2z^2) and v(A+3)/(3A^2) <= 4/9", okz and okc)
    okimp = True; nlt = 0; worst = mp.mpf(0)
    for v0 in [mp.mpf(s) for s in ('0.25', '0.5', '1', '2', '5', '20',
                                   '100', '1000', '100000')]:
        for j in range(1, 100):
            r0 = mp.mpf(j) / 100
            c = v0 * r0; A = 1 + v0; s = v0 * mp.sqrt(1 - r0 ** 2); rho = c / A
            if rho >= mp.sqrt(5) / 3:
                continue
            gain = c ** 2 * (mp.zeta(2, 1 + s) - mp.zeta(2, 1 + v0))
            corr = rho ** 4 / (2 * (1 - rho ** 2)) * (1 + A / 3)
            okimp &= (corr < gain) and (c ** 4 / (2 * v0 * A ** 2) <= gain)
            nlt += 1
            worst = max(worst, corr / gain)
    check("Cor 1.5: m=2 upper bound beats Wisniewska's for rho < sqrt(5)/3",
          okimp and nlt > 100, f"{nlt} pairs, worst corr/gain={mp.nstr(worst,4)}")

# ---------------------------------------------------------------- section 4
def sec4():
    print("\n[4] Proposition 1.4 -- comparison of the two rates")
    ok = True
    for x in GRID:
        for y in GRID:
            if x == y: continue
            xm = mp.mpf(x); ym = mp.mpf(y)
            rho = abs(xm - ym) / (2 + xm + ym)
            Q = abs(xm - ym) / (2 * (1 + mp.sqrt(xm * ym)))
            ok &= (rho <= Q)
            ok &= (rho < 1)
    check("rho <= Q and rho < 1 on the whole grid", ok)
    check("Q is unbounded: Q(1e8,1) > 1e3",
          abs(mp.mpf(10) ** 8 - 1) / (2 * (1 + mp.sqrt(mp.mpf(10) ** 8))) > 1000)

# ---------------------------------------------------------------- section 5
def sec5():
    print("\n[5] Theorem 1.6 -- large-argument asymptotics")
    prev = None
    for (x, y) in [(1e3, 1100), (1e4, 3e4), (1e6, 2e6)]:
        x = mp.mpf(x); y = mp.mpf(y)
        u = (x - y) / 2; v = (x + y) / 2
        rho = abs(u) / (1 + v)                 # Theorem 1.6 states rho = |u|/(1+v)
        pred = (1 + v) * Lam(rho)              # = (u^2/(1+v)) H(rho),  H = Lambda/rho^2
        # the two displayed forms of Theorem 1.6 agree exactly
        assert abs(pred - u * u / (1 + v) * (Lam(rho) / rho ** 2)) < mp.mpf('1e-18') * pred
        act = logGs(x, y)
        rel = abs(act - pred) / abs(act)
        check(f"relative error decreases at ({mp.nstr(x,4)},{mp.nstr(y,4)})",
              prev is None or rel < prev, f"rel={mp.nstr(rel,4)}")
        prev = rel

# ---------------------------------------------------------------- section 6
def sec6():
    print("\n[6] Theorem 1.7 -- multivariate expansion")
    ok = True
    for vv in ('1', '3'):
        v = mp.mpf(vv); e = mp.mpf('0.1')
        x1 = v * (1 + e); x2 = v * (1 + e); x3 = v * (1 - 2 * e)
        lhs = (mp.loggamma(x1) + mp.loggamma(x2) + mp.loggamma(x3)
               - 3 * mp.loggamma(v))
        ser = sum(((-1) ** j) / j * mp.zeta(j, v) * (v ** j * e ** j * (2 + (-2) ** j))
                  for j in range(2, 200))
        ok &= abs(lhs - ser) < mp.mpf(10) ** -18
    check("3-point expansion agrees with log Gamma to 18 digits", ok)

# ---------------------------------------------------------------- section 7
def sec7():
    print("\n[7] Theorem 1.8 -- localization of t")
    # (ii) bracket and the new lower-bound argument
    okbr = True; okmarg = True
    for (x, y) in [(1, 2), (0.5, 3), (2, 7), (1, 10), (20, 100), (0.1, 199.9), (100, 101)]:
        x = mp.mpf(x); y = mp.mpf(y)
        v = (x + y) / 2; c = abs(x - y) / 2; A = 1 + v
        rho = c / A
        t = tval(x, y)
        d = v - t
        z3 = mp.zeta(3, A); z4 = mp.zeta(4, A)
        d0 = c ** 2 * z4 / (4 * z3 + 6 * v * z4)
        hi = c ** 2 * z4 / (4 * z3) * (1 + 2 * rho ** 2 / (3 * (1 - rho ** 2)))
        okbr &= (d0 <= d <= hi)
        S0 = sum((j + 1) * d0 ** j * mp.zeta(j + 2, A) for j in range(3, 400))
        okmarg &= (S0 < 3 * d0 * z4 * (v - d0))
    check("lower/upper bracket eq:locbnd holds", okbr)
    check("key margin U_0 < 3 delta_0 zeta_4 (v-delta_0) of the lower-bound proof", okmarg)
    # (iii) c^4 coefficient
    x, y = mp.mpf(10), mp.mpf(11)
    v = (x + y) / 2; c = abs(x - y) / 2; A = 1 + v
    d = v - tval(x, y)
    z3 = mp.zeta(3, A); z4 = mp.zeta(4, A); z6 = mp.zeta(6, A)
    app = c ** 2 * z4 / (4 * z3) + c ** 4 / (2 * z3) * (z6 / 3 - 3 * z4 ** 3 / (16 * z3 ** 2))
    check("three-term formula eq:locasym matches v-t to O(c^6)", abs(app - d) < mp.mpf(10) ** -8,
          f"residual={mp.nstr(app-d,3)}")
    # Upper bound of the bracket at (0.1, 199.9), quoted as 514.83 in the text.
    x, y = mp.mpf('0.1'), mp.mpf('199.9')
    v = (x + y) / 2
    c = abs(x - y) / 2
    A = 1 + v
    rho = c / A
    z3, z4 = mp.zeta(3, A), mp.zeta(4, A)
    hi = c ** 2 * z4 / (4 * z3) * (1 + 2 * rho ** 2 / (3 * (1 - rho ** 2)))
    check("bracket upper bound at (0.1, 199.9) equals 514.83",
          abs(hi - mp.mpf('514.83')) < mp.mpf('0.005'), f"hi={mp.nstr(hi, 8)}")

# ---------------------------------------------------------------- section 8
def sec8():
    print("\n[8] Theorem 1.9 -- the sharp constant 2/3")
    mn = mp.inf; arg = None
    for x in GRID:
        for y in GRID:
            if x == y: continue
            xm = mp.mpf(x); ym = mp.mpf(y)
            t = tval(xm, ym)
            if t is None: continue
            v = (xm + ym) / 2; s = mp.sqrt(xm * ym)
            th = (t - s) / (v - s)
            if th < mn: mn, arg = th, (x, y)
    check("theta > 2/3 on the grid", mn > mp.mpf(2) / 3, f"min={mp.nstr(mn,10)} at {arg}")
    best = mp.inf
    for V in ('1e2', '1e3', '1e4', '1e5', '1e6', '1e7'):
        V = mp.mpf(V)
        for ei in ('1e-1', '1e-3', '1e-5'):
            e = mp.mpf(ei)
            x = V * (1 + e); y = V * (1 - e)
            t = tval(x, y)
            if t is None: continue
            s = mp.sqrt(x * y)
            best = min(best, (t - s) / (V - s))
    check("infimum path along x=V(1+e), y=V(1-e) stays above 2/3 and approaches it",
          best > mp.mpf(2) / 3 and best - mp.mpf(2) / 3 < mp.mpf('1e-6'),
          f"minimum on the path = {mp.nstr(best,12)}")
    # The value 2/3 is reached only in the joint limit v -> inf, eps -> 0,
    # v*eps -> 0.  For fixed eps the ratio tends to 2/3 + g(eps) with
    # g(eps) ~ eps^2/180; along eps = v^{-2} one has theta = 2/3 + 1/(6v) + O(v^{-2}).
    def theta_path(V, e):
        V = mp.mpf(V); e = mp.mpf(e)
        x = V * (1 + e); y = V * (1 - e)
        t = tval(x, y)
        if t is None:
            return None
        s = mp.sqrt(x * y)
        return (t - s) / (V - s)
    with mp.workdps(50):
        okf = True; vals = []
        for ei in ('1e-1', '1e-3'):
            seq = [theta_path(mp.mpf(10) ** k, ei) for k in (5, 6, 7, 8)]
            okf &= all(seq[i] > seq[i + 1] > mp.mpf(2) / 3 for i in range(len(seq) - 1))
            vals.append(f"eps={ei}: {mp.nstr(seq[-1],10)}")
        check("fixed eps: theta decreases to 2/3 + g(eps) > 2/3", okf, "; ".join(vals))
        # eps = v^{-3/2} is used rather than v^{-2}: both satisfy v*eps -> 0, but the
        # former keeps c = v*eps = v^{-1/2} large enough for the bisection in
        # tval() to bracket reliably at 50 digits.
        # k <= 6 only: beyond that delta = v-t is so small relative to v that the
        # bisection in tval() reaches its precision floor and returns noise.
        okj = True; n = 0; parts = []
        for k in (4, 5, 6):
            V = mp.mpf(10) ** k
            th = theta_path(V, V ** mp.mpf('-1.5'))
            if th is None:
                continue
            r = V * (th - mp.mpf(2) / 3)
            okj &= abs(r - mp.mpf(1) / 6) < mp.mpf('0.02')
            n += 1
            parts.append(f"V=1e{k}: {mp.nstr(r,8)}")
        check("joint limit v*eps->0: theta = 2/3 + 1/(6v) + O(v^-2)",
              okj and n >= 3, f"v(theta-2/3) -> 1/6: " + ", ".join(parts))
        thg = theta_path(mp.mpf('1e7'), mp.mpf('1e-1'))
        g = thg - mp.mpf(2) / 3
        check("g(1e-1) at v=1e7 is 5.60e-5, and g/eps^2 is near 1/180",
              abs(g - mp.mpf('5.60e-5')) < mp.mpf('5e-8')
              and abs(g / mp.mpf('1e-2') - mp.mpf(1) / 180) < mp.mpf('5e-5'),
              f"g={mp.nstr(g,6)}, g/eps^2={mp.nstr(g / mp.mpf('1e-2'),6)}")
    # Product formula: G_* = prod_{m>=0} (A+m)^2/((A+m)^2-c^2).
    okpr = True; parts = []
    for (xi, yi) in [(1, 2), (1000, 999), (3, '3.0001'), (1, 10), (20, 100)]:
        x0 = mp.mpf(xi); y0 = mp.mpf(yi); v0 = (x0 + y0) / 2; c0 = abs(x0 - y0) / 2
        A0 = 1 + v0
        direct = (mp.gamma(x0) * mp.gamma(y0) / mp.gamma(v0) ** 2
                  * 4 * x0 * y0 / (x0 + y0) ** 2)
        N = 40000
        Pr = mp.mpf(1)
        for m in range(0, N):
            u = A0 + m
            Pr *= u * u / (u * u - c0 * c0)
        Pr *= mp.e ** (c0 * c0 / (A0 + N - mp.mpf(1) / 2))   # integral tail of the log-sum
        okpr &= abs(Pr / direct - 1) < mp.mpf('1e-4')
        parts.append(f"({xi},{yi}): rel={mp.nstr(abs(Pr/direct-1),2)}")
    check("product formula for G_*", okpr, "; ".join(parts))
    # Scaling limit: D(rho) = kappa/(1-kappa) - P(rho^2) + 1 > 0 on (0,1),
    # with kappa = (1-sqrt(1-rho^2))/3 and P(rho^2) = 2 artanh(rho)/rho + log(1-rho^2)/rho^2.
    okd = True; wd = mp.inf; ar = None
    for j in range(1, 2000):
        r0 = mp.mpf(j) / 2000
        t0 = (1 - mp.sqrt(1 - r0 ** 2)) / 3
        S0 = 2 * mp.atanh(r0) / r0 + mp.log(1 - r0 ** 2) / r0 ** 2
        dd = t0 / (1 - t0) - S0 + 1
        okd &= (dd > 0)
        if dd < wd:
            wd, ar = dd, r0
    check("scaling-limit inequality D(rho) > 0 on (0,1)",
          okd, f"min D = {mp.nstr(wd,6)} at rho = {mp.nstr(ar,4)} "
               f"(positivity is now proved in Step 2: a_m > b_m for all m >= 2)")
    # NOTE: a "two-term expansion  A*R = D(rho) + E1(rho)/A + O(A^-2)" was checked
    # here in an earlier draft.  It is NOT retained: two independent scripts of
    # mine disagreed on the sign of A*R - D at rho = 0.02, A = 1e4, so the
    # expansion as written is not trustworthy.  What IS verified is the direct
    # positivity of R and of theta - 2/3 (see the checks above).
    # The tempting bound  E(c) <= (c^2/2) sum_n 1/(u_n^2(u_n^2-c^2))  is FALSE:
    # it would require m/(m+1) <= 1/2 termwise, which already fails at m = 2
    # (2/3 > 1/2).  It was used in an earlier draft and is checked to fail here,
    # so that the mistake cannot silently return.
    v0 = mp.mpf(6); c0 = mp.mpf('0.5'); A0 = 1 + v0
    lhs = sum(mp.mpf(m) / (m + 1) * c0 ** (2 * m) * mp.zeta(2 * m + 2, A0)
              for m in range(1, 400))
    rhs = c0 ** 2 / 2 * (sum(1 / ((A0 + n) ** 2 * ((A0 + n) ** 2 - c0 ** 2))
                             for n in range(0, 20000)) + 1 / (3 * (A0 + 20000) ** 3))
    check("the discarded bound E <= (c^2/2) sum 1/(u^2(u^2-c^2)) indeed fails",
          lhs > rhs, f"E={mp.nstr(lhs,6)} > bound={mp.nstr(rhs,6)}")

# ---------------------------------------------------------------- section 9
def sec9():
    print("\n[9] Section 3 -- applications")
    okg = True
    for (x, y) in [(1, 2), (0.5, 3), (2, 7), (1, 10), (0.1, 199.9), (100, 101)]:
        x = mp.mpf(x); y = mp.mpf(y); v = (x + y) / 2
        rho = abs(x - y) / (x + y)
        L = logG(x, y)
        lo = (x - y) ** 2 / 4 * mp.zeta(2, v)
        hi = v * Lam(rho) - mp.log(1 - rho ** 2)
        okg &= (lo <= L <= hi)
    check("Corollary 3.1 (Gurland two-sided bound)", okg)
    okb = True
    for (x, y) in [(1, 2), (1, 1), (0.5, 3), (2, 7), (10, 20)]:
        x = mp.mpf(x); y = mp.mpf(y); v = (x + y) / 2
        rho = abs(x - y) / (x + y)
        B = mp.gamma(x) * mp.gamma(y) / mp.gamma(x + y)
        lo = mp.sqrt(mp.pi) * v ** mp.mpf('-0.5') / 2 ** (2 * v - 1) * mp.e ** ((x - y) ** 2 / 4 * mp.zeta(2, v))
        hi = (mp.sqrt(mp.pi) * (v + 1) ** mp.mpf('0.5') / (2 ** (2 * v - 1) * v)
              * mp.e ** (v * Lam(rho) - mp.log(1 - rho ** 2)))
        okb &= (lo < B < hi)
    check("Corollary 3.2 (beta-function bounds)", okb)

SECTIONS = {1: sec1, 2: sec2, 3: sec3, 4: sec4, 5: sec5, 6: sec6, 7: sec7, 8: sec8, 9: sec9}

if __name__ == "__main__":
    want = [int(a) for a in sys.argv[1:]] or sorted(SECTIONS)
    t0 = time.time()
    print("=" * 74)
    print("Verification suite for the Gurland-ratio paper   (mpmath, %d digits)" % mp.mp.dps)
    print("=" * 74)
    for s in want:
        SECTIONS[s]()
    bad = [r for r in RESULTS if not r[1]]
    print("\n" + "=" * 74)
    print(f"SUMMARY: {len(RESULTS) - len(bad)}/{len(RESULTS)} checks passed "
          f"in {time.time()-t0:.1f} s")
    for tag, _, detail in bad:
        print("   FAILED:", tag, detail)
    print("=" * 74)
    sys.exit(1 if bad else 0)
