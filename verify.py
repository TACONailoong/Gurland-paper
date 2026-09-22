#!/usr/bin/env python3
"""
Numerical verification suite for
  "The exact expansion of Gurland's ratio of the gamma function,
   with sharp remainder estimates and a localization theorem".

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
    print("\n[3] Theorem 1.3 -- sharp two-sided remainder estimate and optimality")
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
    check("lower/upper bracket of (19) = eq:locbnd holds", okbr)
    check("key margin S'_0 < 3 delta_0 zeta_4 (v-delta_0) of the lower-bound proof", okmarg)
    # (iii) c^4 coefficient
    x, y = mp.mpf(10), mp.mpf(11)
    v = (x + y) / 2; c = abs(x - y) / 2; A = 1 + v
    d = v - tval(x, y)
    z3 = mp.zeta(3, A); z4 = mp.zeta(4, A); z6 = mp.zeta(6, A)
    app = c ** 2 * z4 / (4 * z3) + c ** 4 / (2 * z3) * (z6 / 3 - 3 * z4 ** 3 / (16 * z3 ** 2))
    check("three-term formula (20) = eq:locasym matches v-t to O(c^6)", abs(app - d) < mp.mpf(10) ** -8,
          f"residual={mp.nstr(app-d,3)}")
    # (iv) coverage of (A) and (B)
    def condA(x, y):
        v = (x + y) / 2; s = mp.sqrt(x * y); A = 1 + v
        rho = abs(x - y) / (2 + x + y)
        return (v + s) / (2 * (1 + v)) * (1 + 2 * rho ** 2 / (3 * (1 - rho ** 2))) < 1
    def condB(x, y):
        # E(x,y) = int_{v+1/2}^inf K - 1/(A-w) - 1/(2(A-w)^2),  w=(v-s)/2, A=1+v
        v = (x + y) / 2; c = abs(x - y) / 2; s = mp.sqrt(x * y); A = 1 + v
        w = (v - s) / 2
        B0 = v + mp.mpf(1) / 2
        K = lambda t: mp.log(1 / (1 - c ** 2 / t ** 2)) / c ** 2
        closed = (-B0 * mp.log(B0 ** 2 / (B0 ** 2 - c ** 2))
                  + c * mp.log((B0 + c) / (B0 - c))) / c ** 2
        assert abs(mp.quad(K, [B0, mp.inf]) - closed) < mp.mpf(10) ** -12, "closed form for int K"
        return closed - 1 / (A - w) - 1 / (2 * (A - w) ** 2) < 0
    nA = nB = nN = ntot = 0
    for x in GRID:
        for y in GRID:
            if x == y: continue   # rho = 0 on the diagonal; not counted (Table 4)
            ntot += 1
            a, b = condA(mp.mpf(x), mp.mpf(y)), condB(mp.mpf(x), mp.mpf(y))
            nA += a; nB += b; nN += (not a and not b)
    check("neither (A) nor (B) occurs: never", nN == 0 and ntot == 110,
          f"110 off-diagonal pairs of the 121-point grid: A={nA}, B={nB}, neither={nN}")
    # The two additional families of Table 4, so that the whole table is
    # reproducible from this script: fixed v and shape parameter r,
    # x = v(1+r), y = v(1-r), with r = j/400 for j = 1, ..., 399.
    fam = []
    for v0 in (200, 5000):
        v0 = mp.mpf(v0); a = b = nn = tot = 0
        for j in range(1, 400):
            r0 = mp.mpf(j) / 400
            x = v0 * (1 + r0); y = v0 * (1 - r0)
            ca = condA(x, y); cb = condB(x, y)
            a += ca; b += cb; nn += (not ca and not cb); tot += 1
        fam.append((v0, tot, a, b, nn))
    check("Table 4 families: (A) or (B) holds in every case",
          all(f[4] == 0 for f in fam),
          "; ".join(f"v={mp.nstr(f[0], 6)}: {f[1]} cases, A={f[2]}, B={f[3]}, "
                    f"neither={f[4]}" for f in fam))
    # ---- Euler-Maclaurin certificate -------------------------------------
    # For f in C^{2r}[a,inf) with f^{(j)}(inf) = 0 for j <= 2r-1,
    #   sum_{n>=0} f(a+n) = int_a^inf f + f(a)/2
    #                       - sum_{k=1..r} B_{2k}/(2k)! f^{(2k-1)}(a) + R_r ,
    #   R_r = -1/(2r)! int_a^inf B~_{2r}(x) f^{(2r)}(x) dx ,
    #   |R_r| <= |B_{2r}|/(2r)! int_a^inf |f^{(2r)}|
    #          = 2 zeta(2r)/(2 pi)^{2r} int_a^inf |f^{(2r)}| ,
    # because the periodic Bernoulli polynomial B~_{2r} attains its extreme
    # value B_{2r} at the endpoints.  The remainder therefore carries the SAME
    # index 2r as the last coefficient included; an earlier draft paired the
    # k <= r sum with int |f^{(2r+2)}|, which is not the classical statement.
    # NOTE: an earlier draft also divided B_{2k} by (2k-1)! instead of (2k)!;
    # the coefficients used here are -1/12, +1/720, -1/30240 for k = 1,2,3.
    B2K = {1: mp.mpf(1) / 6, 2: mp.mpf(-1) / 30, 3: mp.mpf(1) / 42,
           4: mp.mpf(-1) / 30, 5: mp.mpf(5) / 66, 6: mp.mpf(-691) / 2730}
    # The certificate uses the orders r = 1..5 that it reports as certifying;
    # the reference value below goes up to k = 6 so that its own truncation
    # (~T0^{-15}) is far below the smallest remainder bound (~4e-35) that it has
    # to resolve.

    # Exact derivatives, obtained symbolically with sympy and then evaluated in
    # mpmath.  Numerical differentiation (mp.diff) is not reliable at these
    # orders and silently produced a remainder bound that was far too small.
    import sympy as _sp
    _t, _cc = _sp.symbols('t c', positive=True)
    _Kexpr = -_sp.log(1 - _cc ** 2 / _t ** 2) / _cc ** 2
    _ORD = list(range(14))
    _Kfun = {j: _sp.lambdify((_t, _cc), _sp.diff(_Kexpr, _t, j), modules='mpmath')
             for j in _ORD}
    _Kint = _sp.lambdify((_t, _cc), _sp.integrate(_Kexpr, (_t, _t, _sp.oo)),
                         modules='mpmath')          # int_t^inf K, as a function of (t,c)

    def Kder(j, t, c):
        return _Kfun[j](mp.mpf(t), mp.mpf(c))

    def Lder(j, t, w):                              # L(t) = (t-w)^{-2}
        return (-1) ** j * mp.factorial(j + 1) * mp.mpf(t - w) ** (-(j + 2))

    def Gder(j, t, c, w):
        return Kder(j, t, c) - Lder(j, t, w)

    def em_data(x, y, orders=(1, 2, 3), with_true=False):
        # x and y must be mpf *before* any arithmetic: _Kint evaluates
        # log((t^2-c^2)/t^2), and with binary-float arguments the subtraction
        # 1 - c^2/t^2 is performed in double precision, which destroys ~10 digits
        # when c << t (relative error ~2e-10 at (1000,999)).  This silently
        # corrupted the reference value in earlier drafts.
        x = mp.mpf(x); y = mp.mpf(y)
        v = (x + y) / 2; c = abs(x - y) / 2; s = mp.sqrt(x * y); A = 1 + v
        w = (v - s) / 2
        G = lambda t: Gder(0, t, c, w)
        # int_A^inf G is evaluated in closed form; quadrature of G=A(K-L) only
        # reaches ~1e-10 relative accuracy here, which is far too coarse to
        # validate a 1e-35 remainder bound.  Quadrature is kept as a cross-check.
        I = _Kint(A, c) - 1 / (A - w)
        assert abs(mp.quad(G, [A, mp.inf]) - I) < mp.mpf(10) ** -8 * abs(I), \
            "closed form for int G disagrees with quadrature"
        g0 = G(A)
        out = {}
        for r in orders:
            E = I + g0 / 2
            for k in range(1, r + 1):
                E += -B2K[k] / mp.factorial(2 * k) * Gder(2 * k - 1, A, c, w)
            m = 2 * r
            Rem = (2 * mp.zeta(m) / (2 * mp.pi) ** m
                   * mp.quad(lambda t: abs(Gder(m, t, c, w)), [A, mp.inf]))
            out[r] = (E, Rem)
        true = None
        if with_true:
            # Reference value of sum_{n>=0} G(a+n).  An earlier draft used
            # "sum_{m<=N} G + int_{N+v}^inf G", which is only accurate to
            # |G(N+v)|/2 and therefore could not validate the bound; here the
            # same Euler-Maclaurin formula is applied at T0 = N+v, chosen with
            # N ~ 3c so that c/T0 <= 1/4 and the expansion there converges
            # geometrically.  The result is exact to far more digits than the
            # certificate needs.
            # The reference is built from the two halves separately, in the only
            # way that avoids cancellation:  sum_{m>=1} L(m+v) = zeta(2,1+(v+s)/2)
            # is *exact* in mpmath, and sum_{m>=1} K(m+v) = Phi is a sum of
            # positive decreasing terms, so Euler-Maclaurin (applied at T0 = N+v,
            # where c/T0 is small) converges geometrically and is applied to a
            # quantity with no cancellation at all.  Summing (K-L) directly would
            # cancel ~11 digits and could not validate the certificate.
            N = max(200, 5 * int(c) + 1)
            T0 = N + v
            Phi = sum(Kder(0, m + v, c) for m in range(1, N + 1))
            Phi += _Kint(T0, c) - Kder(0, T0, c) / 2
            for k in range(1, 7):
                Phi += -B2K[k] / mp.factorial(2 * k) * Kder(2 * k - 1, T0, c)
            true = Phi - mp.zeta(2, 1 + (v + s) / 2)
        return out, true

    # The same 11-point grid as the coverage check above, so that "0 of 110" in
    # Table 6 of the paper refers to the identical set of configurations.
    pairs = [(x, y) for x in GRID for y in GRID if x != y]
    RMAX = 6
    ORDERS = tuple(range(1, RMAX + 1))
    # r must not exceed the optimal truncation index: for r = 6 the bound R_r
    # has already started to grow again for some configurations (it is 5.1e-3 at
    # (0.01,0.1)), so it no longer certifies negativity there.  This is the
    # expected behaviour of an asymptotic expansion, not a failure of the bound;
    # the certificate is therefore claimed for the orders that do certify.
    worst = {r: (-mp.inf, None) for r in ORDERS}
    npos = {r: 0 for r in ORDERS}
    for (x, y) in pairs:
        data, _ = em_data(x, y, orders=ORDERS)
        for r in ORDERS:
            val = data[r][0] + data[r][1]
            if val > worst[r][0]:
                worst[r] = (val, (x, y))
            if val >= 0:
                npos[r] += 1
    CERTR = [r for r in ORDERS if npos[r] == 0]
    check("Euler-Maclaurin certificate is negative for all sampled configurations",
          len(CERTR) > 0,
          f"certifying orders {CERTR} of {len(pairs)} configurations; worst "
          + mp.nstr(worst[CERTR[-1]][0], 6) + f" at {worst[CERTR[-1]][1]}"
          + (f"; first failing order: r={min(r for r in ORDERS if npos[r])}"
             if len(CERTR) < len(ORDERS) else ""))
    # The certification itself runs at 35 digits, which is ample there: the
    # smallest remainder bound is ~4e-35.  The validation below must resolve
    # that, so it is carried out at 60 digits -- at 35 digits the reference
    # value is itself uncertain by ~4e-35 (the term K(T0) is computed as
    # log(1/(1-c^2/T0^2)), a cancellation of 1-6e-10).
    with mp.workdps(60):
        okv = True; parts = []
        for (x, y) in [(0.01, 10000), (1, 2), (100, 1), (1000, 999)]:
            data, true = em_data(x, y, orders=ORDERS, with_true=True)
            for r in ORDERS:
                E, Rem = data[r]
                okv &= (abs(E - true) <= Rem)
            parts.append(f"({x},{y}) |E_r-true|<={mp.nstr(max(abs(data[r][0]-true) / data[r][1] for r in ORDERS),4)}*R_r")
        check("EM remainder bound validated against the true sum", okv, "; ".join(parts))

# ---------------------------------------------------------------- section 8
def sec8():
    print("\n[8] Theorem 1.9 -- the sharp constant 2/3 (ex Conjecture 1.9)")
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
    # The sharp value 2/3 is reached only in the JOINT limit v -> inf, eps -> 0,
    # v*eps -> 0.  For a fixed eps the ratio tends to 2/3 + g(eps) with
    # g(eps) > 0 (numerically g(eps) ~ 5.6 eps^2); along eps = v^-2 the approach
    # is theta = 2/3 + 1/(6v) + O(v^-2).  These two statements are what the
    # corrected Table 5 records, and they are checked here.
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
    # Theorem 1.11 (partial resolution): the two estimates it rests on, and the
    # conclusion theta > 2/3 on the range c^2 <= A/6.
    okp = True
    for k in range(0, 4000):
        A = mp.mpf(2) + mp.mpf(k) / 50
        Psi = (A - 1) * A ** 3 / ((A - mp.mpf(1) / 2) ** 3 * (A + 1))
        okp &= (Psi <= 1 - 1 / (3 * A))
    okr = True
    for k in range(0, 2000):
        A = mp.mpf(3) / 2 + mp.mpf(k) / 20
        okr &= (mp.zeta(4, A) / mp.zeta(3, A)
                <= 2 * A ** 3 / (3 * (A - mp.mpf(1) / 2) ** 3 * (A + 1)))
    check("Thm 1.11: Psi(A) <= 1-1/(3A) and the bound (36) = eq:zetabnd for zeta4/zeta3", okp and okr)
    okc = True; nc = 0
    for v in ('2', '3', '5', '10', '50', '200', '1000', '10000'):
        v = mp.mpf(v); A = 1 + v
        for j in range(1, 13):
            c = mp.sqrt(A / 6) * mp.mpf(j) / 12
            if c >= v:
                continue
            x = v + c; y = v - c
            t = tval(x, y)
            if t is None:
                continue
            th = (t - mp.sqrt(x * y)) / (v - mp.sqrt(x * y))
            okc &= (th > mp.mpf(2) / 3)
            nc += 1
    check("Thm 1.11: theta > 2/3 throughout the range (x-y)^2 <= (2+x+y)/3",
          okc and nc > 80, f"{nc} pairs tested, all above 2/3")
    # Remark 1.12 (product formula): G_* = prod_{m>=0} (A+m)^2/((A+m)^2-c^2).
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
    check("Remark 1.12: product formula for G_*", okpr, "; ".join(parts))
    # Remark 1.15 (scaling limit): D(rho) = t/(1-t) - Sigma_0(rho) + 1 > 0 on (0,1),
    # with t = (1-sqrt(1-rho^2))/3 and Sigma_0 = 2 artanh(rho)/rho + log(1-rho^2)/rho^2.
    okd = True; wd = mp.inf; ar = None
    for j in range(1, 2000):
        r0 = mp.mpf(j) / 2000
        t0 = (1 - mp.sqrt(1 - r0 ** 2)) / 3
        S0 = 2 * mp.atanh(r0) / r0 + mp.log(1 - r0 ** 2) / r0 ** 2
        dd = t0 / (1 - t0) - S0 + 1
        okd &= (dd > 0)
        if dd < wd:
            wd, ar = dd, r0
    check("Remark 1.15: the scaling-limit inequality D(rho) > 0 on (0,1)",
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
