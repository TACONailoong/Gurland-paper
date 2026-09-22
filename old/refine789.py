
import mpmath as mp
mp.mp.dps = 30
GRID = [0.001, 0.01, 0.1, 0.5, 1, 2, 5, 20, 100, 1000, 10000]

def condA(x, y):
    v = (x + y) / 2; s = mp.sqrt(x * y); A = 1 + v
    rho = abs(x - y) / (2 + x + y)
    return (v + s) / (2 * (1 + v)) * (1 + 2 * rho ** 2 / (3 * (1 - rho ** 2))) < 1

def condB(x, y):
    v = (x + y) / 2; c = abs(x - y) / 2; s = mp.sqrt(x * y); A = 1 + v
    w = (v - s) / 2
    B0 = v + mp.mpf(1) / 2
    closed = (-B0 * mp.log(B0 ** 2 / (B0 ** 2 - c ** 2))
              + c * mp.log((B0 + c) / (B0 - c))) / c ** 2
    return closed - 1 / (A - w) - 1 / (2 * (A - w) ** 2) < 0

nA = nB = nN = 0; ntot = 0
both = 0
for x in GRID:
    for y in GRID:
        if x == y: continue
        ntot += 1
        a, b = condA(mp.mpf(x), mp.mpf(y)), condB(mp.mpf(x), mp.mpf(y))
        nA += a; nB += b; nN += (not a and not b); both += (a and b)
print("off-diagonal pairs:", ntot, "A:", nA, "B:", nB, "both:", both, "neither:", nN)

# diagonal sanity: rho = 0, condition (A) reduces to v/(1+v) < 1
print("diagonal x=y: condA true?", all(condA(mp.mpf(t), mp.mpf(t)) for t in GRID))

# Remark 1.15 : D(rho) = h(rho^2) - S(rho^2) + 1
def h(y):
    q = mp.sqrt(1 - y); return (1 - q) / (2 + q)
def S(y):
    return mp.nsum(lambda m: y ** m / ((m + 1) * (2 * m + 1)), [0, mp.inf])
def D(rho):
    y = mp.mpf(rho) ** 2
    return h(y) - S(y) + 1
for rho in ['0.001', '0.01', '0.5', '0.9', '0.99', '0.999']:
    r = mp.mpf(rho)
    print("D(%s) = %s   rho^4/360 = %s" % (rho, mp.nstr(D(r), 12), mp.nstr(r**4/360, 12)))
print("t/(1-t) == h(rho^2)?", mp.nstr(mp.mpf('0.7')/(1-mp.mpf('0.7')) - h(mp.mpf('0.7')**2), 5))
# h' - S' > 0
for y in ['0.01', '0.3', '0.7', '0.95', '0.999']:
    yy = mp.mpf(y)
    hp = 3/(2*mp.sqrt(1-yy)*(2+mp.sqrt(1-yy))**2)
    Sp = mp.nsum(lambda m: m*yy**(m-1)/((m+1)*(2*m+1)), [1, mp.inf])
    print("y=%s  h'-S' = %s" % (y, mp.nstr(hp-Sp, 8)))
# series expansion of D
ser = mp.taylor(lambda z: h(z) - S(z) + 1, 0, 6)
print("taylor coeffs:", [mp.nstr(c, 8) for c in ser])
