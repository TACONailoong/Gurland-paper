import mpmath as mp
mp.mp.dps = 40
print("--- product formula:  G_* = prod_{m>=0} (A+m)^2/((A+m)^2-c^2) ? ---")
def Gstar(x,y):
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; c=abs(x-y)/2
    return mp.gamma(x)*mp.gamma(y)/mp.gamma(v)**2 * (4*x*y)/(x+y)**2
for (x,y) in [(1,2),(20,100),(1000,999),(0.1,199.9),(3,3.0001),(1,10)]:
    x=mp.mpf(x); y=mp.mpf(y); v=(x+y)/2; c=abs(x-y)/2; A=1+v
    direct = Gstar(x,y)
    N = 20000
    P = mp.mpf(1)
    for m in range(0,N):
        u=A+m; P *= u*u/(u*u-c*c)
    # tail: prod (1 + c^2/((u^2-c^2))) ~ exp(sum c^2/u^2) over tail
    P *= mp.e**(c*c/(A+N))
    print(f"  ({mp.nstr(x,6)},{mp.nstr(y,6)}): direct={mp.nstr(direct,18)}  product={mp.nstr(P,18)}  rel.diff={mp.nstr(abs(P/direct-1),3)}")
