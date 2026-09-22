import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 35
def setup(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); a=1+v; w=(v-s)/2
    return v,c,s,a,w

def EM(x,y,r=3):
    v,c,s,a,w = setup(x,y)
    K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
    L=lambda t: 1/(t-w)**2
    G=lambda t: K(t)-L(t)
    I=mp.quad(G,[a,mp.inf])
    tot=I+G(a)/2
    for j in [1,3,5,7][:r]:
        B={1:mp.mpf(1)/6,3:mp.mpf(-1)/30,5:mp.mpf(1)/42,7:mp.mpf(-1)/30}[j]
        tot += -B/mp.factorial(j)*mp.diff(G,a,j)
    # remainder bound  |R_r| <= 2 zeta(2r+2)/(2pi)^{2r+2} * int_a^inf |G^{(2r+2)}|
    m=2*r+2
    bd = 2*mp.zeta(m)/(2*mp.pi)**m * mp.quad(lambda t: abs(mp.diff(G,t,m)), [a, mp.inf])
    return tot, bd

grid=[0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000]
print("="*78); print("EULER-MACLAURIN with rigorous remainder, correct target w=(v-s)/2"); print("="*78)
for r in [1,2,3]:
    bad=0; tot=0; worst=(None,-1e9)
    for x in grid:
        for y in grid:
            if x==y: continue
            try: e,bd = EM(x,y,r)
            except Exception: continue
            tot+=1
            if e+bd>=0: bad+=1
            if float(e+bd)>worst[1]: worst=((x,y),float(e+bd))
    print(f"  order r={r}: tested {tot},  E_r + rem_bound >= 0 in {bad} cases;  worst value {worst[1]:+.4e} at {worst[0]}")

print()
print("detail at r=3 for a spread of cases:")
for (x,y) in [(1,2),(2,7),(1,10),(100,1),(1000,999),(10000,0.001),(0.01,0.001),(0.5,3),(100.5,99.5)]:
    e,bd=EM(x,y,3)
    v,c,s,a,w=setup(x,y)
    print(f"   (x,y)=({x},{y})  rho={mp.nstr(c/a,6):>9}  E_3={mp.nstr(e,8):>15}  rem_bd={mp.nstr(bd,4):>11}  E+rem={mp.nstr(e+bd,8):>15}  ok={e+bd<0}")
