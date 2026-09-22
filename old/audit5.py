import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40

def dat(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); a=1+v; w=(v-s)/2
    return v,c,s,a,w

def intK(A,c):
    # int_A^inf K(x) dx ,  K(x)=c^-2 ln(1/(1-c^2/x^2)),  A>c
    return (-A*mp.log(A**2/(A**2-c**2))+c*mp.log((A+c)/(A-c)))/c**2

print("="*78)
print("CORRECT TARGET:  Phi < zeta(2, 1+(v+s)/2)   i.e. t > (v+s)/2")
print("  w=(v-s)/2 ;  sum_{m>=1}[K-L](m+v) should be < 0")
print("="*78)
grid=[0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000]
bad_true=0; bad_bnd=0; tot=0; worst=[]
for xv in grid:
    for yv in grid:
        if xv==yv: continue
        v,c,s,a,w = dat(xv,yv)
        K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
        L=lambda t: 1/(t-w)**2
        true=sum(K(m+v)-L(m+v) for m in range(1,3001))+mp.quad(lambda t: K(t)-L(t),[3000+v,mp.inf])
        # convexity bound: sum K <= int_{v+1/2}^inf K ;  sum L >= 1/(a-w) + 1/(2(a-w)^2)
        bnd = intK(v+mp.mpf(1)/2, c) - 1/(a-w) - 1/(2*(a-w)**2)
        tot+=1
        if true>=0: bad_true+=1
        if bnd>=0: bad_bnd+=1
        worst.append((float(bnd),float(true),xv,yv))
worst.sort(reverse=True)
print(f"  tested {tot}:  true sum >= 0 in {bad_true};  convexity bound >= 0 in {bad_bnd}")
print("  worst 10 (largest convexity bounds):")
for t in worst[:10]:
    print(f"     bound={t[0]:+.6e}   true={t[1]:+.6e}   (x,y)=({t[2]},{t[3]})")

print()
print("="*78)
print("Head-to-head on extreme cases :  true sum  vs  convexity bound")
print("="*78)
for (x,y) in [(10000,0.001),(0.001,10000),(1,10),(100,1),(2,7),(1,2),(1000,999),(100.5,99.5),(1e6,1e6-1)]:
    v,c,s,a,w = dat(x,y)
    K=lambda t: mp.log(1/(1-c**2/t**2))/c**2
    L=lambda t: 1/(t-w)**2
    true=sum(K(m+v)-L(m+v) for m in range(1,3001))+mp.quad(lambda t: K(t)-L(t),[3000+v,mp.inf])
    bnd = intK(v+mp.mpf(1)/2, c) - 1/(a-w) - 1/(2*(a-w)**2)
    print(f"  (x,y)=({x},{y})  rho={mp.nstr(c/a,6):>9}  true={mp.nstr(true,8):>15}  bound={mp.nstr(bnd,8):>15}  bound<0: {bnd<0}")
