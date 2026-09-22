import sys, mpmath as mp
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
mp.mp.dps = 40

def intK(A,c):
    return (-A*mp.log(A**2/(A**2-c**2))+c*mp.log((A+c)/(A-c)))/c**2

def analyse(x,y):
    x=mp.mpf(x); y=mp.mpf(y)
    v=(x+y)/2; c=abs(x-y)/2; s=mp.sqrt(x*y); A=1+v; w=(v-s)/2; rho=c/A
    # Method A  (series / two-term comparison)
    condA = 2*rho**2/(3*(1-rho**2)) < (2+v-s)/(v+s)
    # Method B  (convexity: midpoint rule for K, trapezoid for L)
    bnd = intK(v+mp.mpf(1)/2, c) - 1/(A-w) - 1/(2*(A-w)**2)
    condB = bnd < 0
    return float(rho), bool(condA), bool(condB), float(bnd)

print("="*76)
print("CASE-SPLIT COVERAGE :  condition A (series) OR condition B (convexity)")
print("="*76)
g1=[0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000]
g2=[0.001,0.01,0.1,0.5,1,2,5,20,100,1000,10000,100000]
rows=[]
for x in g1:
    for y in g2:
        if x==y: continue
        rows.append(analyse(x,y))
n=len(rows); onlyA=sum(1 for r in rows if r[1] and not r[2]); onlyB=sum(1 for r in rows if r[2] and not r[1])
both=sum(1 for r in rows if r[1] and r[2]); neither=[r for r in rows if not r[1] and not r[2]]
print(f"  tested {n}:  A only {onlyA},  B only {onlyB},  both {both},  NEITHER {len(neither)}")
for r in neither[:12]: print("    NEITHER rho=%.6f  A=%s B=%s  bound=%+.3e"%(r[0],r[1],r[2],r[3]))

print()
print("  behaviour vs rho  (sorted):")
rows.sort()
import itertools
buckets={}
for rho,A,B,_ in rows:
    k=round(min(max(rho,0.0),1.0),1)
    d=buckets.setdefault(k,[0,0,0])
    d[0]+=1
    if A: d[1]+=1
    if B: d[2]+=1
print("   rho-bucket   n    A-holds  B-holds")
for k in sorted(buckets):
    d=buckets[k]; print(f"   {k:>8.1f}  {d[0]:>4}   {d[1]:>5}   {d[2]:>5}")

print()
print("="*76)
print("fine scan of the boundary:  for each rho, is A or B satisfied?")
print("="*76)
bad=0; tot=0
for i in range(1,400):
    rho=mp.mpf(i)/400
    # choose a configuration with this rho: v large so that w,v free; take c=rho*(1+v)
    V=mp.mpf(200)
    c=rho*(1+V); 
    if c>=V: continue
    s=mp.sqrt(V*V-c*c)
    x=V+c; y=V-c
    rr,A,B,_=analyse(x,y); tot+=1
    if not (A or B): bad+=1
print(f"  V=200, rho in (0,1): tested {tot}, NEITHER in {bad}")
bad=0; tot=0
for i in range(1,400):
    rho=mp.mpf(i)/400
    V=mp.mpf('5000')
    c=rho*(1+V)
    if c>=V: continue
    x=V+c; y=V-c
    rr,A,B,_=analyse(x,y); tot+=1
    if not (A or B): bad+=1
print(f"  V=5000, rho in (0,1): tested {tot}, NEITHER in {bad}")
