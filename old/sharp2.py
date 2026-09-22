import mpmath as mp
mp.mp.dps = 25
def ratio(v, c, N=30000):
    v=mp.mpf(v); c=mp.mpf(c); A=1+v; s=mp.sqrt(v*v-c*c); D=(v-s)/3
    S=mp.mpf(0)
    for n in range(0,N):
        u=A+n; S+=1/(u*u*(u*u-c*c))
    S+=1/(3*(A+N)**3)
    return mp.zeta(3,A-D)/(mp.mpf(3)/2*s*S)
worst=(mp.inf,None,None)
for vi in ['1.01','1.5','2','3','10','100','10000']:
    v=mp.mpf(vi); mn=(mp.inf,None)
    for j in range(1,41):
        c=v*mp.mpf(j)/40
        r=ratio(v,c)
        if r<mn[0]: mn=(r,c)
    print(f"  v={mp.nstr(v,8):>9}: min LHS/RHS = {mp.nstr(mn[0],8)} at c/v={mp.nstr(mn[1]/v,4)}", flush=True)
    if mn[0]<worst[0]: worst=(mn[0],v,mn[1])
print("worst:", mp.nstr(worst[0],8), "at v=", mp.nstr(worst[1],6), "c/v=", mp.nstr(worst[2]/worst[1],4))
