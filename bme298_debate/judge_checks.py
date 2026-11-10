import numpy as np
from scipy import stats
rng=np.random.default_rng(1)
# A/B symmetry
TA=np.array([0.2117,0.161,0.110,0.074]); TB=np.array([0.2095,0.160,0.112,0.072])
print("A/B asym %", np.round(100*abs(TA-TB)/((TA+TB)/2),2))
# TOST power at n=4, margin 5%, alpha=0.05, as function of true bias and sd
def tost_p(d,m):
    n=len(d); mu=d.mean(); se=d.std(ddof=1)/np.sqrt(n); df=n-1
    p1=1-stats.t.cdf((mu+m)/se,df); p2=stats.t.cdf((mu-m)/se,df); return max(p1,p2)
print("TOST power (n=4, margin 5%, alpha .05)")
for sd in [0.5,1.0,2.0]:
    row=[]
    for b in [0,1,2,3,3.5,4]:
        ps=[tost_p(rng.normal(b,sd,4),5) for _ in range(20000)]
        row.append(np.mean(np.array(ps)<0.05))
    print(" sd",sd,"bias 0,1,2,3,3.5,4 ->",np.round(row,2))
# slope CI behaviour for proportional bias with small scatter (synthetic -3% case from math_checks)
BF=(TA+TB)/2
for name,dp in [("bias-3",np.array([-3.1,-2.4,-3.8,-2.0])),("bias-4.3",np.array([-4.6,-3.9,-5.2,-3.5])),("bias0",np.array([0.5,-0.6,0.9,-0.8]))]:
    G=BF*(1+dp/100)
    r=stats.linregress(BF,G); tc=stats.t.ppf(0.975,2)
    lo,hi=r.slope-tc*r.stderr,r.slope+tc*r.stderr
    # trend of %diff vs angle
    ang=np.array([70,75,80,84]); r2=stats.linregress(ang,dp)
    print(name,"slope %.3f 95%%CI [%.3f,%.3f]  R2 %.4f | %%diff trend %.3f pp/deg, 95%%CI +/- %.3f, span over 14deg %.2f pp"%(r.slope,lo,hi,r.rvalue**2,r2.slope,tc*r2.stderr,14*r2.slope))
# Bland-Altman LoA CI (exact-ish, Bland-Altman 1999 approx) n=4
for sd in [1,2]:
    n=4; se_loa=np.sqrt(3*sd**2/n); print("BA n=4 sd",sd,"LoA half-width",1.96*sd,"CI on each LoA +/-",round(stats.t.ppf(.975,n-1)*se_loa,2))
# equivalence margin demonstrable at n=4 given sd: bias + t(.95,3)*sd/2 < m
for sd in [0.5,1,2]:
    print("sd",sd,"max |bias| for TOST pass at 5% margin:",round(5-stats.t.ppf(.95,3)*sd/2,2))
# BH on a family of m tests: per-test thresholds
m=6; print("BH crit for m=6, q=.05:",[round(0.05*i/m,4) for i in range(1,m+1)])
# GCI coarse grid elements across gap if global coarsening r=1.3 twice
for n in [3,4,5,6]:
    print("elements across gap fine",n,"-> after r^2=1.69:",round(n/1.69,2))
