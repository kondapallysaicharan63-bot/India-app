#!/usr/bin/env python3
"""
BME 298 hinge-gap FSI debate -- Round 1 NUMERICAL ANALYST checks.

Every number claimed in HANDOFF_BME298_2026-11-10.md that can be checked by arithmetic
is re-derived here. Each section Cn prints inputs, formula, result, and a verdict:
  PASS  = handoff claim reproduced (or consistent within stated precision)
  FLAG  = claim wrong, unverifiable, or the number has a consequence the handoff misses
  INFO  = new estimate (no handoff claim to compare against)

GLOBAL ASSUMPTIONS (stated once, used everywhere unless overridden):
  blood rho = 1060 kg/m3, mu = 3.5e-3 Pa s (Newtonian), nu = mu/rho
  hinge clearance h = 50 um (handoff Sec. 2)
  1 mmHg = 133.322 Pa
Anything marked ASSUMPTION is not in the handoff and is my choice.
"""
import math
import numpy as np
from scipy import stats, optimize

RHO, MU = 1060.0, 3.5e-3
NU = MU / RHO
H = 50e-6
MMHG = 133.322
DEG = math.pi / 180


def gci3(f1, f2, f3, r21, r32, Fs=1.25):
    """f1 = finest. Returns observed p, extrapolated value, GCI_fine (fraction), asymptotic ratio."""
    e21, e32 = f2 - f1, f3 - f2
    sgn = np.sign(e32 / e21)
    def g(p):
        q = math.log((r21 ** p - sgn) / (r32 ** p - sgn))
        return p - abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
    p = optimize.brentq(g, 0.05, 20)
    fext = (r21 ** p * f1 - f2) / (r21 ** p - 1)
    ea21 = abs((f1 - f2) / f1)
    gci21 = Fs * ea21 / (r21 ** p - 1)
    ea32 = abs((f2 - f3) / f2); gci32 = Fs * ea32 / (r32 ** p - 1)
    asym = gci32 / (r21 ** p * gci21)
    return p, fext, gci21, asym

def hdr(tag, title):
    print("\n" + "=" * 96)
    print(f"{tag}  {title}")
    print("=" * 96)


def verdict(ok, msg):
    print(f"  --> {'PASS' if ok else 'FLAG'}: {msg}")


def info(msg):
    print(f"  --> INFO: {msg}")


print(f"Assumptions: rho={RHO} kg/m3, mu={MU} Pa s, nu={NU:.3e} m2/s, h_gap={H*1e6:.0f} um")

# ---------------------------------------------------------------------------------------------
hdr("C1", "Penalization: Brinkman layer sqrt(K), penalty time scale tau = rho K / mu, ratio to gap")
Ks = [1e-11, 1e-10, 3e-10, 1e-9]
print(f"  {'K [m2]':>8} {'sqrt(K) [um]':>13} {'sqrtK/h_gap':>12} {'tau [s]':>10}")
for K in Ks:
    print(f"  {K:8.0e} {math.sqrt(K)*1e6:13.2f} {math.sqrt(K)/H:12.3f} {RHO*K/MU:10.2e}")
verdict(abs(RHO * 1e-11 / MU - 3e-6) < 0.2e-6, "tau(K=1e-11) = 3.0 us matches handoff '~3 us'; "
        "observed frozen step 8.5e-6 s = %.1f tau (step locked to penalty scale)." % (8.5e-6 / (RHO * 1e-11 / MU)))
verdict(abs(math.sqrt(1e-9) * 1e6 - 32) < 0.5, "sqrt(1e-9) = 31.6 um, i.e. 0.63 of the 50 um gap "
        "(handoff 'comparable to gap' is correct -- and that is the problem, see C3).")

# ---------------------------------------------------------------------------------------------
hdr("C2", "Penalization: leaflet thickness from CAD volume, Darcy leakage through the ghost leaflet")
V_leaf = 208.5e-9  # m3, handoff
# Planform area is not in the handoff. Bound it two ways:
#  (a) user range t = 0.6-0.8 mm -> A = V/t
#  (b) ASSUMPTION: leaflet planform ~ half-disk of radius R = 10.8 mm (ear end faces at z=+/-10.8 mm)
for t in [0.6e-3, 0.7e-3, 0.8e-3]:
    print(f"  t = {t*1e3:.1f} mm -> planform A = V/t = {V_leaf/t*1e6:6.1f} mm2")
A_half = 0.5 * math.pi * 10.8e-3 ** 2
print(f"  half-disk R=10.8 mm: A = {A_half*1e6:.1f} mm2 -> t = V/A = {V_leaf/A_half*1e3:.2f} mm")
info("CAD volume + half-disk planform implies t ~1.1 mm; 0.6-0.8 mm implies A = 260-350 mm2 (> half-disk). "
     "Use t = 0.7 mm as nominal, 1.1 mm as upper bound. Leakage scales as 1/t.")

# Darcy velocity through a penalized slab thick compared with sqrt(K): u = (K/mu) * dp / t
print("\n  Darcy leakage u = K dp / (mu t), t = 0.7 mm (in parentheses: t = 1.1 mm)")
cases = [("leaflet-normal dp at peak systole (from torque, see below)", 140.0),
         ("systolic transvalvular 5 mmHg", 5 * MMHG), ("systolic transvalvular 10 mmHg", 10 * MMHG),
         ("diastole/closed 80 mmHg", 80 * MMHG), ("diastole/closed 100 mmHg", 100 * MMHG)]
A_nom = V_leaf / 0.7e-3
print(f"  {'case':>55} {'dp[Pa]':>7} " + " ".join(f"{'K='+format(K,'.0e'):>17}" for K in Ks))
for name, dp in cases:
    row = []
    for K in Ks:
        u1 = K * dp / (MU * 0.7e-3)
        u2 = K * dp / (MU * 1.1e-3)
        row.append(f"{u1:7.3g}({u2:6.3g})")
    print(f"  {name:>55} {dp:7.0f} " + " ".join(f"{r:>17}" for r in row) + "  m/s")
# leaflet-normal dp estimated from body-fitted torque: T ~ dp_n * A * r_c (ASSUMPTION r_c = 5 mm lever arm)
T70 = 0.2117e-3
dpn = T70 / (A_nom * 5e-3)
print(f"\n  Leaflet-normal dp at 70 deg estimated from T = dp*A*r_c (A={A_nom*1e6:.0f} mm2, r_c=5 mm): {dpn:.0f} Pa")
for K in Ks:
    Q = K * 100 * MMHG / (MU * 0.7e-3) * A_nom * 2  # both leaflets
    print(f"  closed valve, 100 mmHg, both leaflets, K={K:.0e}: leak Q = {Q*6e4:8.2f} L/min")
Kreq = 0.01 * 0.7e-3 * MU / (100 * MMHG)  # K for <=1 cm/s Darcy leak at 100 mmHg
print(f"  K needed for Darcy leak <= 1 cm/s at 100 mmHg (t=0.7 mm): K <= {Kreq:.1e} m2 "
      f"-> tau = {RHO*Kreq/MU:.1e} s, sqrt(K) = {math.sqrt(Kreq)*1e6:.2f} um")
verdict(False, "K=1e-9 is acceptable only for the OPEN leaflet (leak ~6 cm/s at ~1 mmHg normal load, a few % "
        "of bulk velocity -- consistent with the -8.6% torque deficit). For closure/diastole (80-100 mmHg) "
        "the ghost leaflet leaks 4-6 m/s (K=1e-9) or 0.04-0.05 m/s even at K=1e-11 (~2 L/min through both leaflets): a closed ghost valve is not closed; K<=~2e-12 needed (tau ~0.6 us). "
        "Diastolic hinge leakage -- the clinically relevant hinge flow -- cannot be done with O1/O2 at feasible K.")

# ---------------------------------------------------------------------------------------------
hdr("C3", "Penalization next to the 50 um gap: 1D slot flux with one Brinkman wall (exact solution)")
# Slot 0<y<h fluid, no-slip at y=0 (housing, body-fitted); y>h penalized leaflet (Brinkman, semi-infinite).
# mu u'' = -G in fluid;  mu u'' - (mu/K) u = -G in porous. Continuity of u, u' at y=h.
# Gap flux counted as fluid part + decaying (non-Darcy) Brinkman part; far-field Darcy part is C2's leak.
def slot_flux_ratio(K, h=H):
    s = math.sqrt(K)
    # unknowns A (fluid u = -G y^2/(2mu) + A y) and C (porous u = K G/mu + C exp(-(y-h)/s)); G=mu=1
    # u: -h^2/2 + A h = K + C ;  u': -h + A = -C/s
    M = np.array([[h, -1.0], [1.0, 1.0 / s]])
    b = np.array([K + h * h / 2, h])
    A, C = np.linalg.solve(M, b)
    Qf = -h ** 3 / 6 + A * h ** 2 / 2
    Qp = C * s  # integral of decaying part
    Q0 = h ** 3 / 12  # both walls no-slip (real hinge)
    return (Qf + Qp) / Q0, A - h  # ratio, wall shear on the leaflet side proxy
print(f"  {'K':>8} {'sqrtK/h':>8} {'Q/Q_true':>9} {'WSS_housing/true':>17}")
for K in [1e-12] + Ks:
    r, _ = slot_flux_ratio(K)
    s = math.sqrt(K); h = H
    A = np.linalg.solve(np.array([[h, -1.0], [1.0, 1.0 / s]]), np.array([K + h*h/2, h]))[0]
    print(f"  {K:8.0e} {math.sqrt(K)/H:8.3f} {r:9.3f} {A/(H/2):17.3f}")
verdict(False, "with K=1e-9 the penalized leaflet acts as a slip wall (slip length ~sqrt(K)=32 um): the gap "
        "flux is ~3.9x and housing WSS ~1.9x the no-slip value (K=1e-10: 1.9x / 1.2x; K=1e-11: 1.2x / 1.07x). "
        "Need K<=1e-12 for ~6% flux error. "
        "=> O2 (ghost alone) cannot resolve the hinge gap; O1 must take gap flow from the film, not the 3D ghost.")

# ---------------------------------------------------------------------------------------------
hdr("C4", "K-study plan: error ~ sqrt(K)? Richardson in sqrt(K), 3 points, runtime extrapolation")
s = np.sqrt(np.array([1e-9, 3e-10, 1e-10]))
print(f"  sqrt(K) = {s*1e6} um; ratios {s[0]/s[1]:.3f}, {s[1]/s[2]:.3f} (non-constant -> use generalized formula)")
# Two-point Richardson with assumed order p in sqrt(K): T0 = T2 + (T2-T1)/((s1/s2)^p - 1)
def rich2(T1, T2, s1, s2, p=1.0):
    return T2 + (T2 - T1) / ((s1 / s2) ** p - 1)
# Three-point observed order with non-uniform ratio (Celik-type fixed point)
def observed_order(f, s):
    """f, s ordered coarse->fine in sqrt(K). Uses the Celik/Roache 3-grid formula (gci3), sqrt(K) as 'h'."""
    return gci3(f[2], f[1], f[0], s[1] / s[2], s[0] / s[1])[0]
# Synthetic demo: T(s) = T0*(1 - a s) with a chosen so that K=1e-9 gives -5% (halo removed, ASSUMPTION)
T0 = 0.2117; a = 0.05 / s[0]
f = T0 * (1 - a * s) + np.array([0, 0.0005, -0.0003]) * 0.0  # exact model
print(f"  synthetic T(K) = {f}  (T0 = {T0}, -5% at 1e-9)")
print(f"  observed order (exact data) p = {observed_order(f, s):.3f}; 2-pt Richardson p=1 on finest pair: "
      f"{rich2(f[1], f[2], s[1], s[2]):.4f}")
noisy = f + np.array([0.0, 0.0008, -0.0008])  # +/-0.4% noise (mesh/time-step) -- ASSUMPTION
print(f"  with +/-0.4% noise: T = {noisy}; observed p = {observed_order(noisy, s):.2f}; "
      f"p=1 extrapolation = {rich2(noisy[1], noisy[2], s[1], s[2]):.4f}")
verdict(False, "3 points determine (T0, c, p) exactly -> zero residual, no check of the asymptotic range. "
        "Fix p=1 (theory for Brinkman velocity error ~ sqrt(K): needs literature check) and use the 3rd point "
        "as a consistency check. Noise of +/-0.4% already moves observed p from 1.0 to ~2-3.")
# Mesh vs Brinkman layer: mean element edge
Vtot = 54249.7e-9; Ne = 2050445
v_el = Vtot / Ne; a_el = (6 * math.sqrt(2) * v_el) ** (1 / 3)
print(f"  mean element volume {v_el*1e9:.4f} mm3 -> regular-tet edge {a_el*1e3:.3f} mm")
verdict(False, f"sqrt(K) (10-32 um) << local element size (tens-hundreds of um): the Brinkman layer is "
        "unresolved, interface error is O(h_mesh) not O(sqrt(K)); a K study on a fixed mesh will plateau. "
        "K study must be coupled with local refinement (h_mesh <~ sqrt(K)) or reported as K-sensitivity only.")
# runtime scaling from the two runs
t_1e9 = 1 + 20 / 60  # h for full Study 2
t_1e11 = 10 * 0.21 / 0.0867  # h, ASSUMPTION Study 2 end time 0.21 s, constant (frozen) rate
alpha = math.log(t_1e11 / t_1e9) / math.log(100)
print(f"  wall time: K=1e-9 {t_1e9:.2f} h; K=1e-11 extrapolated {t_1e11:.1f} h -> wall ~ K^-{alpha:.2f}")
for K in [3e-10, 1e-10]:
    print(f"    predicted K={K:.0e}: {t_1e9*(1e-9/K)**alpha:.1f} h (ASSUMPTION: power-law between the 2 runs)")
info("K=3e-10 and 1e-10 fit in one night each on the desktop.")

# ---------------------------------------------------------------------------------------------
hdr("C5", "Halo: VA = 299.3 / VB = 300.1 vs 208.5 mm3 -> equivalent halo thickness vs mesh size")
for name, Vg in [("A", 299.3), ("B", 300.1)]:
    dV = Vg - 208.5
    print(f"  {name}: excess {dV:.1f} mm3 = {dV/208.5*100:.1f}%")
dV = 90.8
for t in [0.7, 1.1]:
    Awet = 2 * 208.5 / t  # both faces, edges neglected
    print(f"  t={t} mm: wetted area ~{Awet:.0f} mm2 -> halo thickness {dV/Awet*1e3:.0f} um per face")
verdict(True, "+43.5%/+43.9% matches '+44%'. Halo ~150-240 um per face = ~3-5 gap widths, i.e. ~1/4-2/5 of a "
        f"mean element ({a_el*1e3:.2f} mm): consistent with nearest-point snapping over one element layer. "
        "Halo inflates thickness by ~44% -> plausible cause of FyA +25.6% (form drag on thicker plate). "
        "With a halo of 3-5 gap widths, the ghost leaflet in v22 FILLS the 50 um hinge gap: all K=1e-9 results "
        "near the hinge are meaningless until FixChi is run.")
print(f"  volume bookkeeping: 53832.7 + 2*208.5 = {53832.7+2*208.5:.1f} mm3 (handoff 54,249.7)")
dof_per_el = 1575704 / Ne
print(f"  DOF/element = {dof_per_el:.3f}; P1+P1 NS (4 DOF/node), tets ~5.5-6 el/node -> "
      f"{4/6:.2f}-{4/5.5:.2f} DOF/el")
verdict(abs(53832.7 + 417 - 54249.7) < 0.1, "volume and DOF counts are internally consistent.")

# ---------------------------------------------------------------------------------------------
hdr("C6", "Body-fitted torque vs angle: trend, fit, A/B symmetry")
ang = np.array([70, 75, 80, 84.0])
TA = np.array([0.2117, 0.161, 0.110, 0.074]); TB = np.array([0.2095, 0.160, 0.112, 0.072])
FyA = np.array([49.205, 33.4, 20.9, 13.1])
for nm, y in [("TzA", TA), ("TzB", TB), ("FyA", FyA)]:
    p1 = np.polyfit(ang, y, 1); r2 = 1 - np.sum((y - np.polyval(p1, ang))**2) / np.sum((y - y.mean())**2)
    pe = np.polyfit(ang, np.log(y), 1)
    print(f"  {nm}: linear slope {p1[0]:.4f}/deg R2={r2:.4f}, zero at {-p1[1]/p1[0]:.1f} deg; "
          f"exp decay {pe[0]*100:.2f}%/deg; monotonic={bool(np.all(np.diff(y)<0))}")
sym = np.abs(TA - TB) / ((TA + TB) / 2) * 100
print(f"  A/B asymmetry %: {dict(zip(ang.astype(int), np.round(sym,2)))}")
lo = np.abs(0.0735 - 0.0725) / 0.073 * 100; hi = (0.0745 - 0.0715) / 0.073 * 100
print(f"  84 deg with 2-sig-fig rounding: asymmetry could be {lo:.1f}%..{hi:.1f}%")
verdict(False, f"max A/B asymmetry is {sym.max():.1f}% at 84 deg, not 'within 2.2%' (2.2% fits 70-80 only); "
        "torques at 75-84 deg are reported to 2 sig figs, too coarse for a 2% symmetry or 5% bias gate. "
        "Report 4 sig figs. Torque decays ~8%/deg; extrapolated to 85 deg ~0.064-0.066 mN m (still opening).")

# ---------------------------------------------------------------------------------------------
hdr("C7", "Statistics with n = 4 angles x 2 leaflets: what Bland-Altman / regression / paired t can detect")
# (a) R2 gate: a pure proportional bias of -8.6% gives R2 = 1 -> gate non-discriminating
ghost = 0.914 * TA
r2g = stats.linregress(TA, ghost).rvalue ** 2
print(f"  ghost = 0.914 x body-fitted (8.6% bias everywhere): regression R2 = {r2g:.4f}, slope 0.914")
verdict(False, "R2 > 0.95 gate passes even with a large systematic bias; use slope/intercept with CI and "
        "max |%diff| instead.")
# (b) paired t power: minimum detectable mean %diff at alpha=0.05, power 0.8
def mde(n, sd, alpha=0.05, power=0.8):
    df = n - 1
    tc = stats.t.ppf(1 - alpha / 2, df)
    # two-sided power; lower tail is negligible and nct.cdf can return NaN for large nc, so use sf only
    f = lambda d: stats.nct.sf(tc, df, d / sd * math.sqrt(n)) - power
    return optimize.brentq(f, 1e-3 * sd, 6 * sd)
print(f"  {'n':>3} {'sd of %diff':>12} {'MDE (80% power)':>16} {'95% CI half-width':>18}")
for n in [4, 8]:
    for sd in [0.5, 1.0, 2.0]:
        hw = stats.t.ppf(0.975, n - 1) * sd / math.sqrt(n)
        print(f"  {n:3d} {sd:12.1f} {mde(n, sd):15.2f}% {hw:17.2f}%")
print("  n=8 treats A and B as independent; they are mirror images (r~1), so effective n ~ 4.")
# (c) Bland-Altman limits of agreement CI with n=4
for sd in [1.0, 2.0]:
    n = 4; se_loa = math.sqrt(3 * sd**2 / n); tq = stats.t.ppf(0.975, n - 1)
    print(f"  BA n=4, sd={sd}%: LoA = bias +/- {1.96*sd:.2f}%, each LoA uncertain by +/-{tq*se_loa:.2f}%")
verdict(False, "paired t answers 'is bias != 0' (with a deterministic model error it eventually always is), "
        "not 'is |bias| < 5%'. Use an equivalence test (TOST) or simply the CI of the bias vs the 5% margin; "
        "with n=4 and sd~1% the 95% CI half-width is ~1.6%, so bias <= ~3.4% is needed to demonstrate <5%. "
        "Better: ASME V&V-20-style comparison E = S-D vs u_val from GCI (needs literature check).")
# (d) Benjamini-Hochberg example
p = np.array([0.004, 0.019, 0.030, 0.21]); names = ["TzA bias", "TzB bias", "FyA bias", "slope!=1"]
m = len(p); order = np.argsort(p); crit = 0.05 * (np.arange(1, m + 1)) / m
passed = p[order] <= crit; kmax = np.max(np.where(passed)[0]) + 1 if passed.any() else 0
print("  BH (q=0.05) example with hypothetical p-values:")
for i, j in enumerate(order):
    print(f"    rank {i+1}: {names[j]:9s} p={p[j]:.3f} crit={crit[i]:.4f} -> {'reject' if i < kmax else 'keep'}")
info("BH is mechanically fine but with 3-4 tests on the same 4 angles it adds little; report it only if "
     "the committee asks for p-values.")

# ---------------------------------------------------------------------------------------------
hdr("C8", "Lubrication (Reynolds film) validity: eps, Re_gap, eps*Re, wall shear, unsteadiness")
L_vals = [1e-3, 1.5e-3, 2e-3]
print(f"  eps = h/L: " + ", ".join(f"L={L*1e3:.1f} mm -> {H/L:.4f}" for L in L_vals))
verdict(abs(H / 2e-3 - 0.025) < 1e-9, "eps = 0.025 corresponds to L = 2 mm.")
# Implied gap velocities behind the handoff's eps*Re
for lab, eR in [("systole", 0.03), ("closure", 0.24)]:
    Re = eR / 0.025; U = Re * MU / (RHO * H)
    print(f"  handoff {lab}: eps*Re={eR} -> Re={Re:.1f} -> U={U:.3f} m/s")
# Closure kinematics: 15 deg in T_c, ear radius 2 mm (ASSUMPTION T_c = 2-10 ms)
print("  closure Couette velocity of ear surface:")
for Tc in [2e-3, 5e-3, 10e-3]:
    om = 15 * DEG / Tc
    for r in [2e-3, 2.2e-3]:
        U = om * r
        print(f"    Tc={Tc*1e3:4.0f} ms omega={om:6.1f} rad/s r={r*1e3:.1f} mm U={U:.3f} m/s "
              f"Re={RHO*U*H/MU:.2f} Couette tau=mu U/h={MU*U/H:5.1f} Pa")
print("  pressure-driven (Poiseuille) part across the hinge:")
for dpm in [5, 10, 80, 100]:
    for L in L_vals:
        dpds = dpm * MMHG / L
        Um = H**2 * dpds / (12 * MU); Re = RHO * Um * H / MU
        print(f"    dp={dpm:3d} mmHg L={L*1e3:.1f} mm: dp/ds={dpds:.2e} Pa/m  tau_P=(h/2)dp/ds={H/2*dpds:6.1f} Pa"
              f"  Umean={Um:.3f} m/s Re={Re:5.2f} eps*Re={H/L*Re:.3f}")
verdict(True, "'~270 Pa at closure' is reproduced only as the POISEUILLE term with 80-100 mmHg across "
        "L ~1-1.5 mm (167-333 Pa); the Couette term from ear motion is only 5-40 Pa. Note: 80-100 mmHg is the "
        "CLOSED-valve (diastolic) load; during the closing swing the back pressure is a few mmHg, so 270 Pa "
        "belongs to diastolic leakage, not to the closing motion itself (except water-hammer transients).")
verdict(True, "eps*Re = 0.01-0.08 in systole and 0.12-0.60 in diastolic leakage (L=2..1 mm): film valid in systole, "
        "marginal (O(eps*Re) ~10-60% inertia correction) for 80-100 mmHg leakage, esp. at mouths.")
# Unsteady (squeeze/Womersley) parameter h^2 omega / nu
print("  unsteady film parameter h^2*omega/nu:")
for lab, om in [("cardiac 72 bpm", 2*math.pi*1.2), ("closure 10 ms", 1/10e-3), ("closure 2 ms", 1/2e-3),
                ("impact/water-hammer 0.2 ms", 1/0.2e-3)]:
    print(f"    {lab:28s} omega={om:8.1f} 1/s -> {H**2*om/NU:.4f}")
info("Film is quasi-steady (<0.2) except at leaflet impact (sub-ms), where h^2*omega/nu ~ O(1): "
     "the Reynolds film needs a transient (dh/dt squeeze) term and is unreliable at impact.")

# ---------------------------------------------------------------------------------------------
hdr("C9", "ALE: hinge shear gamma = R dtheta / h_gap; end-face sliding; inversion angle; tip displacement")
R_needed = 0.67 * H / DEG
print(f"  gamma = 0.67/deg requires R = 0.67*h/(1 deg) = {R_needed*1e3:.2f} mm")
verdict(abs(R_needed - 1.9e-3) < 0.1e-3, "R = 1.92 mm = ear radius ~2 mm (also |pivot x| = 1.90 mm). "
        "At the observed inversion site r = 2.1-2.3 mm gamma = %.2f-%.2f per deg." % (2.1e-3*DEG/H, 2.3e-3*DEG/H))
for r in [2.1e-3, 2.3e-3]:
    print(f"  r={r*1e3:.1f} mm: end-face slide {r*DEG*1e6:.1f} um/deg = {r*DEG/H:.2f} gap widths/deg")
print("  angle at which a gap-spanning element of height h/n is sheared to gamma = 1 (45 deg skew),")
print("  if the whole slide is taken by that element (no smoothing):")
for n in [1, 2, 3, 4]:
    for r in [2.1e-3, 2.3e-3]:
        th = (H / n) / r / DEG
        print(f"    n={n} r={r*1e3:.1f} mm -> {th:.2f} deg")
print(f"  observed history: Laplace 1.1, slip 2.3, Yeoh 3.1 deg; gamma at 3 deg, r=2.2 mm = "
      f"{2.2e-3*3*DEG/H:.2f}")
verdict(True, "Laplace failure at 1.1 deg ~ gamma=1 over the full gap (h/r = 1.25-1.36 deg); Yeoh lasting to "
        "~3 deg means the smoother spread the slide over ~2.3 gap widths. Note: uniform simple shear has "
        "det F = 1 and never inverts; inversion is created where a SLIDING face node meets a PINNED corner "
        "(du/dx < -1 inside one element). Remeshing at distortion 1.5 cannot cure this unless the remesh "
        "interval is < ~(h/n)/r = 0.3-0.65 deg at the end faces, i.e. >20-45 remeshes over 14 deg.")
tip = 10.8e-3 * 14 * DEG; chord = 2 * 10.8e-3 * math.sin(7 * DEG)
print(f"  tip displacement 10.8 mm x 14 deg: arc {tip*1e3:.3f} mm, chord {chord*1e3:.3f} mm")
verdict(True, "2.6 mm reproduced. FLAG-lite: 10.8 mm is the ear END-FACE z-coordinate (span direction); the "
        "in-plane pivot-to-tip radius must be read from the CAD (rotation is about z). If tip radius differs, "
        "expected displacement scales linearly.")

# ---------------------------------------------------------------------------------------------
hdr("C10", "ALE run timing: sim-seconds per wall hour and finish-time extrapolation")
import datetime as dt
clk = [dt.datetime(2026,11,8,20,53), dt.datetime(2026,11,8,22,46), dt.datetime(2026,11,9,1,38),
       dt.datetime(2026,11,9,5,59), dt.datetime(2026,11,9,7,54), dt.datetime(2026,11,9,9,58),
       dt.datetime(2026,11,9,14,54), dt.datetime(2026,11,9,20,34), dt.datetime(2026,11,9,22,7)]
ts = [0, 0.066, 0.118, 0.2175, 0.320, 0.384, 0.502, 0.603, 0.654]
angs = [0, 0.9, 1.65, 3.0, 4.5, 5.4, 7.0, 8.4, 9.2]
rates = []
for i in range(1, len(ts)):
    hrs = (clk[i] - clk[i-1]).total_seconds() / 3600; rate = (ts[i] - ts[i-1]) / hrs; rates.append(rate)
    print(f"  interval {i}: {hrs:5.2f} h, dt={ts[i]-ts[i-1]:.4f} s -> {rate:.4f} s/h "
          f"({rate*14:.2f} deg/h); angle check 14*t={14*ts[i]:.2f} vs table {angs[i]}")
tot_h = (clk[-1] - clk[0]).total_seconds() / 3600
rem = 1 - ts[-1]
for lab, rate in [("overall mean", ts[-1] / tot_h), ("last interval", rates[-1]),
                  ("slowest (out-of-core) interval", min(rates)), ("mean of last 3", np.mean(rates[-3:]))]:
    fin = clk[-1] + dt.timedelta(hours=rem / rate)
    print(f"  {lab:32s} {rate:.4f} s/h -> remaining {rem/rate:5.1f} h -> finish {fin:%b %d %H:%M}")
verdict(True, "65.4% done (handoff 66%); finish Nov 10 ~08:00-17:30, consistent with 'midday-evening'. "
        "The run has likely finished (or stalled) by the time this is read: check the Log.")

# ---------------------------------------------------------------------------------------------
hdr("C11", "GCI (Roache, 3 meshes, Fs = 1.25) -- function + worked example; 3D DOF/memory growth")
# hypothetical torques on fine/medium/coarse (ASSUMPTION, illustration only)
f1, f2, f3 = 0.2117, 0.2152, 0.2231
for r in [1.3, 1.5]:
    p, fe, g, asr = gci3(f1, f2, f3, r, r)
    print(f"  example f=({f1},{f2},{f3}), r={r}: p={p:.2f}, f_ext={fe:.4f}, GCI_fine={g*100:.2f}%, "
          f"asymptotic ratio={asr:.3f} (should be ~1)")
info("3D refinement ratio from element counts: r = (N_fine/N_coarse)^(1/3); Celik et al. recommend r >= 1.3 "
     "(needs literature check).")
N0 = 1.58e6
print("  refining FROM the current 1.58 M DOF mesh (memory M = M0 (N/N0)^b, b ASSUMED 1.3-1.5):")
for M0, lab in [(28.7, "M0=28.7 GB (current run fills the 28.7 GB box)"), (15.0, "M0=15 GB (optimistic)")]:
    print(f"   {lab}")
    for r in [1.3, 1.5, 2.0]:
        N = N0 * r ** 3
        print(f"    r={r}: DOF {N/1e6:6.2f} M  memory {M0*(N/N0)**1.3:6.0f}-{M0*(N/N0)**1.5:6.0f} GB")
print("  COARSENING from 1.58 M instead (current mesh = fine grid):")
for r in [1.3, 1.5]:
    print(f"    r={r}: medium {N0/r**3/1e6:.2f} M DOF, coarse {N0/r**6/1e6:.2f} M DOF")
verdict(False, "any refined grid (r>=1.3 -> 3.5 M DOF, 2.2x) exceeds 28-32 GB with a direct solver under "
        "either calibration; the GCI must use the current mesh as the FINEST grid (coarsen r=1.3: 0.72 M, "
        "0.33 M DOF) or run on HPC / with an iterative solver. But the coarse grid must still put >=2-3 "
        "elements across the 50 um gap, so global coarsening by r^2=1.69 may break the gap: refine/coarsen the "
        "bulk and the gap separately and report two GCIs.")

# ---------------------------------------------------------------------------------------------
hdr("C12", "Added-mass ratio for two-way FSI (leaflet rotation)")
I_s = 9.5e-9
# Added rotational inertia of a 2D flat plate (chord c) rotating about one edge, per unit span:
#   I_a' = (9 pi / 128) rho c^4   (= pi rho b^4/8 about mid-chord + pi rho b^2 * b^2, b = c/2)
# times effective span s_eff. ASSUMPTION: c = pivot-to-tip 8-11 mm, s_eff = 0.6-0.8 x 21.6 mm (D-shape).
print(f"  {'c [mm]':>7} {'s_eff [mm]':>10} {'I_a [kg m2]':>12} {'I_a/I_s':>8}")
for c in [8e-3, 10e-3, 11e-3]:
    for se in [0.6 * 21.6e-3, 0.8 * 21.6e-3]:
        Ia = 9 * math.pi / 128 * RHO * c ** 4 * se
        print(f"  {c*1e3:7.1f} {se*1e3:10.1f} {Ia:12.2e} {Ia/I_s:8.2f}")
m_leaf = 208.5e-9 * 1.9e3  # ASSUMPTION pyrolytic carbon ~1.9 g/cm3 (needs check of CAD material)
print(f"  consistency of I_s: m = {m_leaf*1e3:.3f} g (rho_s 1.9e3) -> radius of gyration "
      f"{math.sqrt(I_s/m_leaf)*1e3:.2f} mm (plausible for pivot inside the leaflet)")
verdict(False, "I_added/I_leaflet ~ 1.3-6 (>1; confinement in the housing raises it further): a plain "
        "partitioned (Dirichlet-Neumann) coupling of torque and angle is unstable without strong relaxation "
        "(Aitken) or an artificial added inertia ~I_a. COMSOL's fully-coupled rigid-rotation ODE in the same "
        "Newton system is monolithic and immune -- but in the ghost model the torque depends on chi via General "
        "Extrusion with NO d(chi)/d(theta) Jacobian, so O1/O2 two-way FSI is effectively partitioned w.r.t. "
        "the geometric coupling and will need the Phase E stabilization. Also: with a moving ghost leaflet the "
        "penalty integral misses the momentum of fluid inside chi; add -d/dt int(rho chi r x u) (Uhlmann-type "
        "correction, needs literature check) or the torque is biased by ~rho*V_leaf*r^2*alpha.")
# size of internal-fluid inertia term vs leaflet inertia
I_int = RHO * 208.5e-9 * (math.sqrt(I_s / m_leaf)) ** 2
print(f"  internal-fluid inertia rho_f V r_g^2 = {I_int:.2e} kg m2 = {I_int/I_s*100:.0f}% of I_s")

# ---------------------------------------------------------------------------------------------
hdr("C13", "Other numbers in the handoff")
dT_A = (0.1934 - 0.2117) / 0.2117 * 100; dT_B = (0.1912 - 0.2095) / 0.2095 * 100
dF = (61.81 - 49.205) / 49.205 * 100
print(f"  ghost vs BF: TzA {dT_A:.2f}% (claim -8.6), TzB {dT_B:.2f}% (claim -8.8), FyA {dF:.2f}% (claim +25.6)")
verdict(abs(dT_B + 8.8) < 0.1, "TzB is -8.74% (handoff -8.8%, rounding).")
print(f"  ghost A/B symmetry: {abs(0.1934-0.1912)/0.1923*100:.2f}% (claim 1%)")
print(f"  ALE rotation 14 deg/s x 1 s = 14 deg -> 70+14 = 84 deg: never reaches the 85 deg stop (by design?)")
# sqrt(K) at 3e-10 vs element thickness of the leaflet: number of Brinkman lengths across t
for K in Ks:
    print(f"  K={K:.0e}: leaflet thickness/sqrt(K) = {0.7e-3/math.sqrt(K):.0f} (Darcy regime if >>1)")
info("Torque expression TzA_pen uses (x+pin_x) but (y-pin_y): consistent only if pin_x stores +|x_A| "
     "(1.90) and pin_y stores the signed -6.85. Verify the B-leaflet expression uses (x-pin_x). UNVERIFIED "
     "(script not available here).")
print("\nEnd of Round-1 checks.")


# #############################################################################################
# ROUND 2 sections (R2-C1 ...). Appended Nov 10. Respond to Round-1 openings of other agents.
# #############################################################################################
print("\n\n" + "#" * 96 + "\n# ROUND 2 CHECKS\n" + "#" * 96)

# ---------------------------------------------------------------------------------------------
hdr("R2-C1", "Leakage conflict: Challenger 29.7x vs Numerical 3.9x -- same physics, different control volumes")
def brinkman_layer(K, T, bc, h=H):
    """1D composite: fluid 0<y<h (no-slip at y=0, housing), penalized layer h<y<h+T.
    Far side y=h+T: bc='noslip' (layer backed by a wall: Challenger setup) or 'sym' (du/dy=0).
    T=inf -> semi-infinite layer (Numerical R1 setup). Same dp/dx in fluid and layer.
    Returns flux ratios to the true no-slip slot Q0 = G h^3/(12 mu):
      Qf  = flux in the fluid gap itself
      Qbl = flux in the decaying Brinkman boundary layer (thickness ~sqrt(K)) inside the 'solid'
      Qd  = Darcy core flux K*T*G/mu inside the 'solid' (grows linearly with T)"""
    sK = math.sqrt(K)
    if T == math.inf:
        A, C = np.linalg.solve(np.array([[h, -1.0], [1.0, 1 / sK]]), np.array([K + h * h / 2, h]))
        Qf = -h ** 3 / 6 + A * h * h / 2
        return Qf / (h ** 3 / 12), C * sK / (h ** 3 / 12), math.inf
    E = math.exp(-T / sK)
    M = np.zeros((3, 3)); b = np.zeros(3)
    M[0] = [h, -1, -E]; b[0] = K + h * h / 2           # u continuity at y=h
    M[1] = [1, 1 / sK, -E / sK]; b[1] = h              # u' continuity at y=h
    if bc == "noslip":
        M[2] = [0, E, 1]; b[2] = -K                    # u(h+T)=0
    else:
        M[2] = [0, -E / sK, 1 / sK]; b[2] = 0          # u'(h+T)=0
    A, C1, C2 = np.linalg.solve(M, b)
    Q0 = h ** 3 / 12
    Qf = -h ** 3 / 6 + A * h * h / 2
    Qbl = (C1 + C2) * sK * (1 - E)
    return Qf / Q0, Qbl / Q0, K * T / Q0
print(f"  {'setup':>38} {'K':>7} {'Q_gap':>7} {'Q_BL':>7} {'Q_Darcy':>8} {'TOTAL':>8}")
for lab, T, bc in [("Numerical R1: semi-inf, gap+BL only", math.inf, None),
                   ("Challenger: T=300um, wall-backed", 300e-6, "noslip"),
                   ("T=300um, symmetric far side", 300e-6, "sym"),
                   ("T=700um (leaflet t), wall-backed", 700e-6, "noslip")]:
    for K in [1e-9, 1e-10, 1e-11]:
        qf, qbl, qd = brinkman_layer(K, T, bc)
        tot = qf + qbl + (0 if qd == math.inf else qd)
        qd_s = "excl." if qd == math.inf else f"{qd:8.2f}"
        tot_s = f"{qf+qbl:8.3f}*" if qd == math.inf else f"{tot:8.2f}"
        print(f"  {lab:>38} {K:7.0e} {qf:7.3f} {qbl:7.3f} {qd_s:>8} {tot_s:>8}")
print("  (* = gap + Brinkman layer only; the Darcy core is excluded by construction)")
qf, qbl, qd = brinkman_layer(1e-9, 300e-6, "noslip")
verdict(abs(qf + qbl + qd - 29.7) < 0.1 and abs(qf - 3.63) < 0.01,
        f"Challenger's 3.63x / 26.1x / 29.7x reproduced exactly (gap {qf:.2f}x; solid = BL {qbl:.2f} + Darcy {qd:.1f} "
        f"= {qbl+qd:.1f}x). Numerical R1 3.93x = gap 3.63x + decaying layer 0.29x of a semi-infinite layer.")
print("  Darcy core / true gap flux = 12 K T / h^3:")
for K in [1e-9, 1e-10, 1e-11, 1e-12]:
    print(f"    K={K:.0e}: T=0.3 mm {12*K*0.3e-3/H**3:7.2f}x   T=0.7 mm {12*K*0.7e-3/H**3:7.2f}x")
Kmax = 0.05 * H ** 3 / (12 * 0.3e-3)
print(f"  K for Darcy-core leak < 5% of gap flux: T=0.3 mm K <= {Kmax:.1e}; T=0.6 mm K <= {Kmax/2:.1e} "
      "(Challenger 8.7e-13 = the T=0.6 mm case: PASS)")
verdict(True, "NO CONFLICT. Three different leak paths: (a) slip-enhanced flux in the gap itself (3.6x at 1e-9; "
        "controls gap velocity and housing WSS 1.9x) -- Numerical's number; (b) gap-parallel Darcy flow through the "
        "penalized ear body (26x for a 0.3 mm ear, scales with ear size) -- Challenger's extra; (c) through-thickness "
        "Darcy flow across the leaflet driven by the transvalvular dp (C2: 4.4-5.4 m/s, ~1.3 L/s per leaflet at "
        "80 mmHg) -- both agree (Challenger 1298 mL/s = C2 with t=0.7 mm). For HINGE leakage the relevant total is "
        "(a)+(b) = Challenger's ~30x; for gap WSS it is (a); for valve closure it is (c).")
u_c = 1e-9 * 80 * MMHG / (MU * 0.7e-3)
print(f"  check (c): u = {u_c:.2f} m/s x 297.9 mm2 = {u_c*297.9e-6*1e6:.0f} mL/s per leaflet at 80 mmHg")

# ---------------------------------------------------------------------------------------------
hdr("R2-C2", "K-study runtime conflict: Feasibility (6.5 h / 19 h) vs Numerical R1 (2.8 h / 5.7 h)")
# Data: K=1e-9 Study 2: 80 min wall. K=1e-11: 10 h, frozen dt=8.5e-6 s, 'reached t=0.0867'.
# Numerical R1 ASSUMED: whole 0..0.0867 s at constant rate, full window 0.21 s -> 24.2 h, power law between 2 totals.
# Feasibility ASSUMED: window starts at 0.08 s (0.13 s simulated), so the 1e-11 run advanced 0.0067 s = 788 steps
#   -> 45.7 s/step; dt = min(dt_1e-9 mean, 2.8 tau).
steps_1e11 = 0.0067 / 8.5e-6; sps = 10 * 3600 / steps_1e11
steps_1e9 = 80 * 60 / sps; dt_mean_1e9 = 0.13 / steps_1e9
print(f"  Feasibility reading: K=1e-11 did {steps_1e11:.0f} steps -> {sps:.1f} s/step; "
      f"K=1e-9 80 min -> {steps_1e9:.0f} steps, mean dt {dt_mean_1e9:.2e} s")
alt_sps = 10 * 3600 / (0.0867 / 8.5e-6)
print(f"  Numerical-R1 reading implies {0.0867/8.5e-6:.0f} steps in 10 h = {alt_sps:.1f} s/step at 1.58 M DOF MUMPS "
      "-- implausibly cheap (a single MUMPS factorization at this size is tens of seconds).")
verdict(False, "my R1 K^-0.63 law extrapolated from the TOTAL wall time of the 1e-11 run assuming the step was frozen "
        "from t=0; Feasibility's per-step model is self-consistent (both runs give ~46 s/step) and I adopt it.")
print(f"  {'K':>7} {'tau':>9} {'dt=min(dt9,2.8tau)':>19} {'0.13 s window [h]':>18} {'0.025 s restart [h]':>20} "
      f"{'R1 power law [h]':>17}")
rec = {}
for K in [1e-9, 3e-10, 1e-10, 3e-11, 1e-11]:
    tau = RHO * K / MU; dtk = min(dt_mean_1e9, 2.8 * tau)
    full = 0.13 / dtk * sps / 3600; win = 0.025 / dtk * sps / 3600
    pl = (80 / 60) * (1e-9 / K) ** 0.63
    rec[K] = (full, win)
    print(f"  {K:7.0e} {tau:9.2e} {dtk:19.2e} {full:18.1f} {win:20.2f} {pl:17.1f}")
# range: lower bound if the lock factor grows (dt/tau=5 instead of 2.8), upper if per-step cost +50% (more Newton its)
for K in [3e-10, 1e-10]:
    tau = RHO * K / MU
    lo = 0.13 / min(dt_mean_1e9, 5 * tau) * sps / 3600; hi = 0.13 / min(dt_mean_1e9, 2.8 * tau) * sps * 1.5 / 3600
    print(f"  RECOMMENDED K={K:.0e}: full window {rec[K][0]:.1f} h (range {lo:.1f}-{hi:.1f} h); "
          f"0.025 s restart window {rec[K][1]:.1f} h (range {rec[K][1]*lo/rec[K][0]:.1f}-{rec[K][1]*hi/rec[K][0]:.1f} h)")
print(f"  model check at K=1e-9: predicts {rec[1e-9][0]:.1f} h vs measured 1.33 h -> lock factor dt/tau there is "
      f"{dt_mean_1e9/(RHO*1e-9/MU):.1f}, not 2.8: the lock factor is not constant, hence the lower bound uses 5 tau.")
info("Use 0.025 s restart windows around peak (t=0.11-0.135 s): ~1.2 h at 3e-10 and ~3.7 h at 1e-10; "
     "all 3 K x 2 angles fit in one night. Full-window 1e-10 (~19 h) only if needed.")

# ---------------------------------------------------------------------------------------------
hdr("R2-C3", "Inventor: frozen-angle transpiration flux error ~3 e dtheta/h; back-coupling; leakage/beat")
# Reynolds: q_p = -(h^3/12mu) grad p, q_c = h U/2. A rotation axis offset e from the arc centre moves the arc by
# e*dtheta, so the local gap changes by dh = e*dtheta*cos(phi) (|dh| <= e dtheta). Frozen geometry keeps h_k:
#   dQ_p/Q_p = (1 + dh/h)^3 - 1 ~ 3 dh/h ;  dQ_c/Q_c = dh/h
print(f"  {'e [um]':>7} " + " ".join(f"{'dth='+str(d)+'deg':>22}" for d in [0.25, 0.5, 1.0]) + f" {'dth for 5%':>11}")
for e in [10e-6, 25e-6, 50e-6]:
    row = []
    for d in [0.25, 0.5, 1.0]:
        dh = e * d * DEG
        row.append(f"lin {3*dh/H*100:5.2f}% exact {((1+dh/H)**3-1)*100:5.2f}%")
    dmax = 0.05 * H / (3 * e) / DEG
    print(f"  {e*1e6:7.0f} " + " ".join(f"{r:>22}" for r in row) + f" {dmax:9.2f} deg")
verdict(True, "Inventor's 3 e dtheta/h table reproduced (e=10: 0.26/0.52/1.05%; e=50: 1.31/2.62/5.24%; 5% limits 4.8 "
        "and 0.95 deg). It is a LOCAL upper bound: dh ~ e dtheta cos(phi) changes sign along the arc, so the net flux "
        "error over a full arc is second order unless the arc is short. Bound also applies to the WSS (Poiseuille part "
        "~h: error e dtheta/h, one third of the flux error).")
# back-coupling: viscous friction torque T_f = mu*Omega*R^2*A/h  (Couette on ear area A at radius R)
R_ear = 2.2e-3; A_gap = 6e-3 * 1e-3  # Inventor UNVERIFIED w=6 mm x L=1 mm per hinge
D = MU * R_ear ** 2 * A_gap / H
print(f"  D = mu R^2 A / h = {D:.2e} N m s per hinge (Inventor 2.2e-9); x2 for two ears per leaflet")
for lab, Om in [("mean closing 15deg/30ms", 15 * DEG / 0.03), ("peak 100 rad/s", 100.0)]:
    for nh in [1, 2]:
        print(f"    {lab:24s} hinges={nh}: T_f={nh*D*Om:.2e} N m -> T_f/0.2 mN m = {nh*D*Om/0.2e-3:.1e}")
verdict(True, "ratio 1e-4 (mean) to 1.1e-3 (peak) reproduced per hinge; 2e-4-2.2e-3 with both ears. Still negligible. "
        "Caveat: 0.2 mN m is the systolic 70-deg torque; the closing (reverse-flow) torque is not in the handoff "
        "(UNVERIFIED) -- if larger, the ratio is even smaller. Pressure on concentric arcs/z-normal end faces gives "
        "zero torque about z, so friction is the only hinge torque if e~0.")
def slot_Q(w, L, dp, h=H):
    return w * h ** 3 * dp / (12 * MU * L)
q_inv = 4 * slot_Q(6e-3, 1e-3, 80 * MMHG); q_cli = 4 * slot_Q(4e-3, 2e-3, 100 * MMHG)
print(f"  leakage, 4 hinges: Inventor (w=6mm,L=1mm,80mmHg) {q_inv*1e6:.2f} mL/s x 0.5 s = {q_inv*0.5e6:.2f} mL/beat; "
      f"Clinical (w=4mm,L=2mm,100mmHg) {q_cli*1e6:.2f} mL/s x 0.5 s = {q_cli*0.5e6:.2f} mL/beat")
print(f"  gap sensitivity (h^3): 40 um {(40/50)**3-1:+.1%}, 60 um {(60/50)**3-1:+.1%}")
verdict(True, "0.38 mL/beat (Inventor) and 0.16 mL/beat (Clinical) both reproduced; factor 2.4 comes only from the "
        "UNVERIFIED w/L (6/1 vs 4/2 mm). True range 0.15-0.4 mL/beat until w and L are measured in CAD. Clinical's "
        "+/-10 um -> -49%/+73% leakage also reproduced. Inventor's 0.5% of stroke volume holds for SV~75-80 mL.")
print(f"  other Inventor numbers: R dth/h at 2.2 mm = {R_ear*DEG/H:.3f}/deg, 15 deg -> {R_ear*15*DEG/H:.1f}, slide "
      f"{R_ear*15*DEG*1e6:.0f} um; h^2/nu = {H**2/NU*1e3:.3f} ms; entrance loss Ke=1.5 at 100 mmHg L=1 mm: "
      f"{1.5*RHO*(H**2*100*MMHG/1e-3/(12*MU))**2/2/(100*MMHG)*100:.1f}%")

# ---------------------------------------------------------------------------------------------
hdr("R2-C4", "Clinical: stress accumulation per pass, passes to threshold, power-law damage, Carreau")
L_s = 2e-3; dp = 100 * MMHG; G = dp / L_s
tau_w = H / 2 * G; Um = H ** 2 * G / (12 * MU); t_tr = L_s / Um
print(f"  tau_w = {tau_w:.1f} Pa, U_mean = {Um:.3f} m/s, transit L/U_mean = {t_tr*1e3:.2f} ms, "
      f"tau_w * t = {tau_w*t_tr:.3f} Pa s  (= 6 mu L / h = {6*MU*L_s/H:.3f}, INDEPENDENT of dp)")
# flux-weighted mean SA over streamlines: tau(y)=G|h/2-y|, time L/u(y), weight u(y)dy/Q -> L*int(tau)dy/Q = 3 mu L/h
SA_fw = 3 * MU * L_s / H
y = np.linspace(1e-9, H - 1e-9, 200001); u = G * y * (H - y) / (2 * MU); tau = G * np.abs(H / 2 - y)
SA_y = tau * L_s / u; SA_num = np.trapezoid(SA_y * u, y) / np.trapezoid(u, y)
print(f"  flux-weighted mean SA = 3 mu L/h = {SA_fw:.3f} Pa s (numerical {SA_num:.3f}); "
      "(centre streamline: tau=0 -> SA=0; near-wall streamlines -> large SA but carry little flux)")
thr = 3.5  # 35 dyn s/cm2
print(f"  passes to 3.5 Pa s: wall-shear x mean-time {thr/(tau_w*t_tr):.1f}; flux-weighted {thr/SA_fw:.1f}; "
      f"L=1 mm halves SA -> {thr/(SA_fw/2):.0f} passes")
verdict(True, "Clinical's 0.84 Pa s and ~4 passes reproduced as tau_w x mean transit time = 6 mu L/h. This is a "
        "conservative (upper) estimate: the flux-weighted mean over streamlines is 3 mu L/h = 0.42 Pa s (~8 passes). "
        "Key consequence: linear SA per pass in a Poiseuille slot does NOT depend on dp, only on L/h -- narrower gap "
        "=> more SA per pass (consistent with the Jun 2014 trend Clinical quotes; needs literature check). "
        "Multiple passes require recirculation, so residence time in the pocket is the real driver.")
# Giersiepen power laws (constants as commonly quoted; needs literature check): tau [Pa], t [s], result in %
HI = 3.62e-5 * tau_w ** 2.416 * t_tr ** 0.785
PLT = 3.31e-6 * tau_w ** 3.075 * t_tr ** 0.77
print(f"  Giersiepen HI = 3.62e-5 tau^2.416 t^0.785 = {HI:.3f}%; platelet LDH = 3.31e-6 tau^3.075 t^0.77 = {PLT:.3f}%")
verdict(abs(HI - 0.13) < 0.03 and abs(PLT - 0.38) < 0.08, "Clinical's HI ~0.13% and LDH ~0.38% reproduced (with the "
        "constants above, which need literature check). Ranking use only, as Clinical says.")
# Carreau (Cho & Kensey constants as commonly quoted; needs literature check)
mu0, muinf, lam, nC = 0.056, 0.00345, 3.313, 0.3568
def carreau(g):
    return muinf + (mu0 - muinf) * (1 + (lam * g) ** 2) ** ((nC - 1) / 2)
print(f"  {'gamma_dot':>10} {'mu_Carreau':>11} {'vs mu_inf':>10} {'vs 3.5 mPa s':>13}  (Clinical claim)")
for g, cl in [(4.8e4, 0.7), (7.6e3, 2.2), (1e3, 8), (250, 20), (100, 37), (50, 57)]:
    m = carreau(g)
    print(f"  {g:10.0f} {m*1e3:11.3f} {(m/muinf-1)*100:9.1f}% {(m/MU-1)*100:12.1f}%  ({cl}%)")
print(f"  gap wall shear rate at 100 mmHg, L=2 mm: {tau_w/MU:.2e} 1/s")
verdict(True, "Clinical's Carreau percentages are reproduced relative to mu_inf = 3.45 mPa s. Relative to the model's "
        "3.5 mPa s the gap error is -0.8% (Newtonian slightly OVER-predicts gap viscosity), and the low-shear pocket "
        "errors become +18%/+35%/+55% at 250/100/50 1/s. Conclusion unchanged: Newtonian is fine in the gap, Carreau "
        "sensitivity is needed for residence time.")

# ---------------------------------------------------------------------------------------------
hdr("R2-C5", "Inventor N4: mid-plane sliding split of the end gaps; 5 um concentricity; mesh count")
# Geometric argument: rigid rotation Q(theta) about the pivot axis (z) maps every plane z=const onto itself and every
# cylinder r=const (about the pivot axis) onto itself. A subdomain bounded only by such surfaces is invariant, so
# its mesh can rotate RIGIDLY (zero distortion) and a sliding identity pair on its boundary stays geometrically
# conforming. Check numerically: rotate points of the mid-plane and of a cylinder, measure departure.
rng = np.random.default_rng(0)
pts = np.c_[rng.uniform(-2.3e-3, 2.3e-3, (1000, 2)), np.full(1000, 10.825e-3)]
th = 15 * DEG; Rz = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1]])
p2 = pts @ Rz.T
print(f"  max |dz| of mid-plane points after 15 deg rotation: {np.max(np.abs(p2[:,2]-pts[:,2])):.1e} m (plane invariant)")
print(f"  max in-plane slide of those points (r up to 3.25 mm): {np.max(np.linalg.norm(p2[:,:2]-pts[:,:2],axis=1))*1e6:.0f} um "
      "(carried by the sliding pair, not by element distortion)")
verdict(True, "N4 removes end-gap SHEAR exactly (11.5 -> 0) because z=const planes are invariant under rotation about "
        "z. CONDITION: the leaflet-side rotating slab must be bounded laterally by a CYLINDER about the pivot axis "
        "(rotor-disc topology), not by the ear outline -- otherwise its side faces sweep through the fixed slab "
        "(Inventor's 'swept-footprint problem' is exactly this). So N4 = a thin rotor disc (radius >= ear extent "
        "~2.3 mm) with planar + cylindrical sliding pairs, i.e. rotating-machinery topology.")
# Concentricity: a cylindrical sliding interface at radius r_c inside a lateral gap must stay strictly between the two
# walls for the whole sweep: needs e (axis offset, moves wall by up to e*... actually by e over large angles: the
# wall arc centre moves on a circle of radius e, displacement 2 e sin(theta/2) ~ e theta) + 2 sigma_r (surface
# scatter both walls) < h/2 - delta_min (room for >=1 element layer on each side).
for e in [0, 5e-6, 10e-6, 25e-6]:
    for sig in [5e-6, 10e-6]:
        need = e * 15 * DEG + 2 * sig
        room = H / 2 - H / 6
        print(f"    e={e*1e6:4.0f} um sigma_r={sig*1e6:4.0f} um: wall excursion {need*1e6:5.1f} um vs room "
              f"{room*1e6:.1f} um -> {'OK' if need < room else 'FAILS'}")
verdict(True, "sigma_r ~5 um is an adequate (not unique) tolerance: with 15-deg sweep the cylinder survives if "
        "e*theta + 2 sigma_r < h/2 - h/6 = 16.7 um, i.e. sigma_r = 5 um allows e up to ~25 um. Strain in the lateral "
        "gaps e*theta/h: e=10/50/200 um -> %.2f/%.2f/%.2f (Inventor 0.05/0.26/1.05: PASS)." %
        (10e-6 * th / H, 50e-6 * th / H, 200e-6 * th / H))
# mesh cost of a 25 um slab, 3 layers
A_face = math.pi * 2.3e-3 ** 2; slab = A_face * 25e-6; a = 25e-6 / 3
n_tet = slab / (a ** 3 / (6 * math.sqrt(2)))
tri = 75e-6; n_prism = 3 * A_face / (math.sqrt(3) / 4 * tri ** 2)
print(f"  one 25 um slab, disc r=2.3 mm: isotropic tets (edge {a*1e6:.1f} um) ~{n_tet/1e6:.1f} M; "
      f"swept prisms (3 layers, 75 um triangles) ~{n_prism/1e3:.0f} k")
verdict(False, f"N4 is only affordable with SWEPT/prism slabs: isotropic tets need ~{n_tet/1e6:.0f} M elements per "
        "slab x 8 slabs (2 per end gap x 4 end gaps). Add to Inventor's risk list.")

# ---------------------------------------------------------------------------------------------
hdr("R2-C6", "Corrected validation statistics: TOST equivalence, CI of bias, slope CI, Bland-Altman, BH-FDR")
def validation_stats(bf, gh, margin=5.0, alpha=0.05, label=""):
    """bf, gh: arrays per angle (A/B already averaged). Works on % differences d = 100 (gh-bf)/bf.
    TOST: H0 |mu_d| >= margin. Two one-sided t tests, df=n-1. Equivalent iff the (1-2 alpha) CI lies inside +/-margin."""
    d = 100 * (gh - bf) / bf; n = len(d); m = d.mean(); sd = d.std(ddof=1) if n > 1 else float("nan")
    se = sd / math.sqrt(n); df = n - 1
    t_lo = (m + margin) / se; t_hi = (m - margin) / se
    p_lo = stats.t.sf(t_lo, df); p_hi = stats.t.cdf(t_hi, df); p_tost = max(p_lo, p_hi)
    ci90 = (m - stats.t.ppf(1 - alpha, df) * se, m + stats.t.ppf(1 - alpha, df) * se)
    ci95 = (m - stats.t.ppf(0.975, df) * se, m + stats.t.ppf(0.975, df) * se)
    lr = stats.linregress(bf, gh); ts = stats.t.ppf(0.975, n - 2)
    loa = (m - 1.96 * sd, m + 1.96 * sd); loa_u = stats.t.ppf(0.975, df) * math.sqrt(3 * sd ** 2 / n)
    print(f"  [{label}] d% = {np.round(d,2)}  mean {m:.2f}%  sd {sd:.2f}%")
    print(f"    TOST (margin +/-{margin}%): p = {p_tost:.4f}; 90% CI [{ci90[0]:.2f}, {ci90[1]:.2f}] -> "
          f"{'EQUIVALENT' if p_tost < alpha else 'NOT shown equivalent'}")
    print(f"    95% CI of bias [{ci95[0]:.2f}, {ci95[1]:.2f}]%;  slope {lr.slope:.4f} +/- {ts*lr.stderr:.4f} (95%), "
          f"intercept {lr.intercept:.5f} +/- {ts*lr.intercept_stderr:.5f}, R2 {lr.rvalue**2:.4f}")
    print(f"    Bland-Altman LoA [{loa[0]:.2f}, {loa[1]:.2f}]%, each +/-{loa_u:.2f}% (95%)")
    return p_tost
# (a) real data at 70 deg only: K=1e-9 ghost vs body-fitted
d70 = 100 * (np.array([0.1934, 0.1912]) - np.array([0.2117, 0.2095])) / np.array([0.2117, 0.2095])
print(f"  (a) 70 deg, K=1e-9 (halo present): d = {np.round(d70,2)}% (A, B). A and B are mirror images, not "
      "replicates: n_eff = 1 angle -> no t-based test is possible. Point estimate -8.69% lies OUTSIDE +/-5%: "
      "Phase A FAILS deterministically, no statistics needed.")
F70 = 100 * (61.81 - 49.205) / 49.205
print(f"      FyA at 70 deg: {F70:+.1f}% -> fails any force gate.")
# (b) synthetic 4-angle example (ASSUMPTION: post-FixChi bias ~ -3% with 0.8% scatter)
bfm = (TA + TB) / 2
d_syn = np.array([-3.1, -2.4, -3.8, -2.0])
gh_syn = bfm * (1 + d_syn / 100)
p_T = validation_stats(bfm, gh_syn, label="synthetic torque, bias~-3%")
d_syn2 = np.array([-4.6, -3.9, -5.2, -3.5])
p_T2 = validation_stats(bfm, bfm * (1 + d_syn2 / 100), label="synthetic torque, bias~-4.3%")
Fy_syn = FyA * (1 + np.array([6.0, 4.5, 8.0, 3.0]) / 100)
p_F = validation_stats(FyA, Fy_syn, margin=10.0, label="synthetic FyA, margin 10%")
# BH-FDR across outputs (TOST p-values)
outs = {"Tz (bias~-3%)": p_T, "Tz (bias~-4.3%)": p_T2, "FyA (10% margin)": p_F}
pv = np.array(list(outs.values())); order = np.argsort(pv); m = len(pv)
crit = 0.05 * np.arange(1, m + 1) / m; ok = pv[order] <= crit
k = (np.max(np.where(ok)[0]) + 1) if ok.any() else 0
print("  BH-FDR (q=0.05) over TOST p-values:")
for i, j in enumerate(order):
    print(f"    {list(outs)[j]:18s} p={pv[j]:.4f} crit={crit[i]:.4f} -> {'equivalent' if i < k else 'not shown'}")
verdict(True, "TOST with n=4 can demonstrate |bias|<5% only if mean bias is <~3-3.5% with sd<~1%: the -4.3% "
        "example is NOT shown equivalent even though every point but one is inside 5%. With numerical uncertainty, "
        "replace sd-based CI by the V&V-20-style interval E +/- 2 u_val, u_val = sqrt(GCI_ghost^2 + GCI_bf^2)/Fs "
        "(needs literature check).")

# ---------------------------------------------------------------------------------------------
hdr("R2-C7", "Other Round-1 numeric conflicts: halo thickness, added inertia")
for t, lab in [(0.7, "Numerical nominal"), (0.9, "Advocate"), (1.1, "half-disk upper bound")]:
    print(f"  halo = dV t / (2 V) with t={t} mm ({lab}): {90.8*t/(2*208.5)*1e3:.0f} um")
verdict(True, "Advocate ~200 um (t=0.9 mm) and Numerical 150-240 um are the SAME formula; only t differs. "
        "Halo = 3-5 gap widths for any t in 0.7-1.1 mm: conclusion (hinge sealed) is robust.")
for c, s_, lab in [(10e-3, 21.6e-3, "Feasibility c=10, full span"), (10e-3, 0.8 * 21.6e-3, "c=10, 0.8 span"),
                   (8e-3, 0.6 * 21.6e-3, "c=8, 0.6 span")]:
    Ia = 9 * math.pi / 128 * RHO * c ** 4 * s_
    print(f"  I_a ({lab}) = {Ia:.2e} -> {Ia/9.5e-9:.1f}x")
verdict(True, "1.3-6x (Numerical), 1.5-6.3x (Challenger), 5.3x (Feasibility) are one formula with different "
        "chord/span. All > 1 -> same conclusion. Resolve by measurement: in the body-fitted model, prescribe a "
        "small angular acceleration with flow at rest and read I_a = -T/alpha (one short run).")

# ---------------------------------------------------------------------------------------------
hdr("R2-C8", "Proposed numeric gate table for the Judge")
gates = [
    ("G0 solver", "LinRes (direct) at every step", "< 1e-10"),
    ("G1 ghost geometry", "|V_A - 208.5|/208.5 and same for V_B", "< 1%"),
    ("G2 symmetry", "|TzA-TzB|/mean at each angle (4 sig figs)", "< 3%"),
    ("G3 bias (torque)", "TOST, n=4 angles (A/B averaged), margin", "+/-5%, p<0.05 (90% CI inside)"),
    ("G3b bias (force)", "TOST FyA margin", "+/-10%"),
    ("G4 proportional bias", "95% CI of slope ghost-vs-BF", "within [0.95, 1.05]"),
    ("G5 mesh (bulk)", "GCI_fine torque, current mesh = finest", "< 5%, asympt. ratio 0.9-1.1"),
    ("G5b mesh (gap)", "GCI_fine gap flux / peak WSS (body-fitted submodel)", "< 5% / < 10%"),
    ("G6 K sensitivity", "|T(1e-10) - T(3e-10)| / T", "< 2% (else report as K-limited)"),
    ("G7 ghost not used in gap", "Darcy ratio 12 K T/h^3 and slip ratio", "must be < 5% or gap taken from film/body-fitted"),
    ("G8 film validity", "eps*Re, h^2 omega/nu at evaluation points", "< 0.3 and < 0.3"),
    ("G9 film vs 3D", "gap flux / peak WSS (uniform-gap region)", "within 10% / 15%"),
    ("G10 frozen-angle window", "3 e dtheta / h (e from O6 fit)", "< 5% -> sets dtheta"),
    ("G11 one-way coupling", "hinge friction torque / blade torque", "< 1e-2"),
    ("G12 mass balance (D)", "|sum of in/out fluxes| / inflow", "< 1%"),
    ("G13 stabilization (E)", "change of peak torque / opening time with artificial inertia", "< 2%"),
    ("G14 N4 geometry", "e*theta_sweep + 2 sigma_r", "< h/3 = 16.7 um"),
    ("G15 rheology", "Carreau vs Newtonian change in residence time / WSS", "report; > 5% => Carreau in results"),
]
for g in gates:
    print(f"  {g[0]:26s} | {g[1]:58s} | {g[2]}")
print("\nEnd of Round-2 checks.")
