# R1 — Advocate opening: O1 (ghost leaflet + Reynolds thin-film hinge)

Scratch calculations: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/advocate/calc.py` (blood rho = 1060 kg/m3, mu = 3.5 mPa s, h = 50 um). Every number below comes from that script or from the handoff, unless I mark it UNVERIFIED.

## 1. The thesis: split the problem by scale, solving each scale where its physics is valid

The real hinge spans roughly three orders of magnitude, from a 50 um gap to a 23 mm valve. Every single-mesh method fails in that 50 um region. ALE produced inverted elements at the ear end faces from about 3 deg onward. The shear there is gamma ~ R*dtheta/h = 0.66 per degree, which I checked. Overset meshing still needs overlap inside the gap. A ghost leaflet alone would need sqrt(K) << h. O1 does not try to fix the gap with more mesh. It removes the gap from the 3D problem and solves it with the asymptotically correct equation for a thin gap. That equation is the Reynolds lubrication equation, which is exact to leading order in eps = h/L. It runs on a 2D surface mesh with essentially no 3D cost, and its error has a known form.

## 2. Why the K = 1e-9 result is promising, and why it is not yet evidence

Encouraging facts (handoff §5): torque signs are correct, A/B symmetry is 1%, the time trend is smooth, and the run took 1 h 20 min instead of 10+ h.

The issue we have to face first is the halo. The ghost volume is 299.3 mm3 against 208.5 mm3, an excess of 90.8 mm3 per leaflet. Spread over an assumed plate area of about 232 mm2 per face (thickness 0.9 mm, UNVERIFIED), the halo is **about 200 um thick, four times the hinge gap**. The current ghost run therefore almost certainly has the **hinge sealed**. Its -8.6% torque and +25.6% FyA belong to a fattened blade with a closed hinge, so they say nothing about how accurate the method is. The useful consequence is that the error budget is dominated by a known geometric bug that can be fixed (FixChi, not yet run), not by penalization physics.

Penalization physics by itself is benign for the blade:
- The Brinkman penetration length is sqrt(K) = 31.6 / 17.3 / 10.0 um for K = 1e-9 / 3e-10 / 1e-10.
- The blade boundary layer at peak systole is delta ~ c/sqrt(Re_c) ~ 166 um (U = 1.2 m/s and c = 10 mm assumed, UNVERIFIED). That gives sqrt(K)/delta = 0.19 at K = 1e-9 and 0.06 at K = 1e-10.
- The penalization error is known to decay like a power of the penalty parameter: Angot, Bruneau & Fabrie, *Numer. Math.* 81:497–520 (1999), DOI 10.1007/s002110050401. The exact exponent depends on norm and region (1/2, 1/4 and 1 are all reported), so we need the K study to measure it rather than assume it.
- The stiffness cost is tau = rho*K/mu = 3.0e-4, 9.1e-5 and 3.0e-5 s for those three K values. K = 1e-10 is about 10x stiffer than the run that worked, which is still affordable. K = 1e-11 (3 us) is where the earlier run stalled.

My prediction, which is falsifiable: with sharp chi (V = 208.5 mm3), the fixed-angle torque bias at 70 deg drops below 5%, and a three-K Richardson extrapolation K -> 0 lands within the GCI band of the body-fitted reference. If that fails at Phase A, O1 loses and we fall back to O4. That is a clean gate and does not depend on O1 succeeding.

## 3. Thin-film physics in the 50 um gap

**Equation.** With s the in-surface coordinates on the pocket wall, h(s,t) the local clearance, and U_1, U_2 the ear and pocket surface velocities (pocket U_2 = 0):

  d/dt(h) + div_s( -h^3/(12 mu) grad_s p + h (U_1+U_2)/2 ) = 0

where h(s,t) = d0(R(theta)^-1 s) is taken from the CAD clearance map, so the squeeze term dh/dt carries the eccentric-motion effect.

**Wall shear** on the ear (upper sign) and the pocket (lower sign):

  tau_w = -/+ (h/2) dp/ds + mu (U_1 - U_2)/h

**Validity**, from my own numbers, with dp the hinge pressure drop and L the path length:

| Case | Re_h = rho*U_mean*h/mu | eps = h/L | eps*Re | tau_pressure |
|---|---|---|---|---|
| systole-like dp = 10 mmHg, L = 2 mm | 0.6 | 0.025 | 0.015 | 17 Pa |
| diastolic leakage dp = 80 mmHg, L = 2 mm | 4.8 | 0.025 | 0.12 | 133 Pa |
| dp = 100 mmHg, L = 1 mm (worst short path) | 12 | 0.05 | 0.60 | 333 Pa |

The unsteady parameter is h^2*omega/nu = 0.16 even for a 30 ms closure transient, so quasi-steady lubrication holds. Couette shear from leaflet rotation is small: 0.6–1.2 Pa at 15 deg per 60/30 ms with an ear radius of 2 mm. **Hinge WSS is therefore pressure-driven**, and the film is valid (eps*Re << 1) for most of the cycle. It is marginal (about 0.5) only for the shortest paths at peak diastolic pressure. Shear rates near 200 Pa are about 6e4 1/s, so the Newtonian assumption is reasonable there.

A supporting scaling result: a hinge-flow study (Simon et al., Ann. Biomed. Eng. 2004, found via search, details needs literature check) reportedly estimated Re about 150 in a clinical-size gap. Pressure-driven Re scales as h^3, so 150 at a 150 um gap becomes **5.6 at 50 um**, consistent with my table. If their gap was 100 um, it becomes 19, which is where the film starts to break down. **We must pin down the gap they used.**

## 4. Coupling to the 3D ghost model

1. **Partition.** Divide the pocket into three parts: the uniform-gap film region Gamma_f (2D shell), and "mouth" and "corner" regions, about 3–5 h deep, that stay in the 3D mesh. The mouths are where eps-validity breaks down, and keeping them in 3D is the honest choice.
2. **Pressure BCs.** At each mouth curve dGamma_f, set the film pressure p_f equal to the 3D pressure averaged over a short arc. This is a one-way Dirichlet condition in Phase C and a fixed-point exchange in Phase D.
3. **Flux matching.** The film edge flux q_n = -h^3/(12mu) dp/dn + h*U_n/2 is returned to the 3D model as a boundary mass source or sink on the mouth faces. The check is that the net film flux equals the 3D mouth flux to within 1% (Phase D gate).
4. **No double counting.** Inside Gamma_f, the 3D pocket volume is penalized as solid (chi = 1, u_s = 0) or removed from the 3D geometry. Gap flow then exists only in the film. The torque from the hinge is computed from film tractions only, M = integral of r x (-p n + tau_w t) dA, and is added to the ghost-blade torque.
5. **Size of the hinge torque.** Shear torque is about 200 Pa x 2 mm2 x 1 mm = 4e-4 mN m, roughly **0.2% of the blade torque** (order-of-magnitude, UNVERIFIED geometry). Pressure torque on a near-concentric ear is also small. **One-way coupling of the film is therefore justified for leaflet kinematics.** The film's contribution is WSS, leakage and washout, not dynamics. This greatly reduces coupling risk.
6. **Tooling.** COMSOL has a Thin-Film Flow (shell) interface that solves the Reynolds equation on boundaries. Whether the COE license includes it is UNVERIFIED and needs to be checked today. If it does not, a Coefficient Form Boundary PDE implements the equation above directly.

A leakage sanity number: per hinge, q is about 5.7 mL/min at 80 mmHg (width 3 mm assumed, UNVERIFIED). Clinical leakage volumes per beat are a reference for comparison (needs literature check).

## 5. What is novel

- To my knowledge (needs literature check; my searches found micro-PIV, immersed-boundary and CFD hinge studies, but no lubrication or Reynolds-film model of an MHV hinge coupled to blade FSI), the first **asymptotically justified reduced-order hinge model** coupled to a moving-leaflet simulation, with explicit validity maps of eps*Re(s,t).
- A **closed-form hinge WSS and leakage map** across the cycle at a physiological 50 um clearance, usable for blood-damage indices, which 3D methods cannot resolve on our hardware.
- A **quantified ghost-leaflet error model** (K exponent, halo correction, GCI) against a body-fitted reference. This is publishable even if Phases D–E slip.

## 6. Honest weaknesses of O1

1. **Unproven as a combination.** Each part works in isolation, but the coupled hinge does not exist yet. The thesis must not depend on Phases D–E.
2. **The current ghost numbers are not validation.** The halo probably seals the hinge, and the FyA error (+25.6%) has the opposite sign to the torque error (-8.6%), which is not yet explained.
3. **Mouths, corners and the ear end faces** (z = +/-10.8 mm, the same place ALE fails) fall outside lubrication validity. The 3D treatment there inherits the ghost resolution limit: sqrt(K) is 0.2–0.63 h.
4. **Closure transients.** Water-hammer pressure spikes and possible cavitation at closure can exceed the 100 mmHg I assumed, which pushes eps*Re past 0.5 (UNVERIFIED magnitude).
5. **Continuum limit.** h/d_RBC is about 6, so cell-scale effects (plasma skimming) are not captured by any continuum method, O1 included.
6. **Moving sharp chi** gives no Jacobian with respect to theta and produces force noise. Two-way FSI (Phase E) is the riskiest step.
7. **The reference is not mesh-converged.** All "-x%" claims are provisional until the body-fitted GCI is done.

## 7. The other options from my perspective

- **O2 (ghost alone):** necessary (it is O1's Phase A/B) but not sufficient. Resolving a 50 um gap needs sqrt(K) of about 3 um (K = 1e-11), which stalled at 8.5 us steps. It cannot give gap WSS honestly.
- **O3 (ALE + remesh):** valuable as a kinematics and remeshing result. The inverted elements at the ear ends are structural, because sliding end faces cross a 50 um gap, and 5–6 h restart intervals with flow OFF make flow-on FSI to 15 deg unrealistic in 4 weeks.
- **O4 (two-scale):** the strongest competitor and the right fallback. My 0.2% hinge-torque estimate supports decoupling. O1 is essentially O4 with an asymptotic hinge model in place of a separate 3D hinge solve, and the two share Phases A–B.
- **O5 (manual overset):** not native to COMSOL, about 2x DOF, still needs overlap inside the gap. Too risky in 4 weeks.
- **O6 (pivot optimization):** cheap and complementary. It may explain the end-face inversions and gives the correct h(s,t) for the film. It should run in parallel, not be chosen instead.

## Bottom line (ranking)

1. **O1**: the only option that gives accurate hinge WSS/leakage at 50 um on desktop hardware with a quantifiable error, staged behind gates so it cannot sink the thesis.
2. **O4**: the safest publishable result and the built-in fallback for O1. It shares O1's validation phases.
3. **O2**: the mandatory foundation (K study, GCI), but cannot resolve the gap on its own.
4. **O6**: cheap, runs in parallel, and improves the film's h(s,t) and the ALE diagnosis.
5. **O3**: report as a kinematics/remeshing contribution; not viable for flow near the ear ends.
6. **O5**: deferred to future work (no native support, cost, overlap still inside the gap).
