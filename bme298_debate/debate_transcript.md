# Debate transcript — BME 298 real 50 um hinge-gap FSI: method selection

Orchestrated Nov 10, 2026. Input: HANDOFF_BME298_2026-11-10.md. Common rules: BRIEF_common.md.

## Team
| Role | Brief |
|---|---|
| Advocate | Argues FOR hybrid ghost leaflet + Reynolds thin-film hinge (O1) |
| Challenger | Skeptical top-journal reviewer; attacks O1, proposes strongest alternative |
| Numerical analyst | Verifies every number in Python (math_checks.py / math_checks_output.txt) |
| Clinical/physiology expert | Which outputs matter for hinge thrombosis; can each method deliver them |
| Feasibility engineer | COMSOL 6.4 implementability, hardware, 4-week schedule |
| Inventor | Proposes new approaches with math |
| Judge (separate, did not argue) | Writes verdict.md |

Literature note: journal full-text sites (PMC, Springer) and doc.comsol.com were blocked in this session; citations marked as checked were verified only via search results/abstracts unless stated otherwise. Treat all as 'needs literature check' before citing in the thesis.

---

# ROUND 1 — Opening arguments

---

## Round 1 — Advocate

### R1 — Advocate opening: O1 (ghost leaflet + Reynolds thin-film hinge)

Scratch calculations: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/advocate/calc.py` (blood rho = 1060 kg/m3, mu = 3.5 mPa s, h = 50 um). Every number below comes from that script or from the handoff, unless I mark it UNVERIFIED.

#### 1. The thesis: split the problem by scale, solving each scale where its physics is valid

The real hinge spans roughly three orders of magnitude, from a 50 um gap to a 23 mm valve. Every single-mesh method fails in that 50 um region. ALE produced inverted elements at the ear end faces from about 3 deg onward. The shear there is gamma ~ R*dtheta/h = 0.66 per degree, which I checked. Overset meshing still needs overlap inside the gap. A ghost leaflet alone would need sqrt(K) << h. O1 does not try to fix the gap with more mesh. It removes the gap from the 3D problem and solves it with the asymptotically correct equation for a thin gap. That equation is the Reynolds lubrication equation, which is exact to leading order in eps = h/L. It runs on a 2D surface mesh with essentially no 3D cost, and its error has a known form.

#### 2. Why the K = 1e-9 result is promising, and why it is not yet evidence

Encouraging facts (handoff §5): torque signs are correct, A/B symmetry is 1%, the time trend is smooth, and the run took 1 h 20 min instead of 10+ h.

The issue we have to face first is the halo. The ghost volume is 299.3 mm3 against 208.5 mm3, an excess of 90.8 mm3 per leaflet. Spread over an assumed plate area of about 232 mm2 per face (thickness 0.9 mm, UNVERIFIED), the halo is **about 200 um thick, four times the hinge gap**. The current ghost run therefore almost certainly has the **hinge sealed**. Its -8.6% torque and +25.6% FyA belong to a fattened blade with a closed hinge, so they say nothing about how accurate the method is. The useful consequence is that the error budget is dominated by a known geometric bug that can be fixed (FixChi, not yet run), not by penalization physics.

Penalization physics by itself is benign for the blade:
- The Brinkman penetration length is sqrt(K) = 31.6 / 17.3 / 10.0 um for K = 1e-9 / 3e-10 / 1e-10.
- The blade boundary layer at peak systole is delta ~ c/sqrt(Re_c) ~ 166 um (U = 1.2 m/s and c = 10 mm assumed, UNVERIFIED). That gives sqrt(K)/delta = 0.19 at K = 1e-9 and 0.06 at K = 1e-10.
- The penalization error is known to decay like a power of the penalty parameter: Angot, Bruneau & Fabrie, *Numer. Math.* 81:497–520 (1999), DOI 10.1007/s002110050401. The exact exponent depends on norm and region (1/2, 1/4 and 1 are all reported), so we need the K study to measure it rather than assume it.
- The stiffness cost is tau = rho*K/mu = 3.0e-4, 9.1e-5 and 3.0e-5 s for those three K values. K = 1e-10 is about 10x stiffer than the run that worked, which is still affordable. K = 1e-11 (3 us) is where the earlier run stalled.

My prediction, which is falsifiable: with sharp chi (V = 208.5 mm3), the fixed-angle torque bias at 70 deg drops below 5%, and a three-K Richardson extrapolation K -> 0 lands within the GCI band of the body-fitted reference. If that fails at Phase A, O1 loses and we fall back to O4. That is a clean gate and does not depend on O1 succeeding.

#### 3. Thin-film physics in the 50 um gap

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

#### 4. Coupling to the 3D ghost model

1. **Partition.** Divide the pocket into three parts: the uniform-gap film region Gamma_f (2D shell), and "mouth" and "corner" regions, about 3–5 h deep, that stay in the 3D mesh. The mouths are where eps-validity breaks down, and keeping them in 3D is the honest choice.
2. **Pressure BCs.** At each mouth curve dGamma_f, set the film pressure p_f equal to the 3D pressure averaged over a short arc. This is a one-way Dirichlet condition in Phase C and a fixed-point exchange in Phase D.
3. **Flux matching.** The film edge flux q_n = -h^3/(12mu) dp/dn + h*U_n/2 is returned to the 3D model as a boundary mass source or sink on the mouth faces. The check is that the net film flux equals the 3D mouth flux to within 1% (Phase D gate).
4. **No double counting.** Inside Gamma_f, the 3D pocket volume is penalized as solid (chi = 1, u_s = 0) or removed from the 3D geometry. Gap flow then exists only in the film. The torque from the hinge is computed from film tractions only, M = integral of r x (-p n + tau_w t) dA, and is added to the ghost-blade torque.
5. **Size of the hinge torque.** Shear torque is about 200 Pa x 2 mm2 x 1 mm = 4e-4 mN m, roughly **0.2% of the blade torque** (order-of-magnitude, UNVERIFIED geometry). Pressure torque on a near-concentric ear is also small. **One-way coupling of the film is therefore justified for leaflet kinematics.** The film's contribution is WSS, leakage and washout, not dynamics. This greatly reduces coupling risk.
6. **Tooling.** COMSOL has a Thin-Film Flow (shell) interface that solves the Reynolds equation on boundaries. Whether the COE license includes it is UNVERIFIED and needs to be checked today. If it does not, a Coefficient Form Boundary PDE implements the equation above directly.

A leakage sanity number: per hinge, q is about 5.7 mL/min at 80 mmHg (width 3 mm assumed, UNVERIFIED). Clinical leakage volumes per beat are a reference for comparison (needs literature check).

#### 5. What is novel

- To my knowledge (needs literature check; my searches found micro-PIV, immersed-boundary and CFD hinge studies, but no lubrication or Reynolds-film model of an MHV hinge coupled to blade FSI), the first **asymptotically justified reduced-order hinge model** coupled to a moving-leaflet simulation, with explicit validity maps of eps*Re(s,t).
- A **closed-form hinge WSS and leakage map** across the cycle at a physiological 50 um clearance, usable for blood-damage indices, which 3D methods cannot resolve on our hardware.
- A **quantified ghost-leaflet error model** (K exponent, halo correction, GCI) against a body-fitted reference. This is publishable even if Phases D–E slip.

#### 6. Honest weaknesses of O1

1. **Unproven as a combination.** Each part works in isolation, but the coupled hinge does not exist yet. The thesis must not depend on Phases D–E.
2. **The current ghost numbers are not validation.** The halo probably seals the hinge, and the FyA error (+25.6%) has the opposite sign to the torque error (-8.6%), which is not yet explained.
3. **Mouths, corners and the ear end faces** (z = +/-10.8 mm, the same place ALE fails) fall outside lubrication validity. The 3D treatment there inherits the ghost resolution limit: sqrt(K) is 0.2–0.63 h.
4. **Closure transients.** Water-hammer pressure spikes and possible cavitation at closure can exceed the 100 mmHg I assumed, which pushes eps*Re past 0.5 (UNVERIFIED magnitude).
5. **Continuum limit.** h/d_RBC is about 6, so cell-scale effects (plasma skimming) are not captured by any continuum method, O1 included.
6. **Moving sharp chi** gives no Jacobian with respect to theta and produces force noise. Two-way FSI (Phase E) is the riskiest step.
7. **The reference is not mesh-converged.** All "-x%" claims are provisional until the body-fitted GCI is done.

#### 7. The other options from my perspective

- **O2 (ghost alone):** necessary (it is O1's Phase A/B) but not sufficient. Resolving a 50 um gap needs sqrt(K) of about 3 um (K = 1e-11), which stalled at 8.5 us steps. It cannot give gap WSS honestly.
- **O3 (ALE + remesh):** valuable as a kinematics and remeshing result. The inverted elements at the ear ends are structural, because sliding end faces cross a 50 um gap, and 5–6 h restart intervals with flow OFF make flow-on FSI to 15 deg unrealistic in 4 weeks.
- **O4 (two-scale):** the strongest competitor and the right fallback. My 0.2% hinge-torque estimate supports decoupling. O1 is essentially O4 with an asymptotic hinge model in place of a separate 3D hinge solve, and the two share Phases A–B.
- **O5 (manual overset):** not native to COMSOL, about 2x DOF, still needs overlap inside the gap. Too risky in 4 weeks.
- **O6 (pivot optimization):** cheap and complementary. It may explain the end-face inversions and gives the correct h(s,t) for the film. It should run in parallel, not be chosen instead.

#### Bottom line (ranking)

1. **O1**: the only option that gives accurate hinge WSS/leakage at 50 um on desktop hardware with a quantifiable error, staged behind gates so it cannot sink the thesis.
2. **O4**: the safest publishable result and the built-in fallback for O1. It shares O1's validation phases.
3. **O2**: the mandatory foundation (K study, GCI), but cannot resolve the gap on its own.
4. **O6**: cheap, runs in parallel, and improves the film's h(s,t) and the ALE diagnosis.
5. **O3**: report as a kinematics/remeshing contribution; not viable for flow near the ear ends.
6. **O5**: deferred to future work (no native support, cost, overlap still inside the gap).

---

## Round 1 — Challenger

### R1 — Challenger opening: the hybrid ghost-leaflet + Reynolds hinge (O1) is not reviewer-safe

Role: skeptical reviewer (J Biomech / Ann Biomed Eng / CMBBE level). All numbers are from
`/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/challenger/calc.py` and `calc2.py`
unless marked UNVERIFIED. Blood: mu = 3.5e-3 Pa s, rho = 1060 kg/m3. Gap h = 50 um.

#### 1. What is genuinely good about O1 (concessions)
- A fixed blood mesh removes the failure seen in O3 (inverted elements at the ear end faces z = +/-10.8 mm).
- Each run is cheap (K = 1e-9: 1 h 20 min). Signs are correct, A/B symmetry is 1%, and the time trend is smooth.
- The Reynolds film is physically sound where the gap is uniform: eps = h/L ~ 0.025, and the gap Womersley number
  is alpha = h*sqrt(omega/nu) = 0.074 (computed), so the hinge flow is quasi-steady.
- The phased plan with go/no-go gates and a two-scale fallback is honest and well designed.

The problem is that the ghost-leaflet half of O1 cannot represent a 50 um gap at any K this desktop can run.
O1 therefore cannot deliver the "real 50 um hinge" claim that the thesis title makes.

#### 2. Specific attacks

##### 2.1 Penalty-layer resolution: three length scales in the wrong order
| Scale | Value |
|---|---|
| sqrt(K) at K = 1e-9 | 31.6 um (0.63 h) |
| Mean regular-tet edge (54,250 mm3 / 2.05 M elements) | ~608 um (domain mean; the local size at the hinge is not stated in the handoff) |
| Halo thickness implied by the +91 mm3 excess volume (leaflet t = 0.6–1.0 mm, UNVERIFIED) | 131–219 um |

Brinkman error estimates assume the mesh resolves the sqrt(K) layer. The ordering here is
h_mesh (near the leaflet ~150 um) > h_gap (50 um) > sqrt(K) (32 um). The penalty layer is unresolved, and the
gap is two to four times narrower than the cell that has to hold the chi transition. The volume-penalization
literature reports a model error of order eta^(1/2), with the dominant error coming from a displacement length
(Hester, Vasil & Burns, arXiv:1903.11914). In a 50 um gap, a displacement of ~sqrt(K) is an O(1) change in the
gap height.

##### 2.2 Mass leakage: the penalized ear is a better channel than the real gap
I solved the 1-D Brinkman composite problem: a 50 um fluid gap next to a 300 um penalized layer, driven by the
same dp/dx. Results relative to the exact no-slip Poiseuille flux:

| K | Flux in the gap | Extra flux through the "solid" | Total |
|---|---|---|---|
| 1e-9 | 3.63x | 26.1x | 29.7x |
| 1e-10 | 1.70x | 2.97x | 4.7x |
| 1e-11 | 1.20x | 0.31x | 1.5x |

A Darcy estimate of the flow parallel to the gap inside an ear of thickness T, K*T/(h^3/12), agrees with this
(57.6x at K = 1e-9 with T = 0.6 mm). Keeping the porous-ear leak below 5% of the true gap flux needs
**K <= 8.7e-13 (sqrt(K) = 0.9 um, tau = rho*K/mu = 0.26 us)**. That is ten times stiffer than the K = 1e-11 run,
which already froze at dt = 8.5e-6 s.

The transleaflet leak at closure is worse. With dp = 80 mmHg across one ~298 mm2 face (UNVERIFIED geometry), the
Darcy leak is 1,298 mL/s at K = 1e-9 and 13 mL/s at K = 1e-11, against a cardiac output of ~83 mL/s. **The ghost
leaflet cannot hold closure pressure.** O1's closure-phase hinge results (eps*Re ~ 0.24, tau_w ~ 270 Pa) therefore
rest on a bulk model that leaks more than the hinge carries. Systole is milder: Darcy velocity is 18–240 mm/s for
dp = 50–500 Pa (UNVERIFIED dp), which is a few percent to ~25% of the bulk velocity.

##### 2.3 Force and torque from the penalization integral are not the surface traction
The penalization integral equals the surface traction plus the momentum of the fluid inside the porous body.
Halo size and K change the two parts differently. That fits the observed pattern of **opposite-signed errors:
|Tz| -8.6% and FyA +25.6%**. A uniform halo would bias both in the same direction. One case is reported where
drag from the penalization term came out at about half the reference value ("Asymptotic analysis of hydrodynamic
forces in a Brinkman penalization method...", J Fluid Mech 2024; authors and pages need literature check).
A reviewer will ask for a traction-based check: an integral over a surface slightly outside chi, using the
momentum-flux form. A 25.6% force error at the CAD angle, on the model meant to be validated, is a reject-level
number unless it is shown to disappear when K -> 0 and h -> 0.

##### 2.4 Halo, and sharp versus smoothed chi
The FixChi fix (sharp labels) has not been run. Even if it works, the plan uses sharp chi for validation and a
smoothed, volume-calibrated chi for motion. **That validates a different operator from the one used in the
results.** A reviewer will treat the Phase D numbers as unvalidated.

##### 2.5 Lubrication at mouths, corners and closure
Reynolds needs pressure boundary conditions at the film mouths. In O1 those pressures come from the ghost model
at exactly the place where it is least accurate (a ~150 um cell holding the chi transition). The Couette term is
small: about 1.3 Pa at 15 deg / 30 ms closure with R = 2.2 mm (UNVERIFIED kinematics). So the hinge WSS is
dominated by the (h/2)dp/ds term and inherits the mouth-pressure error directly. Corners, the ear end faces (the
same region where O3 inverted) and any contact or near-contact at closure break h << L. The film is valid where
it matters least.

##### 2.6 Coupling and double counting
Phase D "blocks" the 3D pocket. If it is blocked with penalization, Section 2.2 applies: the blocked pocket
leaks 3–77x the film flux at K = 1e-9 to 1e-10. If the film flux is not subtracted from the 3D flow, the gap flow
is counted twice. There is no conservative interface formulation. "Mass balance < 1%" is a gate with no design
behind it.

##### 2.7 No mesh-converged reference
The body-fitted fixed-angle "truth" has no GCI. Every Phase A gate (bias < 5%) is measured against an unconverged
target. This is the first comment any reviewer will make.

##### 2.8 Two-way FSI without the d(chi)/d(theta) Jacobian
General Extrusion is not differentiated with respect to theta, so the torque-angle coupling is missing from Newton.
A potential-flow flat-plate estimate puts the added rotational inertia at **1.5–6.3x the leaflet inertia
(9.5e-9 kg m2)**. Chord, span and axis offset are UNVERIFIED and this is an order-of-magnitude estimate only. At
an added-mass ratio above 1, a scheme that is effectively partitioned is known to be unstable or to need heavy
artificial inertia or damping. The Phase E gate "stabilization effect < 2%" will then fail or cannot be checked.

##### 2.9 Timeline (24 days to Dec 4)
- The K study needs K <= 1e-12 to control the leak (Section 2.2). Even K = 1e-10 is ~13 h per run if dt scales
  with tau (UNVERIFIED), and K = 1e-11 already failed.
- A 3-mesh GCI at r = 1.3 means 3.5 M and 7.6 M DOF. Direct-solver memory grows about 2.9x and 8.2x
  (N^(4/3) scaling, UNVERIFIED constant), which is beyond 32 GB. Iterative solvers have to be tuned for the
  penalized operator.
- Phases C, D and E each contain an unsolved methods problem.

Realistic: A and possibly B. C–E will not pass their gates by December.

#### 3. Proposed alternative (O7): body-fitted two-scale with quasi-steady moving-wall hinge submodel (O4 + O6)
1. **Global (blade scale):** leaflet dynamics from a reduced-order rigid-body ODE,
   I*theta'' = T_bf(theta, t) - I_a*theta'' - c*theta'. T_bf comes from the existing body-fitted fixed-angle
   torque table (70/75/80/84 deg, interpolated). I_a is the added inertia, estimated or computed. This runs in
   minutes. Cross-check against the v21 pin-joint model and the ALE kinematics.
2. **O6 pivot fit:** fit rotation plus sliding to the CAD pocket so the ears do not intersect the housing. This
   explains the O3 inversions and fixes the wall kinematics.
3. **Local (real 50 um hinge):** a small body-fitted hinge-pocket domain at a set of angles. Apply moving-wall
   (u = omega x r + sliding) velocity boundary conditions, and take mouth pressures from the global model.
   Quasi-steadiness is justified by the gap Womersley number of 0.074. A small domain makes a 3-mesh **GCI on the
   real gap** affordable on 4 cores.
4. **Novelty kept:** the Reynolds film becomes a **validated reduced model against 3D body-fitted truth** (O1's
   Phase C, done properly). Hinge WSS, flux and residence time are reported at 50 um, compared with the 150 um
   used by Simon et al.

Why a reviewer accepts this more readily:
- Every surface is body-fitted with no-slip, so there is no porous leak.
- It is mesh-converged where the claim is made.
- It follows a published precedent. Simon et al. drove pulsatile hinge CFD with leaflet motion and boundary
  conditions extracted from bulk FSI (Simon, Ge, Sotiropoulos, Yoganathan, "Numerical investigation of the
  performance of three hinge designs of bileaflet mechanical heart valves", Ann Biomed Eng 38(11):3295–3310, 2010,
  doi:10.1007/s10439-010-0086-3; also Simon et al., "Simulation of the three-dimensional hinge flow fields of a
  bileaflet mechanical heart valve under aortic conditions", Ann Biomed Eng, doi:10.1007/s10439-009-9857-0, gap
  ~150 um; volume and year need check).

Limitations to state openly:
- No two-way coupling at the hinge.
- dh/dt is captured only through the normal wall velocity on a frozen geometry, which is first order.

#### 4. Brief evaluation of the other options
- **O2 ghost alone:** fine for blade-scale torque (-8.6%) at mid-systole, provided GCI and K -> 0 are done.
  Disqualified for the gap and for closure (Section 2.2).
- **O3 ALE + remesh:** proves the kinematics but the mesh is invalid at the ear ends. Usable for blade-scale
  motion, not for hinge flow. Pairing it with the O6 pivot fit may remove the inversions (worth one run).
- **O4 two-scale:** the safest publishable core. This is the basis of my O7.
- **O5 manual overset:** no COMSOL built-in, overlap is needed inside 50 um, and DOF roughly doubles. Not
  feasible in 4 weeks.
- **O6 pivot/kinematics:** cheap, high value, and a prerequisite for any correct hinge wall velocity.

#### Bottom line (ranking)
1. **O7 (O4 + O6 + quasi-steady body-fitted hinge submodel with GCI):** accurate at the real gap, has a
   precedent, fits in 4 weeks.
2. **O4:** the same safety, with less hinge physics.
3. **O6:** cheap, and needed for correct hinge kinematics in every option.
4. **O1:** keep only Phase A/B as a blade-scale cross-check. The porous ear leaks 3–77x the gap flux at feasible K.
5. **O2:** blade-scale only.
6. **O3:** kinematics demonstration only. Its end-face inversions invalidate hinge flow.
7. **O5:** not feasible in this time or with this tool.

Literature found this session (pmc.ncbi.nlm.nih.gov was not reachable, so only abstracts and metadata were checked):
Simon et al. 2010, Ann Biomed Eng, doi:10.1007/s10439-010-0086-3; Simon et al., Ann Biomed Eng,
doi:10.1007/s10439-009-9857-0; Hester, Vasil & Burns, "Improving convergence of volume penalised fluid-solid
interactions", arXiv:1903.11914; Angot, Bruneau & Fabrie, "A penalization method to take into account obstacles
in incompressible viscous flows", Numer Math 81:497–520, 1999, doi:10.1007/s002110050401. The 2024 JFM Brinkman
force paper's authors need literature check.

---

## Round 1 — Numerical

### R1 — Numerical Analyst: checks of the Nov 10 handoff

Source: `bme298_debate/math_checks.py` (sections C1–C13), output in `bme298_debate/math_checks_output.txt`.
Assumptions: blood rho = 1060 kg/m3, mu = 3.5e-3 Pa s (Newtonian), h_gap = 50 um. Assumptions I made that are not in the handoff: leaflet thickness t = 0.7 mm nominal (1.1 mm is the upper bound if the planform is a half-disk of R = 10.8 mm); closure time 2–10 ms; Study 2 end time 0.21 s; memory scaling exponent b = 1.3–1.5. All literature points are marked "needs literature check" because I did not search the literature this round.

#### 1. Results table

| # | Check | Result | Verdict | Implication |
|---|---|---|---|---|
| C1 | sqrt(K), tau = rho K/mu | sqrt(K) = 3.2 / 10 / 17.3 / 31.6 um; tau = 3.0 us / 30 us / 91 us / 0.30 ms. The frozen step of 8.5 us is 2.8 tau | PASS | O1/O2: the time step is tied to tau, so K cannot go much below 1e-10 on the desktop |
| C2 | Darcy leakage through the ghost leaflet, u = K dp/(mu t) | Open leaflet (about 140 Pa normal load, back-calculated from torque): 6 cm/s at K=1e-9, a few % of the bulk velocity, which fits the −8.6% torque. **Closed valve at 80–100 mmHg: 4.4–5.4 m/s at K=1e-9 and 4–5 cm/s at K=1e-11 (about 2 L/min through both leaflets).** Keeping the leak under 1 cm/s needs K ≤ 1.8e-12 (tau ≈ 0.6 us) | **FLAG** | O1/O2 cannot model diastolic hinge leakage or closure with a ghost leaflet at any K that is practical. They are usable only for the open or opening phase |
| C3 | Exact 1D slot solution: 50 um gap with one Brinkman wall | Gap flux ratio vs. no-slip: 3.9x (K=1e-9), 1.9x (1e-10), 1.22x (1e-11), 1.06x (1e-12). Housing WSS ratio: 1.88 / 1.23 / 1.07 / 1.02 | **FLAG** | **O2 cannot resolve the real 50 um gap.** In O1 the gap flux and WSS must come from the Reynolds film, with the 3D pocket blocked |
| C4 | K-study design and Richardson extrapolation in sqrt(K) | The ratios of sqrt(K) are not constant (1.83, 1.73), so the generalized formula is needed. Three points fix T0, c and p exactly, leaving no residual. A ±0.4% numerical noise moves the observed p from 1.0 to 3.1. The mean element edge is 0.61 mm, much larger than sqrt(K), so the Brinkman layer is unresolved and the error is controlled by the mesh, O(h), not by O(sqrt K). Runtime goes as K^-0.63: about 2.8 h at 3e-10 and 5.7 h at 1e-10 | **FLAG** | O1/O2 Phase B: fix p = 1 and use the 3rd point as a check. Report it as a "K sensitivity", not a K convergence, unless the interface is refined locally to h ≤ sqrt(K) |
| C5 | Halo | +43.5% / +43.9% (claim +44%). Equivalent halo is 150–240 um per face, i.e. **3–5 gap widths**, about 1/4–2/5 of the mean element | PASS (number), **FLAG** (meaning) | The v22 ghost leaflet **filled the hinge gap**. Every K=1e-9 result near the hinge is void until FixChi is run. The +25.6% Fy fits a plate whose thickness grew by about 44% |
| C5b | Volume and DOF bookkeeping | 53,832.7 + 2×208.5 = 54,249.7 mm3; 0.77 DOF per element (P1P1 is consistent) | PASS | — |
| C6 | Body-fitted torque vs. angle | Monotonic. Linear R² = 0.999, slope −0.0099 mN m/deg, about −7.5%/deg. **A/B asymmetry is 1.0 / 0.6 / 1.8 / 2.7%.** At 84° the true value could be anywhere from 1.4% to 4.1% because only 2 significant figures are given | **FLAG** | The handoff's "within 2.2%" is wrong at 84°. A 5% bias gate cannot be tested on 2-sig-fig data. Report 4 significant figures |
| C7 | Statistics with n = 4 angles × 2 leaflets | A pure −8.6% bias gives R² = 1.000, so **the R² > 0.95 gate cannot discriminate.** Paired t at n=4 and sd=1%: minimum detectable effect 2.1%, CI half-width ±1.6%. Bland–Altman LoA each uncertain by ±2.8% (sd 1%) to ±5.5% (sd 2%). Leaflets A and B are mirror images, so the effective n is about 4, not 8. The BH example works mechanically but adds little here | **FLAG** | Phase A gate: replace "paired t + R²" with an equivalence check, i.e. the 95% CI of the bias falls inside ±5% (which needs bias ≲ 3.4% at sd 1%), plus the slope with its CI. Better still, use a V&V-20-style E vs. u_val built on the GCI (needs literature check) |
| C8 | Lubrication validity | eps = 0.025 means L = 2 mm. Systole: eps·Re = 0.008–0.06. Diastolic leakage at 80–100 mmHg: **eps·Re = 0.12–0.60** for L = 2…1 mm. Poiseuille WSS (h/2)dp/ds = 133–333 Pa, so **"~270 Pa" is reproduced (L ≈ 1–1.2 mm).** Couette WSS from ear motion (closure in 2–10 ms, r = 2 mm) is only 4–20 Pa. Unsteadiness h²ω/ν: 0.006 (cardiac), 0.08–0.38 (closure), **3.8 at impact (0.2 ms)** | PASS with caveats | O1 film is valid in systole and marginal for diastolic leakage, mostly at the mouths. It is **not valid at impact** (it needs the squeeze term, and even then it is unreliable). The 270 Pa belongs to the **closed-valve leakage phase**: almost all of it is pressure-driven, and during the closing swing itself dp is only a few mmHg |
| C9 | ALE shear and inversion | 0.67 per degree requires R = 1.92 mm, which matches the ear radius (PASS). At r = 2.1–2.3 mm the end faces slide 37–40 um per degree, i.e. 0.73–0.80 gap widths per degree. gamma = 1 across the full gap is reached at 1.25–1.36° (Laplace failed at 1.1°). The onset at 3° matches gamma ≈ 2.3, meaning Yeoh spread the slide over about 2.3 gap widths. One element of height h/n reaches gamma = 1 after 0.62–0.68° (n=2), 0.42–0.45° (n=3) or 0.31–0.34° (n=4). Note that uniform shear has det F = 1, so inversion comes from sliding nodes meeting pinned nodes at the corners | PASS (numbers) | O3: curing the end faces by remeshing needs a remesh every 0.3–0.65°, i.e. 20–45 remeshes over 14°. That is not realistic at about 5 h per out-of-core interval. Reducing the sliding instead (end-face slip with a gap-local mesh that rotates rigidly, O6 kinematics) is the only cheap fix |
| C9b | Tip displacement | 10.8 mm × 14° = 2.64 mm (arc) | PASS (numerically) | 10.8 mm is the **z-coordinate of the ear end face**, not the in-plane distance from pivot to tip. Read the tip radius from the CAD before using this as a pass/fail check |
| C10 | ALE run timing | Rates of 0.018–0.054 s/h (0.25–0.75°/h); mean 0.026 s/h; 65.4% complete. **Projected finish: Nov 10 between 08:40 and 17:30** | PASS | It has probably finished or stalled. Check the Log and the LinRes column |
| C11 | GCI and memory | Roache/Celik 3-grid function implemented (worked example: p = 3.1, GCI = 1.6%, asymptotic ratio 0.98). Refining from 1.58 M DOF: r = 1.3 gives 3.5 M DOF and **42–93 GB**; r = 2 gives 12.6 M DOF and 224–649 GB (direct solver, b = 1.3–1.5) | **FLAG** | On the desktop the GCI is only possible if the **current mesh is the finest grid** (coarsen to 0.72 M and 0.33 M DOF), or with HPC or an iterative solver. Global coarsening by r² = 1.69 can leave fewer than 2 elements across the 50 um gap, so run the GCI separately for the bulk and the gap |
| C12 | Added mass | Flat plate rotating about its edge: I_a = (9π/128)ρc⁴s_eff = 1.2–5.9e-8 kg m², so **I_a/I_s = 1.3–6** (more with confinement). The fluid inside chi has inertia equal to 56% of I_s | **FLAG** | Partitioned coupling is unstable without relaxation. COMSOL's fully coupled rigid-body solver is monolithic, so O3 and O4 are safe. O1/O2 Phase E is effectively partitioned, because chi comes from extrusion and has no Jacobian. It needs artificial inertia or Aitken relaxation **and** an internal-fluid momentum correction to the penalty torque (Uhlmann-type, needs literature check). Leaving that correction out biases the torque during acceleration by an amount of order 0.56·I_s·α |
| C13 | Miscellaneous | Torque differences −8.64% / −8.74% (claimed −8.6 / −8.8); Fy +25.6% PASS; ghost A/B 1.14% PASS. ALE runs only to 84°, 1° short of the 85° stop. Torque-expression signs (x+pin_x)/(y−pin_y) are UNVERIFIED, since I cannot see the script | minor | Check the sign convention in the B-leaflet torque before Phase A |

#### 2. Where the handoff is wrong or misleading
1. **A/B symmetry "within 2.2%" is wrong.** It is 2.7% at 84°, and the 2-significant-figure data can only bound it to 1.4–4.1%.
2. **The R² > 0.95 gate always passes.** Torque varies threefold across the angles, so any proportional bias still gives R² ≈ 1. A paired t-test also answers the wrong question: it tests "bias ≠ 0", not "|bias| < 5%".
3. **The "~270 Pa at closure" label is misleading.** It is the Poiseuille term at the closed-valve load (80–100 mmHg). The closing motion itself contributes only about 4–20 Pa of Couette shear.
4. **The handoff misses that the ghost leaflet cannot hold a closed-valve pressure** (C2), and that a penalized wall next to the gap multiplies the gap flux by 1.2–3.9 for K = 1e-11 to 1e-9 (C3).
5. **The K study will not converge on a fixed mesh.** sqrt(K) is 10–32 um while the elements are hundreds of um. Three K values give no redundancy for Richardson extrapolation.
6. **A GCI on 3 refined meshes does not fit in 28–32 GB.** The current mesh has to be the finest grid.
7. **Phase E's added-mass concern is real, but the larger risk is that the moving-chi torque has no internal-inertia correction and no Jacobian.**

#### 3. Implications per option
- **O1 (hybrid ghost + film):** The numbers support it for **systolic and opening** bulk flow, with gap flux and WSS from the film (C3 makes the film mandatory, not optional). Diastolic leakage is outside the ghost's reach (C2) but inside the film's, with eps·Re ≤ 0.6, provided the leaflet is treated as rigid and impermeable in the film model. Impact is outside both.
- **O2 (ghost alone):** Numerically disqualified for the gap at K ≥ 1e-11 (C3) and for the closed valve (C2). Fine as a bulk-flow tool.
- **O3 (ALE):** Kinematics and remeshing are proven. End-face inversion is geometric (C9) and will not go away with more remeshing at this cost. Coupling is monolithic, so added mass is not a problem.
- **O4 (two-scale):** The most robust numerically. The fixed or prescribed real-hinge model is body-fitted (no K error, no leakage), and a GCI on its local gap mesh is affordable. It needs a consistent hand-off of dp across the hinge, which C8 shows is the quantity that controls WSS.
- **O5 (overset):** No numbers favour it. It still needs at least 2–3 cells across 50 um in the overlap and roughly doubles the DOF, which already exceeds memory (C11).
- **O6 (pivot or sliding kinematics):** Cheap, and directly attacks the 37–40 um/deg end-face sliding (C9). It should come before any further ALE work.

#### Bottom line (ranking)
1. **O4**: body-fitted, monolithic, and GCI-feasible locally; the safest accurate result.
2. **O1**: the best route to novelty, but only for systole and diastolic leakage via the film, with the gates corrected (equivalence test, sig-figs, coarsening-based GCI).
3. **O6**: a cheap enabler for O3/O4 that explains the end-face inversion quantitatively.
4. **O3**: proven kinematics, but the end faces cannot be fixed by remeshing within the timeline.
5. **O2**: fails at the 50 um gap and at closure for any practical K.
6. **O5**: too expensive in DOF and memory, with no numerical advantage.

---

## Round 1 — Clinical

### R1 — Clinical / Physiology Expert

Scope: which outputs matter clinically for hinge thrombosis and blood damage, which cardiac phase dominates, what each option O1–O6 can credibly deliver, and whether the model inputs are physiologically correct.
Numbers come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/clinical/calc.py` (Python, order-of-magnitude slot and Couette estimates). They are not CFD.

#### 1. What the literature says matters (checked this session through search abstracts and snippets only)

I could not open full texts because pmc.ncbi.nlm.nih.gov and link.springer.com are blocked in this session. Each claim below rests on the abstract or snippet returned by search. Before the thesis cites any of them, someone must open the full text.

- **Simon HA, Ge L, Sotiropoulos F, Yoganathan AP**, "Simulation of the three-dimensional hinge flow fields of a bileaflet mechanical heart valve under aortic conditions", *Ann Biomed Eng* 2010 (DOI 10.1007/s10439-009-9857-0, https://link.springer.com/article/10.1007/s10439-009-9857-0). Micro-CT SJM hinge. The abstract reports two strong diastolic leakage jets and peak shear of **1920 dyn/cm² (192 Pa) in systole versus 6115 dyn/cm² (612 Pa) in diastole**, so diastole is about 3.2 times systole. The year in the DOI prefix (2009 vs 2010) still needs checking.
- **Simon HA et al.**, "Numerical investigation of the performance of three hinge designs of bileaflet mechanical heart valves", *Ann Biomed Eng* 2010 (https://experts.umn.edu/en/publications/numerical-investigation-of-the-performance-of-three-hinge-designs/). Gap width matters, and **recess/leaflet shape matters more**, for flow distribution and thromboembolic potential. Elevated shear is linked to hemolysis and platelet activation; stasis and recirculation are linked to poor washout and thrombus.
- **Jun BH, Saikrishnan N, Yoganathan AP et al.**, micro-PIV of SJM valves with three hinge gap widths, *Ann Biomed Eng* 2014 (PubMed 24976188, https://www.pubmed.gov/24976188/; author list and title still need checking). The widest gap gave the highest Reynolds shear stress (about 1000 N/m²) during **leakage flow**. The **narrowest gap gave the poorest washout and the most stasis**. The authors conclude the design must balance shear against washout.
- **Govindarajan V, Udaykumar HS, Chandran KB**, 2D hinge model of the closing phase (https://pmc.ncbi.nlm.nih.gov/articles/PMC2744151; exact title, year and venue still need checking). High shear combined with **long residence time in leakage flow during closure** can activate platelets.
- **Dumont K, Vierendeels J, Kaminsky R et al.**, ATS vs SJM Regent FSI, *J Biomech Eng* 2007 (https://researchconnect.stonybrook.edu/en/publications/comparison-of-the-hemodynamic-and-thrombogenic-performance-of-two/). During forward flow, **no platelets exceeded the 35 dyn·s/cm² stress-accumulation threshold**. During regurgitation, 0.81% (SJM) vs 0.63% (ATS) exceeded it, and the difference was attributed to hinge design.
- Hedayat & Borazjani, APS DFD 2017 abstract (https://archive.aps.org/dfd/2017/d4/5). Bulk-flow platelet activation was several times larger than hinge activation for a 25 mm SJM. This is an abstract, not peer reviewed, but it is a fair counterpoint to examiners who argue "the hinge is everything".
- The Bluestein stress-accumulation and power-law models, and Yun/Dasi lattice-Boltzmann hinge blood-damage studies, appeared in the search results but I did not check them. They need a literature check.

**Clinical conclusion:** four studies (Simon 2010, Jun 2014, Govindarajan, Dumont 2007) agree that hinge risk is dominated by **diastolic leakage plus the closure transient**, not by peak systole. The mechanisms are high gap shear (activation) and poor washout in the pocket (stasis and accumulation). Systolic hinge flow is mostly forward washout.

#### 2. Phase-by-phase magnitudes (my estimates)

Assumptions: h = 50 µm, slot length L ≈ 2 mm (from the handoff's ε = 0.025), ear radius R ≈ 1.92 mm (back-computed from the handoff's 0.67 shear per degree), μ = 3.5 mPa·s.

| Phase | Driver | γ̇ (1/s) | τ (Pa) | Comment |
|---|---|---|---|---|
| Systolic opening (15° in ~40 ms, UNVERIFIED timing) | ear wall motion, Couette | ~250 | ~0.9 | Benign. The pressure-driven part depends on Δp across the hinge, not estimated. |
| ALE test (14°/s) | Couette | ~9 | 0.03 | About 18× slower than physiologic, so kinematics only. |
| Closure (ω = 32–300 rad/s, UNVERIFIED) | ear Couette plus squeeze flow | 1.2e3–1.1e4 | 4–40 | Squeeze and impact transients not included. Simon's 612 Pa suggests much higher local values. |
| Diastolic leakage, Δp = 80–100 mmHg | Poiseuille slot | 3.8e4–4.8e4 | **133–167** | Re ≈ 5–6, transit ≈ 5–6 ms. Same order as the 612 Pa in Simon 2010 (real recess geometry, jets, corners). |

Damage metrics for **one** diastolic transit at Δp = 100 mmHg:
- Linear stress accumulation τ·t = **0.84 Pa·s**. The Hellums-type threshold is about 3.5 Pa·s (35 dyn·s/cm²), so about 4 transits or recirculating re-entries reach it. This is exactly why **residence time and washout in the pocket** are the clinically decisive coupling.
- Giersiepen power law: HI ≈ 0.13% and LDH-release index ≈ 0.38%. These are steady-calibration models applied outside their exposure range, so use them for ranking designs or angles only, never as absolute values.
- Leakage through four ear gaps (assumed 4 mm span each, UNVERIFIED): about 0.32 mL/s, or about 0.16 mL per diastole. This is a small share of total valve leakage, but it is the stream that washes the hinge.

#### 3. Physiologic correctness of the inputs

**Newtonian vs non-Newtonian.** With the Carreau (Cho–Kensey) model, viscosity is +0.7% above μ∞ at 4.8e4 1/s, +2.2% at 7.6e3, +8% at 1e3, +20% at 250, +37% at 100 and +57% at 50 1/s.
- Newtonian is acceptable **inside the gap** during leakage and closure (error under 2%).
- Newtonian is **not** acceptable for stagnation and washout in the pocket or recess, nor for diastolic stasis, where γ̇ is below 100 1/s and the error exceeds 35%.
- Recommendation: run Carreau as a sensitivity case for any residence-time or washout output. The advisor rule allows this because Carreau is more physiologic, not less.

**Angle range.**
- Abbott material and two SJM patents found this session give full opening of **84–85°** (Abbott Regent page: https://cardiovascular.abbott/int/en/hcp/products/structural-heart/surgical-valve-solutions/regent-valve.html; US6176876, https://patents.google.com/patent/US6176876). Opening is search-verified.
- The **closed angle** (about 25–30° from the orifice plane) is UNVERIFIED. No source gave it directly. A Govindarajan snippet mentions an SJM model with about 55° leaflet traverse, and 85 − 55 = 30°. That is an inference, not a source. Ask the advisor or the CAD supplier for it.

**The 70→85° study range excludes closure and the closed position.** This is the most important clinical point in this file. The modeled window covers only the late-opening and fully-open systolic plateau, which by §1–2 is the **lowest-risk** phase for the hinge. Torque at 70–84° is good mechanics validation, but it says almost nothing about thrombosis.
- A thesis titled around the "real 50 µm hinge gap" will be asked: *what happens in the gap in diastole?*
- The answer needs at least **one closed-angle geometry** with a diastolic transvalvular Δp: aortic about 80–100 mmHg (UNVERIFIED waveform source; confirm the valve position with the advisor), mitral about 100–120 mmHg peak systolic LV pressure.

**Gap fidelity.**
- The possible pivot offset (handoff §4) matters clinically. If true motion is rotation plus slide, the gap at the closed position may differ from 50 µm.
- Since τ_w ∝ h·Δp/L and leakage flux ∝ h³, a ±10 µm error changes leakage by about ±50–70%.
- O6 is therefore not cosmetic. It sets the clinically dominant input.

**Waveform and inertia.** The handoff Study 2 waveform peaks at t = 0.13 s. Its source is unverified and needs advisor confirmation. Leaflet inertia of 9.5e-9 kg·m² is fine for opening but irrelevant to diastolic leakage at a fixed closed position.

**Cavitation.** It needs about 95 kPa (≈ 713 mmHg) of local drop from atmospheric absolute pressure to vapor pressure. That happens only in sub-millisecond closure and impact transients. No option here is compressible or multiphase, or resolves leaflet impact. State it as out of scope and flag regions with p_abs < p_vapor as a cavitation-risk indicator only.

#### 4. Option-by-option: what each can credibly deliver

| Output | O1 ghost + Reynolds film | O2 ghost alone | O3 ALE remesh | O4 two-scale | O5 overset | O6 pivot opt. |
|---|---|---|---|---|---|---|
| Leaflet torque/force, 70–85° | Yes (≤9% bias now) | Yes | Not credible (flow OFF) | Yes (body-fitted) | In principle | No (kinematics only) |
| Leaflet-surface WSS | **No**. The diffuse interface (√K ≈ 32 µm, comparable to h) smears the wall. | **No** | No | Yes | Possibly | No |
| Gap WSS, uniform-gap part | Yes, **analytic-quality** in the film region, if h(x) is right | No. The gap is under-resolved and the penalized wall sits about √K from the true wall. | No | Yes (3D, fixed angle) | Only if the overlap resolves 50 µm | Indirectly (true h) |
| 3D recirculation and residence time in the pocket/recess | **No** inside the film. The 3D pocket is blocked or approximated, so washout is missed. | Poor near walls | No | **Yes** (fixed angle, Lagrangian tracers) | Maybe | No |
| Diastolic leakage jets (closed angle) | Only if the film runs at the closed angle and the mouths are treated in 3D | Weak (penalty leakage through a K-dependent "porous" leaflet contaminates the small leakage flux) | No | **Yes**, if a closed-angle body-fitted case is built | — | Sets h at closure |
| Closure transient | Not credible: ε·Re ≈ 0.24 plus a moving sharp χ gives noisy forces, and the film breaks down at the corners | Qualitative only | No | Prescribed-motion 2D/3D sub-model possible | Designed for this, but unproven in 4 weeks | No |
| Cavitation | No | No | No | No (indicator only) | No | No |
| Platelet SA / power-law damage | Gap only (film streamlines, analytic transit) | Unreliable (wall shear wrong) | No | **Yes** (pathline SA in 3D) | — | — |

Notes:
- **O3** (flow-OFF ALE) delivers no hemodynamics, and its mesh is invalid at the ear ends. It is a kinematics and meshing feasibility result only.
- **O5** is clinically the right tool for closure, but it is unproven in COMSOL and not realistic in 4 weeks.
- **O6** is the cheapest way to make the gap physiologically correct across angles.

New options:
- **O7, closed-position body-fitted leakage model.** A fixed-angle case at the closed angle with diastolic Δp, steady or a short transient. It is the single most clinically valuable addition. It reuses the working v20_fixedangles workflow and gives gap WSS, leakage flux, jets, pocket residence time and Lagrangian SA in the dominant phase. Risk: the CAD closed-stop angle is unknown, and rotating the leaflet into the seat may need a geometry fix.
- **O8, 2D hinge cross-section with prescribed closing motion** (Govindarajan-style). Cheap, so it can resolve the gap with ≥10 cells and capture closure squeeze shear. It is a sub-model, not a replacement.

#### 5. Minimum clinically meaningful deliverable set (publishable)

1. Validated fixed-angle torques at 70–84° with mesh GCI. This already exists, and the GCI is still missing.
2. **A closed-angle (O7) leakage case**: hinge gap WSS (max, area-weighted), leakage flow rate per hinge, jet description, and pocket residence time (tracer washout or mean age-of-fluid) with a Newtonian vs Carreau sensitivity.
3. **Particle-based stress accumulation** (linear τ·t plus one power law) along pathlines through the gap, opened vs closed, reported as a distribution and the fraction above 35 dyn·s/cm². Compare it against Dumont 2007's 0.81%, citing it as a published reference to compare against, not as a match.
4. O1 thin-film gap WSS versus 3D (O4/O7) in the uniform-gap region, as the methods-novelty piece.
5. A gap-width sensitivity (for example 40/50/60 µm, motivated by O6 and Jun 2014). This is cheap and clinically framed as the trade-off between shear and washout.

Out of scope, stated explicitly: impact, cavitation, and two-way FSI through closure.

#### Bottom line (clinical ranking)

1. **O4, extended with O7 (closed-angle leakage)**: the only path that credibly delivers diastolic gap WSS, pocket washout and SA in the phase that dominates hinge thrombosis.
2. **O1**: an elegant and novel analytic-quality gap WSS, but it cannot give pocket recirculation, leaflet WSS or closure. It is valuable only when validated against item 1.
3. **O6**: cheap, and it fixes the clinically dominant input (h at each angle, especially closed).
4. **O8 (2D closure sub-model)**: low cost, adds the closure-shear story.
5. **O2**: gives torque, but its wall shear and leakage are K-contaminated, so it is weak clinically.
6. **O5**: right in principle for closure, infeasible in 4 weeks.
7. **O3**: no hemodynamics and an invalid mesh at the ears, so it is a feasibility appendix only.

---

## Round 1 — Feasibility

### R1 — Feasibility Engineer (COMSOL 6.4, 4-core / 28–32 GB desktop, ~3.5 weeks)

Scripts: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/feasibility/sched.py` and `plan.py`.

#### 0. What the measured run times tell us
I took these numbers from the handoff and checked them in Python:
- **Ghost cost per step.** The K = 1e-11 run did 788 steps (0.0067 s of simulated time at dt = 8.5e-6 s) in 10 h, which is **45.7 s per time step** at 1.58 M DOF with MUMPS. The K = 1e-9 Study 2 run took 80 min, so it used about 105 steps over 0.13 s (mean dt about 1.24e-3 s). Step size limits the cost. The cost of a single step does not.
- **Penalty time scale.** τ = ρK/μ is 3.0e-4, 9.1e-5, 3.0e-5 and 3.0e-6 s for K = 1e-9, 3e-10, 1e-10 and 1e-11. In the 1e-11 run the step froze at about 2.8τ. With dt ≈ min(1.24e-3, 2.8τ), the full 0.13 s window would take **1.3 h (measured) / 6.5 h / 19.4 h / 194 h**. A 0.025 s restart window around peak flow (t = 0.11 → 0.135 s) would take **0.3 / 1.2 / 3.7 / 37 h**. K = 1e-11 is infeasible on any calendar. The K study has to use restart windows.
- **ALE (O3).** 0.654 s and 9.2° in 25.2 h, flow OFF. That is **0.026 sim-s/h, or 0.36°/h**, well below the 0.06–0.1 sim-s/h figure quoted in the task. Remeshing cost goes with the angle swept, not with simulated time. A real opening (15° in about 50 ms) still needs roughly 15+ remesh restarts, and those were 2–6 h each once out-of-core set in.
- **GCI meshes (refinement ratio r = 1.3, 3D).** Coarse 0.72 M DOF, fine 3.46 M DOF. With MUMPS, memory grows about N^(4/3) and time per step about N² (UNVERIFIED scaling rule), so the fine mesh needs **about 2.9× the memory and 4.8× the time per step**. The ALE run already went out-of-core on 28.7 GB. Assume the fine mesh runs out-of-core on the desktop (×3 slowdown, as seen in ALE intervals 5–6), or use HPC or an iterative solver (GMRES + SCGS/AMG) for it.
- **Added mass.** A flat-plate estimate (chord about 10 mm, span 21.6 mm, rotating about an edge) gives I_added ≈ 5.1e-8 kg m², which is **about 5.3× the leaflet inertia (9.5e-9 kg m²)**. The geometry here is assumed, so this is UNVERIFIED. With a ratio this far above 1, segregated or fixed-point coupling between the leaflet ODE and the fluid will diverge unless it is relaxed or stabilized. This is the core risk of Phase E.

**Main finding:** for the phased plan, compute is not the bottleneck. Engineer time is (section 3).

#### 1. Option-by-option implementation and risks

##### O1 Hybrid ghost leaflet + Reynolds thin-film hinge
**Steps:**
- **Ghost leaflet:** Laminar Flow, plus a Volume Force −(μ/K)·χ·(u − u_s) on all three domains.
- **Indicator χ:** General Extrusion back-rotation (source = all domains, sharp indA/indB labels; this is the FixChi fix).
- **Film:** Thin-Film Flow interface (CFD Module) on the pocket-wall boundaries. In 6.4 it is boundary-level only (Reynolds equation, dependent variable is pressure). Set the film thickness to h = d0(R⁻¹x), using a second General Extrusion from the CAD gap. Set the wall velocity to user-defined ω×r.
- **Mouth coupling:**
  - (a) 3D pressure → film: a Border / pressure condition on the film edges, using a linear extrusion of the 3D p onto the film edge. One-way coupling is enough for Phase C.
  - (b) Film edge flux → 3D: an inflow/outflow boundary on a cut "mouth" face in 3D. This needs a geometry cut and a blocked 3D pocket so the gap is not counted twice.
  - **UNVERIFIED:** whether 6.4 offers a built-in Thin-Film ↔ Laminar Flow multiphysics coupling. I found none. The search only found the bearing/MEMS-style wall and base coupling. Plan for manual weak or extrusion coupling.
- **FSI (Phase E):** a Global ODEs and DAEs node with I·θ'' = T_pen(u, θ) + T_stop(θ), with T_pen as an intop integral.

**Numerical risks:**
- **Missing Jacobian terms.** General Extrusion gives no ∂χ/∂θ Jacobian term (UNVERIFIED for the destination-map dependence; COMSOL documents `nojac` and Jacobian fill-in for nonlocal couplings, but I could not open the page to confirm moving-map derivatives). The fully coupled Newton solver therefore degrades to a fixed-point iteration on θ. With an added-mass ratio of about 5, that diverges unless I_added is put on the left-hand side and compensated explicitly on the right-hand side (artificial-inertia stabilization), or unless θ is updated with relaxation inside a segregated step.
- **χ noise.** A sharp moving χ jumps when element centroids cross the surface, which shows up as torque noise. Use a smoothed χ (tanh over about 1–2 elements, calibrated to volume) for motion.
- **Film Jacobian.** The h(x,t) extrusion has the same missing-Jacobian issue, but in one-way mode h is prescribed, so it does not matter there.

**Wall-clock estimate:**
- Phase A: 5.3 h
- Phase B: about 14 h ghost, plus about 30 h for body-fitted GCI (assumed)
- Phase C: < 1 h compute, 3–4 days of development
- Phase D: about 5 h, 3–4 days of development
- Phase E: about 10 h if it works, 4–6 days of development, success probability ≤ 50% (judgement)

##### O2 Ghost leaflet alone
Same as O1 without the film; this is the existing v22 model.
- **Risks:** with sqrt(K) ≈ 32 µm against a 50 µm gap, the hinge is *not* resolved. The local mesh needs about 5 elements across the gap, which is unaffordable in a fixed mesh. The gap flux will be wrong, but blade torque within 5–9% is plausible.
- **Cost:** the K study + GCI above, about 45 h total. **Feasible and low-risk, but it says nothing new about the hinge.**

##### O3 Body-fitted ALE with automatic remeshing
Moving Mesh (Yeoh), prescribed mesh velocity + Nitsche, Sliding on 18 faces, remeshing at 1.5. These settings are proven.
- **Blocker:** inverted elements on the ear end faces (z = ±10.8 mm) from about 3° onward. The likely fix is to give the 50 µm end gaps an Identity Pair or a separate swept ring domain that rotates rigidly. That is a geometry rebuild, 3–5 days, uncertain.
- **Cost:** with flow ON, expect 2–3× the flow-off rate (UNVERIFIED), which is ≥ 60–120 h for one opening. Two-way FSI adds the Global ODE plus added-mass instability.
- **Not feasible as a primary result.** Keep it as a kinematics/remeshing demonstration (the existing run).

##### O4 Two-scale
- **(i) Blade-scale model:** the ghost leaflet (O2) or the pin-joint v21 for the FSI opening history θ(t).
- **(ii) Hinge model:** the real-hinge body-fitted model at fixed angles (already exists, v20_fixedangles), with the 3D pressure from (i) as boundary conditions. Optionally, the Thin-Film Flow model in the same pocket for comparison.
- **Risks:** low, because every piece already runs. The weakness is one-way coupling (hinge → blade feedback is ignored, but that is small since the gap flux is ≪ the main flow; UNVERIFIED magnitude).
- **Cost:** reuses Phase A/B runs plus about 1 h per angle of hinge runs. **Highest probability of a finished, defensible result.**

##### O5 Manual overset
- **How it would be built:** two components (a rotating leaflet-body mesh and a fixed background mesh), with fringe interpolation through General Extrusion both ways and hole cutting by χ.
- **Why it fails here:** there is no built-in overset, so the Jacobian across the fringes is missing in both directions, and the overlap must still fit inside the 50 µm gap. DOF is about 2×, so out-of-core is certain.
- **Development time:** ≥ 3 weeks. **Infeasible in this calendar.**

##### O6 Pivot optimization / rotation + sliding
- **Geometry check:** sample the ear and pocket surfaces, rotate them about a candidate pivot, and minimize the maximum penetration and the clearance variance. Python (scipy.optimize) on exported STL/mesh points, or a COMSOL parametric sweep with `minop` distance probes.
- **Cost:** about 1 day, 2 h compute.
- **What it gives:** it explains the ear-face inversion, and it tests whether pure rotation is kinematically admissible. If it is not, the ALE needs a translational DOF (prescribed mesh displacement = R(θ)x + s(θ)). **Cheap, novel, do it.**

##### O7 (new) Quasi-steady stationary ghost solves for the K and mesh trend
- **What:** stationary solves at the peak-flow boundary conditions (t = 0.13 s) avoid the τ step limit entirely. About 10 Newton iterations × 46 s per solve, roughly 15–30 min each (UNVERIFIED that steady converges at Re near peak).
- **Use:** as a *screening* tool for the K → 0 extrapolation (Richardson in sqrt(K)). Confirm only the final pair with transient restart windows.

#### 2. Solver notes (all options)
- **Linear solver:** keep MUMPS and always check LinRes at about 1e-14. Do not use cuDSS: on the T4 it silently gave LinRes = 1.
- **Penalty problems:** fully coupled, Jacobian update "On every iteration" (the ALE fix also applies here).
- **K ≤ 1e-10:** consider excluding the penalty-region velocity from the time-error estimate, or setting a larger absolute tolerance there. Otherwise BDF tracks the τ transient (UNVERIFIED that this relaxes the step limit without hurting the torque; test against the K = 1e-9 run).
- **Fine GCI mesh:** run on HPC, or use iterative GMRES with an AMG/SCGS preconditioner on the desktop.

#### 3. Schedule arithmetic

##### Compute
Totals computed in plan.py:

| Phase | Compute (h) |
|---|---|
| A: 4 angles, sharp χ, full window | 5.3 |
| B1: K study 3e-10, 1e-10 × 2 angles (0.025 s restarts) | 10.0 |
| B2: ghost 3-mesh GCI at 70° (fine out-of-core ×3) | 4.0 |
| B3: body-fitted coarse + fine at 70° (medium assumed 2 h, UNVERIFIED) | 29.7 |
| C: film one-way, 4 angles | 0.8 |
| D: prescribed opening + film, 3 smoothing widths | 5.2 |
| E: two-way FSI, 3 stabilization settings | 10.0 |
| O6: pivot fit | 2.0 |
| **Total** | **67 h (134 h with a 2× rerun factor)** |

Capacity from Nov 10 to Nov 27 (18 days, then writing):
- One desktop at 50–60% utilization: 216–259 h.
- Two COE desktops (handoff shows coe-v1-rdp-5 and COE-V2-8 both in use; availability UNVERIFIED): 432–518 h.

**Compute fits on one desktop even with 2× reruns.** HPC is only needed to make the fine GCI meshes in-core (it saves about 20 h). It is not required.

##### Engineer time
This is the real limit.
- A: 1 day, B: 2 days, C: 3–4 days, D: 3–4 days, E: 4–6 days, O6: 1 day.
- **Total 14–18 working days**, against about 13 working days available (Nov 10–27, minus weekends and Thanksgiving Nov 26).
- **E does not fit with margin.** E is therefore a stretch goal. A+B+C+O6+O4 fits (about 8–9 days), and D fits if C passes by about Nov 19.

#### 4. Gantt (best feasible plan: A → B → C → O4, with O6 in parallel; D is conditional; E is stretch)

| Day | Desktop 1 (interactive) | Desktop 2 / HPC (batch, overnight) | Gate |
|---|---|---|---|
| Tue Nov 10 | Save the ALE log; FixChi + Torque check, VA/VB = 208.5 mm³ | Study 2 at 70°, sharp χ (1.3 h) | Volume within 1% |
| Wed Nov 11 | Script the angles 75/80/84; O6 pivot-fit script | Angles 75/80/84 (4 h) | |
| Thu Nov 12 | Bland–Altman / regression (Phase A); build coarse/fine meshes | K = 3e-10, 1e-10 restarts at 70° (5 h) | **A: bias < 5%, R² > 0.95** |
| Fri Nov 13 | O7 stationary K screening; O6 results | K study at 84° (5 h); ghost GCI (4 h) | |
| Sat–Sun Nov 14–15 | (light) | Body-fitted coarse + fine GCI (about 30 h, HPC preferred) | |
| Mon Nov 16 | GCI and K extrapolation analysis | Spare / reruns | **B: GCI < 5%** → thesis safe |
| Tue–Thu Nov 17–19 | Phase C: Thin-Film Flow on pocket walls, one-way p from 3D, compare with v20 resolved gap (flux, WSS) | Film runs at 4 angles (< 1 h) | **C: within 10%** |
| Fri Nov 20 | O4 assembly: blade θ(t) and loads + real-hinge fixed-angle hinge loads/WSS | Hinge runs at extra angles | Two-scale result done |
| Mon–Wed Nov 23–25 | If C passed: Phase D, smoothed χ + film, prescribed opening. Else: harden O4/GCI | D runs, 3 widths (5 h) | **D: mass balance < 1%** |
| Thu–Fri Nov 26–27 | Thanksgiving / buffer; freeze results. E only if D passed on Nov 24 | E attempt (10 h) | Code freeze Nov 27 |
| Nov 28 – Dec 3 | Thesis writing, figures, slides | none | Defense Dec 4–8 |

#### 5. Uncertainties (UNVERIFIED)
- How Thin-Film Flow connects to 3D Laminar Flow in 6.4 (no built-in coupling found), and how its Border node behaves.
- Whether General Extrusion destination-map derivatives reach the Jacobian when the map depends on a Global ODE variable.
- Whether excluding variables from the time-error estimate relaxes the τ step limit safely.
- The body-fitted reference run time and memory at the fine mesh.
- The added-mass ratio, which needs the real leaflet chord and pivot location.

Documentation found: the COMSOL 6.4 "Thin-Film Flow Interface (3D)" page (https://doc.comsol.com/6.4/doc/com.comsol.help.cfd/cfd_ug_fluidflow_thinfilm.12.03.html) and the 6.4 "Common Settings for Nonlocal Couplings" page (https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_definitions.21.104.html), both seen only as search snippets. doc.comsol.com was blocked in this session, so I could not read the full text. The wall/base velocity options I cite come from the 5.3 Fluid-Film Properties page (https://doc.comsol.com/5.3/doc/com.comsol.help.mems/mems_ug_fluidflow.07.07.html) and may differ in 6.4. No journal literature was checked: needs literature check.

#### Bottom line (feasibility ranking)
1. **O4 Two-scale.** Every component already runs, and it fits with days to spare. This is the guaranteed publishable core.
2. **O1 Hybrid (phases A–D).** Compute is trivial (about 25 h). Phases C/D are 6–8 days of development and the novel part. Phase E (two-way FSI) is a stretch because of the missing ∂χ/∂θ Jacobian and an added-mass ratio of about 5.
3. **O2 Ghost alone.** Cheap and low-risk (about 45 h for K + GCI), but it cannot resolve the 50 µm gap (sqrt(K) ≈ 32 µm). Treat it as the validation layer for O1/O4, not as the result.
4. **O6 Pivot optimization.** About 1 day and almost no compute. It explains the ALE inversion and adds a kinematic novelty point.
5. **O7 Stationary K screening.** A cheap accelerator for Phase B. Not a result on its own.
6. **O3 Body-fitted ALE.** More than 60–120 h per opening with flow, and the ear-face inversion is unsolved. Report it as a kinematics demonstration only.
7. **O5 Overset.** No built-in support, at least 3 weeks of development, and still needs an overlap inside the 50 µm gap. Infeasible before defense.

---

## Round 1 — Inventor

### R1 — Inventor: two new options (O7, O8) for the real 50 um hinge

All numbers below come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/inventor/calc.py` (plus one inline check).
Inputs: rho = 1060 kg/m3, mu = 3.5 mPa s, h = 50 um, ear radius from pivot R = 2.2 mm (where the inverted elements appear).
**UNVERIFIED (need CAD/model):** gap path length L_s ~ 1 mm, wetted gap width w ~ 6 mm, end-face disc radius 1.5 mm, closing peak Omega ~ 100 rad/s, opening and closing times of 40 and 30 ms.

#### 0. Scale facts that shape the rest

| Quantity | Value | Consequence |
|---|---|---|
| End-gap shear in current ALE, R*dtheta/h | **0.77 per degree, 11.5 at 15 deg** (576 um of slide across a 50 um gap) | No smoothing can fix this: the topology has to change |
| Gap viscous time h^2/nu | 0.76 ms | Gap flow is quasi-steady on the cycle time scale |
| Gap Womersley (1.2 Hz / 30 ms event) | 0.08 / 0.40 | Inertance negligible |
| Gap Re (h-based), 5 mmHg systole -> 100 mmHg closure | 0.6 -> 12; eps*Re 0.03 -> 0.60 | Linear (Stokes) in systole, mildly nonlinear at closure |
| Entrance-loss share, Ke*rho*U^2/2/dp (Ke = 1.5) | 0.2% -> 3.8% | A small quadratic correction is enough |
| Pressure-driven WSS (h/2)dp/L_s | 17 Pa (5 mmHg) -> 333 Pa (100 mmHg) | Same order as Simon et al. 2010 (192 Pa systole, 612 Pa diastole, reported as dyn/cm2) |
| Couette WSS mu*Omega*R/h | 1.3 Pa (closing mean) -> 15 Pa (Omega = 100) | Wall motion is a minor WSS source |
| Couette / pressure gap flux | 1.5% (mean closing) -> 14% (peak) | Motion matters only near impact |
| Hinge viscous friction torque / blade torque (0.2 mN m) | 1e-4 (mean) -> 1.1e-3 (peak) | **The hinge does not feed back on leaflet kinematics** |
| Hinge leakage, 4 hinges over 0.5 s at 80 mmHg | ~0.4 mL/beat (~0.5% of stroke volume) | The hinge barely perturbs the bulk flow |
| Fraction of the cycle with a moving leaflet | ~8% (about 70 ms of 830 ms) | Fixed-angle runs are exact for the other ~92% |

The last three rows are the key physical result. **One-way coupling (blade FSI -> real-hinge submodel) is justified quantitatively**, to about 1e-3 in torque and 0.5% in volume. That turns O4 from a "fallback" into a defensible, error-bounded method. It also tells us where the remaining effort should go: the ~70 ms of motion, and the hinge outputs (leakage, WSS, washout).

---

#### O7 (new): frozen-geometry moving-wall (FGMW) real-hinge submodel + calibrated hinge impedance

##### Math
Thin-film mass balance in the gap (lubrication order):

  div( h^3/(12 mu) grad p ) = div( h U_t / 2 ) + dh/dt,   with U_t = tangential wall velocity and dh/dt = relative normal wall velocity.

**FGMW:** on a body-fitted mesh frozen at angle theta_k, set the leaflet wall to the full rigid-body velocity u = theta_dot e_z x (r - r_pivot), normal component included (this is a transpiration condition). In the Reynolds balance above:
- The tangential part supplies the Couette term exactly.
- The normal part supplies dh/dt exactly to first order.
- For a rigid body, the closed-surface integral of (omega x r).n is zero, so mass is conserved.
- The only thing missing is the change of h itself inside a window: dh = (dh/dtheta)(theta - theta_k).

With e = the eccentricity between the ear clearance arcs and the true rotation axis, dh <= e*dtheta. Because pressure flux goes as h^3, the flux error is about 3*e*dtheta/h:

| e | dtheta = 0.25 deg | 0.5 deg | 1 deg | Largest dtheta for 5% error |
|---|---|---|---|---|
| 10 um | 0.3% | 0.5% | 1.0% | 4.8 deg |
| 50 um | 1.3% | 2.6% | 5.2% | 0.95 deg |
| 200 um | 5.2% | 10.5% | 21% | 0.24 deg |

**Exactness result:** where the clearance surfaces are surfaces of revolution about the pivot axis (e = 0), FGMW is exact to lubrication order. The ear end faces z = +/-10.8 mm are planes normal to the pivot (z) axis, so they always qualify. e comes directly from the O6 pivot fit, which makes O6 a required input rather than a side study.

The bulk error is a separate term: delta/L of the leaflet stump inside the submodel box, 35 um per 0.5 deg at r = 4 mm. It is controlled by keeping the box small and taking box boundary conditions from the blade FSI.

**Calibrated hinge impedance (idea a + d).** Per hinge and per angle:

  dp_h = R(theta)*(Q - C(theta)*theta_dot) + K(theta)*rho*Q|Q|/(2 A^2)

The inertance term L*dQ/dt is dropped because the Womersley number is 0.08-0.4; it should be kept as a fit term and shown to be insignificant. Companion WSS law:

  tau_max = (h_eff/2)*dp_h/L_eff + mu*R*theta_dot/h

Friction torque is D*theta_dot, with D ~ 2.2e-9 N m s. It is shown to be negligible rather than included.

Calibration:
- R and K: least squares on dp_h(t), Q(t) from the existing v20 fixed-angle runs (70/75/80/84 deg). Each run spans a range of dp over the waveform, so it gives many points per angle at no new solve cost. Only control surfaces at the two gap mouths need to be added.
- C: from two steady Stokes "unit" solves per angle (dp = 0 with theta_dot = 1, and dp = 1 with theta_dot = 0). This is justified by the linearity numbers above, and each solve takes minutes because there is no time stepping. Alternatively, C can be regressed from the FGMW runs.

Validation:
- Leave-one-angle-out on R(theta).
- A held-out FGMW run at 77.5 deg.
- Check against the FGMW sequence through closure: predicted vs 3D leakage volume and peak WSS, target 10%.

##### COMSOL sketch (6.4)
1. Build a hinge box submodel (about 6 x 6 x 4 mm around one pocket) from v20.
2. Use a parameter th_k in a Rotate geometry node and remesh per angle with a boundary-layer mesh in the gap. This is the same workflow as v20_fixedangles.
3. Set the leaflet wall to Moving wall with u = -thd*(y - yp), v = thd*(x - xp). theta_dot(t) comes from the blade FSI or the prescribed kinematics.
4. Take box inlet/outlet conditions from the blade-FSI solution via Interpolation functions exported on the box faces.
5. Hand over each window's solution as the next window's initial value. Export the field on a grid, then import it as an interpolation function; this is safer than cross-mesh `initsol` after a geometry rebuild (UNVERIFIED which of the two works in 6.4).
6. Run windows only for opening and closing: about 30 windows each at 0.5 deg. Use plain fixed geometry elsewhere.

There is no ALE, no remeshing, no inverted elements, and no penalization parameter K.

##### Novelty and risk
- Simon, Ge, Sotiropoulos & Yoganathan (Ann Biomed Eng 2010, 38(3):841-853, https://pmc.ncbi.nlm.nih.gov/articles/PMC2841271) already simulated 3D SJM hinge flow with an immersed-boundary method and prescribed, experimentally measured leaflet motion. Two-scale hinge simulation is therefore **not** new.
- What is new (needs literature check; my searches found no lumped hinge impedance law):
  - the error-bounded FGMW/transpiration formulation and its exactness for revolution surfaces;
  - a calibrated R/K/C(theta) hinge law that can be reused as a boundary condition in any valve model;
  - a quantitative proof that one-way coupling is valid.
- Risks:
  - Extracting dp_h and Q from the existing runs needs well-defined mouth surfaces.
  - Box-boundary sensitivity, so box size must be part of the convergence study.
  - The reference is not yet GCI-converged, same as for every option.
- **P(works in 4 weeks): ~65% overall.** About 85% for the impedance fit from existing fixed-angle data alone (a publishable minimum), and about 60% for the full FGMW closure sequence.

---

#### O8 (new): split-gap sliding interface, a rescue for O3 (idea b, reformulated)

##### Why not a true rotating frame
In a frame that co-rotates with the leaflet, the NS equations gain rho*(2 Omega x u + Omega x (Omega x r) + Omega_dot x r).
- Coriolis/viscous is 2*rho*Omega*l^2/mu: 0.01 (gap, closing mean) to 0.15 (gap, peak), but **1.3 to 15 at pocket scale (0.5 mm)**. The terms cannot be neglected.
- The butterfly pocket walls are not axisymmetric, so they would move in that frame and still need ALE.
- Rejected. Rigid-rotation ALE (sliding mesh) needs no fictitious forces at all.

##### The idea
The inversion happens because one 50 um mesh layer must absorb 576 um of in-plane slide. Split each end gap at its mid-plane, z = +/-10.825 mm (geometry to be confirmed), with a planar sliding identity pair:
- The housing-side half-layer (25 um) is a separate geometry object with a fixed mesh.
- On the leaflet side, the pair face has prescribed in-plane displacement equal to rigid rotation over the ear footprint and free in-plane sliding elsewhere.
- Laminar Flow "Flow Continuity" on the pair carries velocity and traction across the non-conforming meshes. This is the same mechanism COMSOL uses for rotating machinery. COMSOL docs found: Rotating Machinery, Laminar Flow (https://doc.comsol.com/6.4/doc/com.comsol.help.cfd/cfd_ug_fluidflow_single.06.055.html) and the blog "Deformed Mesh Interfaces: Rotations and Linear Translations" (https://www.comsol.com/blogs/deformed-mesh-interfaces-rotations-and-linear-translations/). Both note the identity pair, Form Assembly, and "destination side finer".

A plane normal to the rotation axis is a surface of revolution, so the slide is pure tangential motion with **zero distortion: the end-gap shear drops from 11.5 to 0**. The remaining distortion moves to the lateral gaps. There, a non-concentric arc gives a normal strain of e*theta/h = 0.05 / 0.26 / 1.05 at 15 deg for e = 10 / 50 / 200 um. Yeoh smoothing handles this for e <= 50 um but not for e = 200 um. The O6 fit again decides.

Concentricity check before any build: sample the pocket-wall and ear-side points in the CAD, fit a circle about the pivot axis, and report the radial scatter sigma_r. Requirement: sigma_r << h, about 5 um, along the thin-gap arcs. If this passes, the lateral thin arcs can also become cylindrical sliding pairs.

##### Risk
- The v20 single blood domain must be re-split into an assembly.
- The 25 um slab needs at least 3 elements across its thickness.
- The non-conforming pair mesh must be finer on the destination side.
- **The swept-footprint problem:** the fixed slab must cover the whole area the ear sweeps.
- **P ~ 40%**, at about 1.5 weeks of geometry and mesh work. It does produce the first fully body-fitted moving real-hinge FSI, which is the most impressive result if it works.

#### Idea (c) alone and other discards
- (c) without (a) is just O4's prescribed-motion part, so I merged it into O7, where its derived error bound is the contribution.
- A literal rotating reference frame is rejected (Section O8).
- Full two-way coupling of the hinge impedance into the blade FSI is unnecessary (friction torque ratio 1e-3). It is optional for completeness.

#### Bottom line (ranking)
1. **O4 + O7**: highest probability of an accurate real-hinge result. The one-way coupling is now quantitatively justified, and FGMW plus the calibrated impedance law add publishable novelty.
2. **O6**: cheap, and it supplies the eccentricity e that sets the error of both O7 and O8. Do it first.
3. **O1**: most novel, but the film is invalid at mouths and corners, and moving sharp chi gives noisy forces and has no Jacobian. Keep it as a stretch goal, and validate it against O7.
4. **O2**: needed as the bulk-flow baseline and for the K/GCI study, but on its own it does not resolve the 50 um gap (sqrt(K) ~ 32 um).
5. **O8**: a targeted fix for O3's end-face inversion; high payoff, ~40% probability.
6. **O3 as is**: proven kinematics, but the mesh at the ear ends is invalid, so its flow there cannot be used.
7. **O5**: no built-in support in COMSOL and about 2x the DOF; leave it for follow-up work.

##### Sources checked this session
- Simon HA, Ge L, Sotiropoulos F, Yoganathan AP. "Simulation of the three-dimensional hinge flow fields of a bileaflet mechanical heart valve under aortic conditions." Ann Biomed Eng 2010;38(3):841-853. https://pmc.ncbi.nlm.nih.gov/articles/PMC2841271 (abstract-level facts only).
- "High-resolution measurements of leakage flow through bileaflet mechanical heart valve hinges" (Klusak & Quinlan; year and venue not verified). https://researchrepository.universityofgalway.ie/entities/publication/4dc93b67-aa9f-4bc9-9f78-d7100adac0fc. It reports that Simon et al. estimated a hinge Re of about 150. That definition differs from my h-based Re of about 10, and the difference needs a literature check.
- The COMSOL 6.4 documentation and blog listed under O8.
- Transpiration BC background: NASA report https://ntrs.nasa.gov/api/citations/19970016811/downloads/19970016811.pdf (only the search summary was read; needs literature check before citing it in the thesis).
---

# ROUND 2 — Rebuttals

## Canonical option labels for Round 2 (orchestrator)
Round-1 agents used "O7"/"O8" for different ideas. From Round 2 on, use these labels:
- O1 Hybrid ghost leaflet + Reynolds thin-film hinge
- O2 Ghost leaflet alone
- O3 Body-fitted ALE + automatic remeshing
- O4 Two-scale (blade FSI + separate real-hinge fixed/prescribed model)
- O5 Manual overset mesh
- O6 Pivot optimization / rotation+sliding kinematics
- N1 Closed-angle body-fitted diastolic-leakage case (Clinical's "O7")
- N2 Real-hinge box submodel: frozen-angle body-fitted meshes with full rigid-body wall velocity + fitted lumped hinge
       pressure-leakage law (Inventor's "O7"; merges Challenger's "O7" = O4+O6 + body-fitted hinge submodel)
- N3 2D hinge cross-section with prescribed closing motion (Clinical's "O8")
- N4 ALE end-gap mid-plane sliding split / Flow Continuity pair (Inventor's "O8")
- N5 Stationary peak-flow solves for cheap K screening (Feasibility's "O7")

Known numeric conflicts to resolve in Round 2:
- K-study run time: Feasibility (0.13 s window: 6.5 h @3e-10, 19 h @1e-10) vs Numerical (~2.8 h @3e-10, ~5.7 h @1e-10, K^-0.63).
- Penalized-leaflet leakage: Challenger (29.7x Poiseuille flux at K=1e-9, ~1300 mL/s at closure) vs Numerical (3.9x gap flux; 4.4-5.4 m/s through leaflet at 80-100 mmHg).
- Halo thickness: Advocate ~200 um (0.9 mm leaflet assumed) vs Numerical 150-240 um.
- Added inertia ratio: 1.3-6x (Numerical), 1.5-6.3x (Challenger), ~5.3x (Feasibility).

---

## Round 2 — Advocate

### R2 — Advocate rebuttal: O1 revised to "O1-R", a validated Reynolds-film hinge law

Scratch scripts: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/advocate/r2.py`. My R1 numbers come from `calc.py` in the same folder. Labels follow OPTION_MAP_R2.

Summary of this round: most of the Round-1 attacks are correct, and they are aimed at the ghost-leaflet half of O1. The Reynolds-film half survives. It also turns out to be the physics-based version of Inventor's N2 hinge law. So the case I now make is for O1-R, which keeps the film and gives up the ghost-film coupling.

#### 1. Challenger: "29.7x Poiseuille flux at K = 1e-9; ~1300 mL/s through the closed leaflet"

**Conceded, with a reconciliation.** I re-solved the 1-D composite slot: a 50 um gap next to a 300 um Brinkman layer. The gap carries 3.63x / 1.70x / 1.20x / 1.06x the no-slip flux and the porous layer adds 26.1x / 2.97x / 0.31x / 0.03x, for K = 1e-9 / 1e-10 / 1e-11 / 1e-12. That gives 29.7x total at K = 1e-9, exactly Challenger's number.

Numerical's "3.9x" is the gap-only part of the same physics (3.63x in my configuration). The two figures are consistent, not in conflict.

Darcy transleaflet velocity at 80 mmHg for a 0.7 mm leaflet: 4.35 m/s (K=1e-9), 0.44 m/s (1e-10), 4.4 cm/s (1e-11). This reproduces Numerical C2.

**Consequences for O1:**
1. The ghost leaflet must never be used to block the pocket, and never at the closed angle. My R1 §4 step 4 offered "penalize as solid" as one option for blocking. **That option is withdrawn.**
2. Wherever the film is active, the pocket is removed geometrically, and the ear-to-housing gap exists only in the film.
3. The leaflet is impermeable in the film by construction: the Reynolds BC on the ear is u·n = wall velocity.
4. Closure and diastole are taken **out of the ghost model altogether** (see §3).

#### 2. Numerical

**"The R² > 0.95 gate always passes."** Conceded. I replace it with an equivalence test (TOST): the 90% CI of the bias must lie inside ±5%. With n = 4 effective pairs the CI half-width is 1.18% at sd = 1% and 2.35% at sd = 2%. The test therefore passes only if |bias| < 3.8% (sd 1%) or < 2.6% (sd 2%). I also add the slope CI. Torques are to be reported to 4 significant figures.

**"The K study plateaus, because mesh 0.61 mm >> sqrt(K)."** Conceded. I checked local refinement. A shell around both leaflet faces (about 930 mm², UNVERIFIED area) costs about 31 M tets at 30 um elements and about 557 M at 10 um. Neither fits on the desktop. So the ghost result becomes a **K sensitivity at fixed mesh** with p fixed at 1, reported honestly as such. It is no longer a convergence claim.

**"The GCI does not fit in memory."** Conceded. For the bulk model, the current mesh is the finest grid and two coarsened grids are added. The real-gap GCI moves to the small N2 hinge box, where it is cheap (Numerical C11, Inventor).

#### 3. Clinical: "70–85° excludes the high-risk diastole and closure; O1 gives no pocket washout"

**Both points conceded.**
- Simon et al. 2010 reports, at abstract level as Clinical checked, a diastolic peak shear of 612 Pa against 192 Pa in systole.
- Our film estimate at 80–100 mmHg is 133–333 Pa (R1 table). The film is valid there: eps·Re is 0.12 at L = 2 mm and 0.60 at L = 1 mm.
- So the film's natural operating point is the closed-valve leakage phase, which is exactly the phase Clinical says matters. The ghost leaflet cannot be in that phase (§1).
- The film is therefore driven by mouth pressures from a **body-fitted closed-angle model (N1)**, not from the ghost model.
- Pocket recirculation and residence time are 3D quantities, and the film cannot give them. They come from N1/N2 with tracers or age-of-fluid.

The film contributes three things:
- WSS and flux maps in the gap.
- The h³ sensitivity. With h = 40/60 um, the flux changes by a factor of 0.51 or 1.73.
- The eccentric-gap correction. For a gap whose centre is offset by e (as a fraction of the mean gap), the flux grows by (1 + 1.5e²): 1.06x, 1.38x and 1.96x for e = 0.2, 0.5 and 0.8. The film evaluates this from the CAD clearance map without any 3D solve.

**Washout belongs to N1/N2. Gap shear vs. gap width belongs to the film. The two are complementary.**

#### 4. Feasibility: "14–18 engineer-days against about 13 available"

**Conceded.**
- Phase E (two-way ghost FSI) is **dropped**. Its missing dχ/dθ Jacobian and added-mass ratio of 1.3–6 do not fit in this timeline.
- Phase D (ghost + film, coupled during motion) is also dropped.
- What remains for O1-R is Phase C: about 3–4 days of film development (Feasibility's estimate), run on geometry and boundary conditions that N1/N2 produce anyway.

#### 5. Inventor: "N2, and the hinge barely back-couples, so one-way O4 is justified"

**Agreed.** My R1 gave 0.2% for the hinge torque relative to the blade, and Inventor's number is 1e-4 to 1.1e-3. One-way coupling is justified. This also removes the main reason for a two-way ghost-film coupling.

N2's lumped law is dp_h = R(θ)(Q − Cθ̇) + Kρ Q|Q|/(2A²). The Reynolds film supplies its coefficients from first principles:
- R(θ) = 12μ ∫ ds / (w h(s,θ)³). The uniform case is 5.6e10 Pa·s/m³, i.e. Q = 0.19 mL/s per hinge at 80 mmHg (L = 1 mm and w = 6 mm UNVERIFIED).
- C(θ) comes from the Couette and squeeze terms.
- K is the only empirical term, and Inventor estimates the entrance loss at 0.2–3.8% of dp.

So **N2 fits R from 3D runs, while O1-R predicts R from the CAD clearance map without a 3D solve and validates it against N2.** Agreement within 10% would turn a fitted law into a predictive one that can be reused at any gap width, eccentricity or θ (including the 40/50/60 um sensitivity) without remeshing a 50 um gap. To my knowledge that predictive capability is the novelty (needs literature check; no search this round or last found a Reynolds-film hinge law for a mechanical heart valve).

#### 6. Challenger's other points, briefly

- **Mouth pressures taken from the least accurate ghost region:** resolved. They now come from body-fitted N1/N2 runs.
- **Sharp χ validated vs. smoothed χ used:** moot. The ghost leaflet is now used only at fixed angles with sharp χ, as a blade-scale cross-check, with no moving χ.
- **Opposite-signed Tz/Fy errors:** still open. It is now a Phase A item. FixChi has to be run, followed by a traction-based check on a surface outside χ. If |ΔFy| stays above 10% after FixChi, the ghost model is reported as unsuitable for forces. That would not affect O1-R.

#### 7. Revised O1 ("O1-R"): scope, phases, gates

**Scope.** Real-gap hinge flux and WSS from a Reynolds film on the CAD clearance map. It is validated against body-fitted 3D (N1 at the closed angle under diastolic dp; N2 frozen-angle boxes in systole) and then used for gap-width and eccentricity sensitivity. The ghost leaflet is optional and limited to a systolic, fixed-angle blade cross-check.

| Phase | Work | Gate |
|---|---|---|
| A′ | FixChi and ghost fixed angles 70/75/80/84 (optional cross-check) | TOST 90% CI inside ±5%; Fy within 10% |
| C1 | Film on CAD h(s); mouth dp from N1 (closed angle, 80–100 mmHg) and from N2 at 2 systolic angles | Film Q and area-mean gap WSS within 10% of body-fitted 3D in the uniform-gap region; eps·Re map reported |
| C2 | Film-predicted R(θ), C(θ) vs. N2-fitted values | Within 10–15%; otherwise the film is reported as a bounded approximation |
| C3 | Gap-width (40/50/60 um) and eccentricity (from O6) sensitivity using the film only | None (an output) |

Dropped: D, E, the K-convergence claim, and the ghost-film coupling. Fallback: if C1 fails, N1/N2 stand alone and the film is reported as a negative result with its validity map. That is still publishable as a methods note.

#### 8. Primary deliverable or add-on?

**Add-on.** The primary deliverable should be **O4 + N2 + N1**: body-fitted, monolithic, local GCI, and covering diastole. **O1-R is the novelty layer on top**: the reduced, predictive hinge law. I now consider it the strongest novelty claim the thesis can make, but the thesis must not depend on it. This is a change from my R1 position.

#### Bottom line (updated ranking)

1. **O4 + N2 (+ N1)**: an accurate, body-fitted real-gap result that covers diastole, with a GCI that fits on the desktop.
2. **O1-R**: a predictive Reynolds hinge law validated against #1, and the main novelty; cheap because it reuses N1/N2 boundary conditions.
3. **O6**: about 1 day, and it supplies the eccentricity e needed by O1-R, N2 and N4.
4. **O2 (fixed-angle, systole only)**: a blade-scale cross-check; it cannot be used in the gap or at closure (29.7x leak, 4.4 m/s Darcy).
5. **N5**: a cheap K-sensitivity screen, useful only if O2 is kept.
6. **N4**: a high-payoff rescue for O3 (about 40% chance of success); a stretch goal.
7. **O3**: kinematics only. **O5**: future work.

---

## Round 2 — Challenger

### R2 — Challenger rebuttal

Labels follow OPTION_MAP_R2. New numbers come from
`/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/challenger/r2.py`.
Blood: mu = 3.5 mPa s; gap h = 50 um.

#### 1. Reconciling 29.7x (mine) with 3.9x (Numerical, C3)

Both numbers are correct. They are fluxes through **different control volumes**:

| Component, K = 1e-9 | Flux / true no-slip slot flux |
|---|---|
| Fluid gap only (0 < y < h) | 3.63 (same value in my finite-slab model and in the semi-infinite analytic) |
| + decaying Brinkman skin inside the penalized wall (~sqrt(K) deep) | +0.29, giving **3.93 = Numerical's C3** |
| + far-field Darcy flow through the penalized slab, K*T/(h^3/12) | +26 at T = 300 um, giving **29.7 = my R1**. It is 64x at T = 700 um and grows linearly with T |

Numerical deliberately moved the Darcy part into C2, the "leak through the leaflet". I counted it as flow
crossing the hinge cross-section.

Which number matters:
- **For gap WSS and the slot flux seen by the 3D model, the right figure is 3.9x.** The housing WSS is 1.9x.
  This applies to O2, and to any 3D hinge region O1 keeps.
- **For hinge mass balance (O1 Phase D) and for the mouth pressures fed to the film, the Darcy shunt in parallel
  with the gap matters.** At K = 1e-10 the shunt alone is still 3–7x the gap flux for T = 300–700 um.

I withdraw the implication that the film itself is contaminated. As Numerical notes, the film treats the ear as
impermeable. My closure-leak figure also agrees with Numerical's: 1298 mL/s per leaflet at 80 mmHg rescales to
194.7 L/min for both leaflets at 100 mmHg, against their 194.5 L/min.

#### 2. Responses to the Advocate

**(a) Falsifiable Phase A gate.** I credit the willingness to fail. As written, though, the gate does not
discriminate:
- R² > 0.95 passes even under a pure −8.6% bias (Numerical C7).
- The gate tests **blade** torque, a quantity O1's hinge claim does not depend on.

Amend it as follows:
- Pass only if the 95% CI of the bias lies inside ±5% (an equivalence test, which needs a bias of about 3.4% or
  less).
- Include **FyA** in the gate, not just Tz.
- Add a traction-based force check on an offset surface.
- Add one hinge-relevant quantity: the **mouth-to-mouth Δp** across one hinge, ghost vs body-fitted.

The last item is what the film consumes.

**(b) "Hinge torque is about 0.2% of blade torque, so one-way coupling is justified."** I accept the number.
Inventor's 1e-4 to 1e-3 agrees with it. But it **undercuts O1's Phases D–E**. If the hinge does not feed back on
the kinematics, nothing is gained by running the film inside a moving-χ ghost FSI. The film only needs θ(t),
θ̇(t) and the mouth pressures, and a body-fitted model supplies those without porous leakage. The Advocate's own
coupling argument reduces O1 to "O4 with a film as the local model". That is the plan I am proposing.

**(c) Scale split.** I agree with the principle. The split fails at its seam:
- The film's WSS is dominated by the pressure term (h/2)dp/ds. Couette shear is only 4–20 Pa (C8).
- So hinge WSS is only as good as the mouth Δp.
- O1 takes that Δp from the 3D side, in cells that are 150–240 um wide (the halo, C5), next to a penalized ear
  with a Darcy shunt. The split is sound only if the 3D side near the hinge is **body-fitted**.

**(d) Novelty.** "First Reynolds-film hinge model coupled to a leaflet simulation" is plausible. It needs a
literature check; nobody in Round 1 found a prior example. My alternative keeps that novelty, as a film model
*validated against body-fitted 3D at 50 um*. The ghost leaflet adds only "we also penalized the blade", which
earns no novelty credit and is the source of most reviewer objections.

#### 3. Attacking my own alternative (N2 = O4 + O6 + frozen-angle body-fitted hinge)

**(i) Novelty.** Inventor is right. Simon, Ge, Sotiropoulos & Yoganathan already simulated 3D SJM hinge flow
under prescribed motion and boundary conditions from bulk FSI:
- Ann Biomed Eng 2010, doi:10.1007/s10439-009-9857-0.
- The three-hinge-design study, Ann Biomed Eng 38(11):3295–3310, doi:10.1007/s10439-010-0086-3, which varied
  the gap.

Jun et al. (2014; details need checking) measured three gap widths with micro-PIV. **Two-scale prescribed-motion
hinge CFD is therefore not novel.** What N2 can still claim:
1. The real CAD 50 um gap. Clinical reports roughly 150 um in Simon 2009; this needs a full-text check.
2. A local GCI at the gap, which the abstracts I found do not mention.
3. The error-bounded frozen-geometry transpiration formulation.
4. A calibrated R/K/C(θ) hinge law.
5. The film as a validated reduced model.

Honest verdict: this is a **methods-and-verification paper** (CMBBE or J Biomech Eng tier), not a discovery
paper. Reviewers will call it incremental unless items 2–5 are carried through.

**(ii) Frozen-angle transpiration error 3eΔθ/h.** The window size Δθ depends entirely on the eccentricity e,
which nobody has measured:

| e | Max Δθ for 5% error | Windows for opening (15°) | Windows for full swing (55°, closed angle UNVERIFIED) |
|---|---|---|---|
| 10 um | 4.8° | 4 | 12 |
| 50 um | 0.95° | 16 | 58 |
| 200 um | 0.24° | 63 | 231 |

Butterfly recesses are **not designed as surfaces of revolution**. Where the ear meets the recess wings, e of
100 um or more is plausible (UNVERIFIED). In that case N2 needs hundreds of remeshed geometries. The bound holds
only where h is about 50 um; it is looser where h is larger. That helps, but it is unproven. **Until O6 measures
e, N2's cost is unknown.**

**(iii) Meshing cost.**
- Estimate for a 6×6×4 mm box with 5 or more prisms across the gap: about 1.9e5 elements and 1.4e5 DOF
  (UNVERIFIED gap area of 12 mm²).
- At 0.25–1 h per remesh plus solve, that is 4–16 h for 16 windows and 57–229 h for 231 windows.
- A 3-grid GCI at 4 angles adds 6–24 h.
- Compute is fine for e ≤ 50 um. The risk is **boundary-layer meshing failures in a sliver gap at some angles**,
  which is the same class of failure that ended O3, and the need to automate the geometry rotation, remesh and
  field hand-off. Inventor marks the cross-mesh initial-value hand-off as UNVERIFIED in 6.4.
- The box faces need conditions from a global model. If that model is the ghost, my own §2(c) objection applies
  to N2. Box faces must sit in the bulk, several mm from the gap, and must be checked against the full body-fitted
  v20 model.

**(iv) The closure-phase gap.** This is N2's real hole:
- The global reduced ODE uses a torque table covering 70–84° only. There is **no closing-phase torque under
  reverse flow** and no closed-stop angle (UNVERIFIED).
- θ(t) through closure would have to come from literature kinematics (needs literature check).
- The frozen-angle approach fails at impact, where h²ω/ν ≈ 3.8 (C8).

So N2 delivers the opening and plateau (low clinical risk) and a closed-angle steady leakage case. It does not
deliver the closing transient. **No option on the table delivers closure dynamics in 4 weeks, O1 included.** The
thesis should state this.

#### 4. Clinical (N1, N3) and Feasibility

- **N1, closed-angle body-fitted leakage:** I endorse it as the highest-value addition. It is steady and reuses
  the v20_fixedangles workflow, and it targets the phase that four studies identify as dominant.
  - Risks: the unknown closed angle, and the leaflet–seat contact. Zero-gap regions cannot be meshed, so a minimum
    clearance has to be imposed and reported as a sensitivity.
  - Leakage scales as h³: ±10 um changes it by about −50% to +70%. That makes the O6 clearance at the closed angle
    a dominant input.
  - Run Carreau there, as Clinical says.
- **N3, 2D closing cross-section:** cheap and resolvable, and it gives the squeeze-shear story. But Simon et al.
  (three-designs abstract) report that recess shape matters more than gap width, and 2D cannot represent the
  butterfly recess. Use it as a supplementary sensitivity only. If the 2D ear is concentric, a rotating sliding
  interface avoids O3's inversion.
- **Feasibility:** I accept the conclusion that engineer time, not compute, is the binding constraint: 14–18 days
  needed against about 13 available. My own plan costs O6 (1 day), N1 (3), N2 box and sweep (4), film vs N2 (3)
  and local GCI (2), which is **13 days, i.e. zero margin**. The fix is to restrict N2 to 4 frozen angles plus the
  closed angle (no motion windows) unless O6 shows e ≤ 50 um. Motion windows become the stretch goal, replacing
  O1's Phase E.

#### Bottom line (updated ranking)

1. **N2-lite + N1** (body-fitted real-hinge box at 70/75/80/84° plus the closed angle, local GCI, film validated
   against it): accurate at 50 um, with clinically relevant diastole. The novelty is methodological, so it must
   be argued carefully.
2. **O6:** one day of work that sets e (and therefore N2's cost) and h at closure (leakage ∝ h³). Do it first.
3. **O4:** the safe core that N2 refines.
4. **O1, Phases A–C only:** a blade-scale cross-check plus the film. D–E add no physics, given the 0.2% coupling
   argument.
5. **N3:** a cheap closure-shear supplement; 2D misses the recess.
6. **N5:** a useful accelerator for whatever K screening survives.
7. **O2:** bulk flow only (gap flux 3.9x, closed-valve leak).
8. **N4:** high payoff, but about 40% probability and 1.5 weeks.
9. **O3:** kinematics demonstration only.
10. **O5:** infeasible in this calendar.

#### The one decisive experiment that would change my mind

Run the **FixChi sharp-χ ghost at 70° with local refinement to h_mesh ≤ sqrt(K) around one hinge, at K = 1e-10**.
Compare the **mouth-to-mouth hinge Δp** and FyA with the body-fitted v20 70° run on two meshes. If both agree
within 5%, with the equivalence CI inside ±5%, the ghost can feed the film honestly. I would then rank O1 Phases
A–D above N2, because it covers motion without hundreds of remeshes. If the Δp is off by more than 10%, O1's hinge
WSS is unusable and N2 stands. The run costs about 6 h of compute (Numerical C4 estimate for K = 1e-10) plus one
day of mesh work.

---

## Round 2 — Numerical

### R2: Numerical Analyst verdicts on the Round-1 numbers

Source: `bme298_debate/math_checks.py`, sections R2-C1 to R2-C8, appended to the Round-1 checks. The whole file runs cleanly, and its output is in `math_checks_output.txt`. Assumptions are as in R1. I did not check any literature this round; every literature item still needs a check.

#### 1. Leakage conflict (Challenger 29.7x vs Numerical 3.9x): no conflict (R2-C1)

I rebuilt both setups in one routine. The routine models a 1D fluid gap with no-slip at the housing, sitting next to a penalized layer of thickness T, and splits the flux into three parts. All values are flux ratios against the true no-slip slot.

| K | Flux in the gap itself | Decaying Brinkman layer | Darcy core inside the "solid" (T = 0.3 mm) | Challenger total | Numerical R1 |
|---|---|---|---|---|---|
| 1e-9 | 3.63 | 0.29 (semi-infinite) / −2.74 (wall-backed) | 28.8 | **29.7** | **3.93** = gap + layer |
| 1e-10 | 1.70 | 0.18 / 0.09 | 2.88 | 4.67 | 1.88 |
| 1e-11 | 1.20 | 0.02 | 0.29 | 1.51 | 1.22 |

Challenger's 3.63 / 26.1 / 29.7 is reproduced exactly. They are three different leak paths:
- **(a) Slip-enhanced flux in the gap.** This is the 3.6x, about √K of effective slip. It sets the gap velocity and the housing WSS (1.9x). It is my number.
- **(b) Darcy flow parallel to the gap through the penalized ear body.** This equals 12KT/h³: 26x for a 0.3 mm ear and 67x for 0.7 mm. It is Challenger's extra flux.
- **(c) Darcy flow through the leaflet thickness, driven by the transvalvular dp.** This is 4.35 m/s, or 1,297 mL/s per leaflet at 80 mmHg (t = 0.7 mm). Challenger's 1,298 mL/s is the same calculation.

Which one applies where:
- **Hinge leakage** is (a)+(b), so **Challenger's ~30x is the right figure for the hinge.**
- **Gap WSS** is (a).
- **Valve closure** is (c).

Challenger's limit of K ≤ 8.7e-13 for a leak under 5% corresponds to T = 0.6 mm. I get 1.7e-12 for T = 0.3 mm. Both are PASS, and both are about 1000x below what can be run. **The verdict stands: no ghost-leaflet number may be used inside the gap or at closure.**

#### 2. K-study run time (R2-C2): I withdraw my R1 estimate

My R1 power law (K^-0.63) assumed the K = 1e-11 run covered 0 to 0.0867 s with the time step frozen throughout. That would mean 10,200 steps in 10 h, or 3.5 s per step, which is implausibly cheap for a MUMPS factorization at 1.58 M DOF.

Feasibility's reading is self-consistent: the window starts at 0.08 s, so the run made 788 steps at 45.7 s per step, and the K = 1e-9 run gives the same per-step cost (about 105 steps in 80 min). The difference between the two estimates is the data point each one extrapolates from:
- Mine used the total wall time of the 1e-11 run.
- Feasibility's used the per-step cost and the locked step size.

**FLAG on my own R1. I adopt Feasibility's model with a range.** The ratio of dt to τ (the "lock factor") is 4.1 at K = 1e-9, not 2.8. Their model therefore over-predicts the 1e-9 run (1.9 h predicted, 1.33 h measured). I take a lock factor of 5τ as the lower bound and add 50% Newton overhead for the upper bound.

| K | Full 0.13 s window | **0.025 s restart window around peak (recommended)** |
|---|---|---|
| 3e-10 | 6.5 h (3.6–9.7) | **1.2 h (0.7–1.9)** |
| 1e-10 | 19.4 h (10.9–29.2) | **3.7 h (2.1–5.6)** |
| 1e-11 | 195 h | 37 h (infeasible) |

#### 3. Inventor (R2-C3)

**Frozen-angle error.** The Reynolds film has a pressure-driven flux q ∝ h³. An axis offset e moves the gap arc locally by dh = e·dθ·cos φ, so δQ/Q = (1 + dh/h)³ − 1 ≈ 3e·dθ/h. **PASS**: the table is reproduced (e = 10 µm: 0.26 / 0.52 / 1.05%; e = 25 µm: 0.65 / 1.31 / 2.62%; e = 50 µm: 1.31 / 2.62 / 5.24%, exact 5.33%). The 5% limits of 4.8°, 1.9° and 0.95° are reproduced. Two caveats:
- It is a local upper bound. cos φ changes sign along the arc, so the net flux error over a long arc is smaller.
- The WSS error, which is linear in h, is one third of the flux error.

**Back-coupling.** D = μR²A/h = 2.0e-9 N m s per hinge, against Inventor's 2.2e-9. The friction torque relative to 0.2 mN m is 0.9e-4 for the mean closing motion and 1.0e-3 at peak. With both ears it is 1.8e-4 to 2.0e-3. **PASS.** One caveat: 0.2 mN m is the systolic torque. The closing torque is not known, and if it is larger the ratio only gets smaller.

**Leakage.** With w = 6 mm, L = 1 mm and 80 mmHg, the leak is **0.38 mL/beat**. Clinical's w = 4 mm, L = 2 mm and 100 mmHg give 0.16 mL/beat. **Both PASS.** The 2.4x difference comes only from the unverified w/L ratio, so I would report 0.15–0.4 mL/beat until w and L are measured in the CAD. The gap sensitivity is −49% at 40 µm and +73% at 60 µm (PASS).

**Other Inventor numbers.** All PASS: 0.77 shear per degree, 11.5 at 15°, 576 µm of slide, h²/ν = 0.76 ms, and a 3.8% entrance-loss share.

#### 4. Clinical (R2-C4)

**Stress accumulation per pass.** At 100 mmHg and L = 2 mm: τ_w = 167 Pa, U = 0.40 m/s, transit time 5.04 ms, so τ_w·t = **0.84 Pa s. PASS.** About 4.2 passes reach 3.5 Pa s. Two things are new here:
- **τ_w·t = 6μL/h, which is independent of dp.** Linear stress accumulation per pass in a slot depends only on L/h. A narrower gap gives more stress accumulation per pass, which matches Clinical's washout argument.
- 0.84 Pa s is the conservative bound. The flux-weighted mean over streamlines is 3μL/h = **0.42 Pa s, or about 8 passes**. Report both, or use pathlines.

**Giersiepen power-law damage.** HI = 0.133% and platelet index = 0.383%. PASS, with the constants still needing a literature check.

**Carreau.** Clinical's percentages (0.7 / 2.2 / 8 / 20 / 37 / 57%) are reproduced, but they are measured against μ∞ = 3.45 mPa s. Against the model's 3.5 mPa s:
- In the gap the Newtonian model is −0.8%, a slight over-prediction.
- In the pocket at 250 / 100 / 50 s⁻¹ the error is +18 / +35 / +55%.

The conclusion does not change: Newtonian is fine in the gap, and a Carreau case is needed for residence time.

#### 5. Inventor N4: split each end gap at its mid-plane with a sliding interface (R2-C5)

**Shear.** **PASS, with a condition.** A rotation about z maps every plane z = const onto itself (numerically |dz| = 0). The end-gap mesh shear therefore drops from 11.5 to 0, and the in-plane slide is carried by the identity pair. The condition is that the rotating slab must be bounded laterally by a **cylinder about the pivot axis**, not by the ear outline. Otherwise its side faces sweep through the fixed slab. In other words, N4 is a thin rotor disc of radius at least 2.3 mm with planar and cylindrical sliding pairs, which is rotating-machinery topology.

**5 µm concentricity tolerance.** PASS as adequate. A cylindrical interface placed inside a lateral gap survives the whole sweep if e·θ + 2σ_r < h/2 − h/6 = 16.7 µm (here e is the pivot-axis offset and σ_r the surface scatter). σ_r = 5 µm then allows e up to about 25 µm. σ_r = 10 µm fails even at e = 0.

**Mesh cost.** **FLAG.** A 25 µm slab with 3 isotropic tet layers needs about 6 M elements per slab, and there are 8 slabs. Swept prisms (3 layers of 75 µm triangles) need about 20 k. **N4 is feasible only with swept or prism meshes.**

#### 6. Corrected validation statistics (R2-C6)

**Real data at 70°, K = 1e-9.** d = −8.64% and −8.74%. A and B are mirror images, not replicates, so this is effectively one sample and no t-test is possible. The point estimate lies outside ±5%, so **Phase A fails deterministically on the current run.** The +25.6% Fy error fails any force gate.

**Synthetic 4-angle cases** (A and B averaged per angle):

| Case | Mean ± sd | TOST p (margin ±5%) | 90% CI | Slope (95% CI) | Bland–Altman limits of agreement |
|---|---|---|---|---|---|
| Bias about −3% | −2.82 ± 0.79% | **0.006, equivalent** | [−3.76, −1.89] | 0.968 ± 0.043 | [−4.4, −1.3] ± 2.2 |
| Bias about −4.3% | −4.30 ± 0.75% | **0.080, not shown** | [−5.19, −3.41] | 0.953 ± 0.041 | [−5.8, −2.8] ± 2.1 |
| FyA, margin ±10% | +5.4 ± 2.1% | 0.011, equivalent | [2.9, 7.9] | 1.061 ± 0.081 | [1.2, 9.6] ± 5.9 |

All three cases have R² = 0.9998, which again shows that an R² gate tells you nothing.

**Benjamini–Hochberg FDR across the TOST p-values:** two outputs are equivalent and the −4.3% case is not. With n = 4, |bias| < 5% can be demonstrated only when the mean bias is at most about 3–3.5% and the sd is about 1% or less.

**Recommendation.** Fold in the numerical uncertainty using a V&V-20-style interval, E ± 2u_val with u_val built from the two GCIs (needs literature check).

#### 7. Other conflicts (R2-C7)

- **Halo thickness.** Halo = ΔV·t/(2V): 152 / 196 / 240 µm for t = 0.7 / 0.9 / 1.1 mm. Advocate's ~200 µm (t = 0.9 mm) and my 150–240 µm come from the same formula. Either way the halo is 3–5 gap widths, so the hinge is sealed. PASS.
- **Added inertia.** I_a = (9π/128)ρc⁴s gives 5.3x for Feasibility's chord and span, and 1.3–6x across plausible geometry. All three agents agree within the geometry spread. **Measure it directly:** in the body-fitted model, prescribe α with the flow at rest and take I_a = −T/α.
- **Challenger spot checks.** All PASS: Couette 1.3 Pa, gap Womersley 0.076, and GCI memory 2.9x under N^(4/3) scaling.

#### 8. Proposed numeric gates for the Judge (R2-C8)

| Gate | Metric | Threshold |
|---|---|---|
| G0 | LinRes on every step | < 1e-10 |
| G1 | Ghost volume vs. CAD 208.5 mm³ | < 1% |
| G2 | A/B torque symmetry (4 significant figures) | < 3% |
| G3 / G3b | TOST on torque / FyA, n = 4 angles | ±5% / ±10%, p < 0.05 |
| G4 | Slope 95% CI | inside [0.95, 1.05] |
| G5 | Bulk GCI_fine with the current mesh as finest; asymptotic ratio | < 5%; 0.9–1.1 |
| G5b | Gap GCI for flux / WSS (body-fitted submodel) | < 5% / < 10% |
| G6 | \|T(1e-10) − T(3e-10)\|/T | < 2%, otherwise report as K-limited |
| G7 | Ghost used inside the gap only if 12KT/h³ and slip excess are both | < 5% (i.e. never at feasible K) |
| G8 | Film validity: ε·Re and h²ω/ν | < 0.3 each |
| G9 | Film vs. 3D body-fitted: flux / peak WSS | 10% / 15% |
| G10 | Frozen-angle window 3e·dθ/h | < 5% (this sets dθ) |
| G11 | Hinge friction torque / blade torque | < 1e-2 |
| G12 / G13 | Mass balance / effect of stabilization | < 1% / < 2% |
| G14 | N4 geometry: e·θ + 2σ_r | < 16.7 µm |
| G15 | Carreau vs. Newtonian change in residence time / WSS | > 5% means Carreau goes in the main results |

#### Bottom line (ranking, updated)

1. **N2 (real-hinge box submodel) together with O4 (two-scale):** body-fitted, so the leak does not arise (R2-C1). Its frozen-angle error is bounded and computable (R2-C3), one-way coupling is justified to about 1e-3, and a local gap GCI is affordable.
2. **N1 (closed-angle leakage case):** the only route to the clinically dominant diastolic numbers. Stress accumulation is 0.42–0.84 Pa s per pass, and leakage is 0.15–0.4 mL/beat.
3. **O6 (pivot optimization):** it supplies e, which sets gates G10 and G14 and the leakage sensitivity (h³).
4. **O1 (hybrid ghost + film):** valid only for bulk flow and the film, with the gap taken from the film. Phase A fails until FixChi is run; the K study is feasible in 0.025 s windows.
5. **N4 (mid-plane split), as an O3 rescue:** geometrically exact, but it needs rotor-disc topology and swept meshes.
6. **N5 (stationary K screening) and N3 (2D closure model):** cheap supporting tools.
7. **O2 (ghost alone):** disqualified in the gap at every K that can be run.
8. **O3 as it stands:** kinematics only.
9. **O5 (overset):** too costly in degrees of freedom.

---

## Round 2 — Clinical

### R2 — Clinical / Physiology Expert: rebuttal

Labels follow OPTION_MAP_R2: my R1 "O7" is now **N1**, my "O8" is now **N3**.
Numbers come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/clinical/r2.py`, plus the R1 `calc.py`.
Literature status is unchanged from R1: I checked abstracts and search snippets only, because the full-text hosts were blocked.

#### 1. To the Advocate: is gap WSS from a Reynolds film clinically useful without pocket washout?

**Partly, and only as half of a trade-off.** The best experimental hinge-gap study I found, the SJM micro-PIV study with three gap widths (PubMed 24976188, abstract), reports that the narrowest gap had the **lowest** shear but the **poorest washout**. Its conclusion is that hinge design has to balance the two.

The lubrication scalings show why WSS alone points the wrong way. At fixed Δp:

- Leakage flux Q ∝ h³, wall shear τ_w ∝ h, transit time ∝ 1/h², so stress accumulation per pass (τ·t) ∝ 1/h.
- Going from 50 to 40 µm gives Q ×0.51, τ_w ×0.80 and τ·t per pass **×1.25**.
- Going from 50 to 60 µm gives Q ×1.73, τ_w ×1.20 and τ·t ×0.83.

If you report only film WSS, a narrower hinge looks safer. Exposure time and washout say it is worse. So film WSS has clinical value only when it comes with the film's own transit time (from which it gives an analytic τ·t per pass), and with a pocket residence or washout measure from a 3D body-fitted model. The film can supply the first, cheaply and with high accuracy, in the uniform-gap region. It cannot supply the second.

I accept the Advocate's validity numbers. ε·Re is 0.12 to 0.60 in diastolic leakage, and **in diastole the film does not need the ghost bulk at all**. With the valve closed, the driving Δp is the transvalvular pressure, which can be prescribed from the waveform. That makes a standalone diastolic film model clinically credible for gap WSS and flux, and it avoids the problem in §2.

#### 2. To the Challenger and Numerical analyst: ghost leaflets leak at closure

The clinical consequence is decisive. Numerical C2 gives a leak through the ghost leaflet of about 2 L/min at K = 1e-11, a K value that cannot be run. That is about **17 mL over a 0.5 s diastole**. At the K that does run (1e-9) the leak scales to about **1.7 L per diastole**.

Compare these with:

- A true hinge leakage of **0.13 to 0.48 mL/beat** for all 4 hinges (slot band: span 4–6 mm, L = 1–2 mm, Δp = 80–100 mmHg).
- A total in-vitro mechanical-valve regurgitant volume of **5 to 11 mL/beat**. This figure is from a search snippet of a multi-valve in-vitro comparison; the exact citation needs a literature check.

Even at the infeasible K, the false leak is **35 to 130 times the hinge leak** and **larger than the whole-valve regurgitation**. The Challenger's flux-ratio and Numerical's velocity-ratio numbers disagree in magnitude, but both say the same thing clinically: the leak is far above the signal at any runnable K.

So **no diastolic or closure output from a ghost bulk model (O1 bulk, O2) is clinically credible.** That includes leakage jets, regurgitant volume, washout of the pocket in diastole, and mouth pressures fed to a film. Penalized flow dominates exactly the region we care about.

The fix is to restrict ghost use to systole, the opening phase and torque, and to drive every diastolic and closure hinge output from body-fitted geometry (N1/N2) or from a prescribed Δp (film). The thesis text should say this explicitly. A reviewer will catch it otherwise.

#### 3. To the Inventor: N2 (frozen-angle body-fitted mesh with transpiration and a lumped law)

- **Leakage phase (about 60% of the cycle): credible.** At the closed angle nothing moves, so transpiration is not needed. N2 reduces to my N1 and is exact apart from mesh error. It delivers gap WSS, leakage Q, jets, pocket residence time and pathline stress accumulation. The prerequisite is the **closed-stop angle**, which is still UNVERIFIED. My estimate is about 30°, inferred from an 85° opening (search-verified) and a 55° traverse (snippet). It must come from the CAD or the advisor.
- **Closure phase: credible for hinge flow, with two clinical caveats.**
  1. Couette (wall-motion) shear is only 1 to 15% of pressure-driven gap shear (4–20 Pa against 133–333 Pa). So closure hinge hemodynamics are governed by the Δp(t) history across the closing valve. N2's credibility at closure therefore depends on where its box boundary conditions come from. If they come from the ghost model, §2 applies and it fails. They must come from the waveform, the v21 pin-joint model, or body-fitted runs.
  2. Closure sweeps 85° → about 30°, but every existing geometry and validation covers 70–85°. N2 needs about 55° of new frozen angles, and the ears must not intersect the pocket along the way. This is O6's job.
- **What N2 does not capture:** leaflet impact, rebound, water-hammer and cavitation. The thesis should state these as out of scope.
- **Clinical implication of negligible back-coupling.** The Inventor's friction-torque ratio of about 1e-3 and leakage of about 0.5% of stroke volume mean the hinge does not change leaflet motion or bulk flow. Hinge thrombogenicity is therefore a **post-processable function of (θ(t), Δp(t), h)**. Gap-width and hinge-shape design studies can then be run on the small box without repeating blade FSI, which is clinically the most useful property in the whole debate.
- **Caveat:** this holds for a clean hinge. The clinical failure mode is thrombus or pannus that *restricts* motion, and then the coupling becomes strong. Mention it as a limitation.
- **The lumped R/K/C law is engineering novelty, not a clinical output.** It is fine if cheap, but it is not part of the clinical minimum.

#### 4. To Feasibility: what survives on about 13 engineering days

| # | Clinical output | Importance (1–5) | Deliverability (1–5) | Days | Score |
|---|---|---|---|---|---|
| a | Closed-angle gap WSS, leakage Q, jets (N1 body-fitted) | 5 | 4 | 2–3 | **20** |
| b | Pathline stress accumulation (τ·t), % of particles > 3.5 Pa·s, closed vs open | 5 | 4 | 1 (post-processing of a, e) | **20** |
| c | Pocket residence time / washout, closed angle, Newtonian vs Carreau | 5 | 3 | 1 | 15 |
| d | Gap sensitivity 40/50/60 µm (parametric on N1) | 4 | 4 | 1 | 16 |
| e | Systolic hinge WSS and washout from the *existing* v20 fixed-angle runs | 3 | 5 | 0.5 | 15 |
| f | Film vs 3D gap WSS in the uniform region (O1 Phase C, diastole-driven) | 3 | 3 | 3–4 | 9 |
| g | Closure sequence (N2 windows, 85 → 30°) | 4 | 2 | 4+ | 8 |
| h | N3 2D closure cross-section | 2 | 4 | 1.5 | 8 |
| i | Cavitation or impact | 3 | 0 | — | 0, out of scope |

Items a to e (with A/B torque validation and GCI already in the plan) need about 6 to 7 days, so they fit. f fits if a passes by about Nov 18. g is a stretch.

#### 5. Minimum clinically defensible deliverable set, with acceptance criteria

1. **N1 closed-angle leakage case (a, c).**
   - Mesh GCI < 5% on Q and on the area-averaged gap WSS, and < 10% on the 99th-percentile WSS (the maximum itself is mesh-noise).
   - Mass imbalance < 1%.
   - Gap WSS in the uniform region within ±30% of the analytic slot value (h/2)Δp/L, which is 133–333 Pa for L = 1–2 mm (computed).
   - Peak WSS of the same order (×0.3 to ×3) as the 612 Pa diastolic peak in the Simon et al. 2010 abstract.
   - Hinge leakage of 0.05 to 1 mL/beat for all 4 hinges (analytic band 0.13–0.48), and below 10% of the literature total of 5–11 mL/beat (the total needs a literature check).
2. **Washout.**
   - Report the mean fluid age in the pocket against the diastolic duration (about 0.5 s).
   - If Carreau changes the residence time by more than 10%, Carreau becomes the reported model.
   - There is no literature threshold for washout (needs literature check), so report it comparatively across 40/50/60 µm.
3. **Stress accumulation (b).**
   - Linear τ·t with the 3.5 Pa·s (35 dyn·s/cm²) threshold, plus one power law (Giersiepen) used for ranking only.
   - The fraction of particles above threshold should be within an order of magnitude of Dumont et al. 2007 (0.81% in regurgitation, 0% in forward flow, abstract checked).
   - **The forward-flow fraction must be lower than the regurgitation fraction.** This is the qualitative clinical check.
4. **Gap sensitivity (d).** Must reproduce the trend in the micro-PIV gap-width study: narrower gap gives lower τ_w and poorer washout.
5. **Stated limitations.** Impact, cavitation, thrombus feedback, the continuum assumption at h/d_RBC of about 6, and the UNVERIFIED closed angle and waveform source.

#### Bottom line (updated clinical ranking)

1. **N1 + O4**: the only path to credible dominant-phase (diastolic) hinge outputs. Every item in the minimum set comes from it.
2. **O6**: needed to fix h at the closed angle and to make the closing sweep admissible. Leakage ∝ h³.
3. **N2**: the right path for closure-phase hinge flow, but only with Δp boundary conditions that do not come from the ghost model. It is a stretch.
4. **O1, film part only (diastole-driven standalone)**: high-accuracy gap WSS and transit time. It is clinically incomplete without washout.
5. **N3**: a cheap qualitative closure story.
6. **O2 / O1 ghost bulk**: systolic torque only. Clinically disqualified for diastole and closure because the porous leak is 35–130 times the hinge leak even at an unrunnable K.
7. **O3, N4**: kinematics only. **O5**: infeasible.

---

## Round 2 — Feasibility

### R2: Feasibility Engineer rebuttal (new options N1–N5, revised O1, revised plan)

Scripts are in `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/feasibility/`: `r2.py` (sizes and K runtime) and `r2plan.py` (plan hours).

Cost model:
- 45.7 s per time step at 1.58 M DOF. This is measured, from the K = 1e-11 run.
- Time per step scales as N^1.7. The exponent is ASSUMED; anything from 1.5 to 2 moves the box-model estimates by at most a factor of 2.
- Memory: 15–28.7 GB × (N/1.58 M)^1.3–1.5. This is the Numerical analyst's bracket.

#### 1. K-runtime conflict (my 6.5 / 19 h vs Numerical's 2.8 / 5.7 h)
Both numbers are extrapolations from one stalled run, and they differ only in the assumed time-step law.

| K | My estimate: dt ≤ 2.8τ | Numerical: power law K^-0.63 | 0.025 s restart window (my basis) |
|---|---|---|---|
| 3e-10 | 6.5 h | 2.8 h | 1.2 h |
| 1e-10 | 19.5 h | 5.7 h | 3.7 h |

- **My estimate (upper bound).** I take the step at 2.8τ, which was observed when the step froze. τ = ρK/μ is the penalty time scale. I assume the Study 2 window runs 0.08 → 0.21 s, which gives 788 steps in 10 h and 45.7 s/step.
- **Numerical's estimate.** They assume the K = 1e-11 run started at t = 0 and kept a constant sim-s/h rate. That implies 3.5 s/step if the step had been 8.5 µs throughout, which is not credible for MUMPS at 1.58 M DOF. So the step froze only late in that run, and the run's rate was not constant.
- **The truth lies between the two.** The K = 1e-9 solver log settles it in 5 minutes: count the steps and read the time per step. Record this as task T0.
- **It does not change the plan.** Numerical's own finding C4 is that the K study plateaus on the 0.61 mm mesh. Only restart windows (≤ 3.7 h) and N5 stationary screening (about 1.5 h for 4 K values) are worth running. All of them fit in one night under either estimate.

#### 2. O1 re-assessed: what would local refinement cost?
For a K study to show genuine convergence, the mesh needs h ≲ sqrt(K) = 10–32 µm around the leaflets. Leaflet wetted area is about 1,190 mm² (both leaflets, t = 0.7 mm). Extra isotropic tets needed:

| Band thickness | Element size h | Extra tets | Approx. extra DOF |
|---|---|---|---|
| 200 µm | 100 µm | 2.0 M | +1.6 M (doubles the model) |
| 200 µm | 60 µm | 9.4 M | +7 M |
| 100 µm | 30 µm | 37 M | +29 M |
| 100 µm | 10 µm | 1,000 M | — |

- For a **moving** leaflet the band must cover the swept volume, about 565 mm³. That is 4.8 M tets at 100 µm or 22 M tets at 60 µm.
- **Conclusion: genuine K-convergence of the moving ghost leaflet is infeasible on any hardware we have.**
- **What is feasible:** at fixed angle the leaflets are already separate geometry domains. A conforming boundary-layer mesh on the leaflet surface (5 prism layers, 0.3 mm tangential size) adds only about 0.15 M elements. That makes a clean, separate K-sensitivity study at 70° possible, 0.5 engineer-day.
- **O1 therefore shrinks to:**
  - the ghost leaflet as a blade-scale cross-check (Phase A plus K-sensitivity);
  - the Reynolds film validated against a body-fitted 3D hinge (N2), which is the novelty;
  - Phases D/E (moving χ with the film, two-way FSI) dropped from the thesis.

#### 3. New options: COMSOL steps, engineer-days, compute

##### N2: frozen-angle body-fitted hinge box with moving-wall velocity and a fitted lumped law (core)
**COMSOL steps:**
1. Copy v20. Intersect the blood domain with a 6×6×4 mm Block around one hinge pocket. The hinge region is mirror-symmetric (A/B and ±z), so one box per hinge type is enough; that symmetry is UNVERIFIED and should be checked with O6.
2. Make the leaflet stump a Rotate node with parameter `th_k`. Do a **Parametric Sweep over th_k**: Sweep rebuilds geometry and mesh per value. Use "Free tetrahedral + Boundary Layers" with 5 layers across the 50 µm gap. This is the robust choice; a swept mesh fails at odd angles.
3. Laminar Flow:
   - leaflet walls: Moving Wall `u = -thd*(y-yp)`, `v = thd*(x-xp)`; this is Inventor's transpiration argument;
   - box faces: Open Boundary or Pressure condition with p from an Interpolation function, exported from the global model at the 4 existing angles and interpolated in θ.
4. Lumped law:
   - Mouth cut-planes with `intop` flux Q and mean Δp.
   - Two unit Stokes solves per angle for C(θ): set thd = 1 with Δp = 0, then Δp = 1 with thd = 0.
   - Fit R(θ), K(θ), C(θ) with scipy least squares.
5. GCI with 3 box meshes (0.18 / 0.4 / 0.88 M DOF) at 70° and 84°.

**Memory:** 0.18 M DOF is my optimistic count and about 0.4 M is realistic, which needs 2–8 GB. The fine 0.88 M mesh needs 6–14 GB. **The real-gap GCI fits on the desktop.** That answers the Challenger's §2.7 point that the reference has no GCI.

**How many angles:** the moving-wall error is 3eΔθ/h, where e is the eccentricity of the ear about the pivot. If O6 gives e ≤ 10 µm, 1° spacing is enough: 16 angles over 70–85°. If e is about 50 µm, use 0.5° (31 angles). The moving-wall windows, 30 × 0.5°, are each only 2.7 ms long at dt = 1e-4.

**Cost:** 5–6 engineer-days, about 12 h compute, all of it overnight work. **P(success) about 75% (judgement).**

##### N1: closed-angle body-fitted leakage (core, clinically dominant)
**COMSOL steps:**
1. Rotate the v20 leaflets to the closed angle. This angle is unknown, so **ask Dr. Bellofiore today**.
2. Check the rim gap at the housing seat. If the CAD shows contact, a seat gap must be set. That is a model value, so it needs advisor approval.
3. Cut along the x = 0 and z = 0 symmetry planes to a quarter model, about 0.6 M DOF. Steady flow at Re_h ≈ 10 makes symmetry plausible (UNVERIFIED).
4. Stationary solve with Auxiliary Sweep continuation on Δp from 0 to 100 mmHg.
5. Age-of-fluid: add a Transport of Diluted Species with source = 1 to get pocket residence time.

**Cost:** 2.5–3 engineer-days, ≤ 3 h compute. **P about 70%.** The main risk is geometry: seat contact and rotating the leaflet by a large angle.

##### N3: 2D hinge cross-section with prescribed closing motion (optional)
**COMSOL steps:** a 2D xy section through the ear. If the ear outline is concentric with the pivot, use a Rotating Domain with a circular identity pair (Flow Continuity) in the mid-gap. Otherwise use 2D Moving Mesh with remeshing, which is cheap in 2D. Put at least 10 cells across the gap and use dt of about 1e-5 over 30 ms.

**Cost:** 2–3 engineer-days, about 2 h compute. **P about 80%.**

**Limit:** it misses the end-face and axial leakage paths. It is a good closure-squeeze appendix but **not in the critical path**.

##### N4: mid-plane sliding split in the ALE model (stretch, drop for thesis)
**COMSOL steps:**
1. Rebuild v20 with **Form Assembly** so each end gap splits into two 25 µm slabs: a housing half (fixed) and a leaflet half (rigid rotation).
2. On the mid-plane Identity Pair, apply Laminar Flow **Flow Continuity**. Make the destination side finer.
3. In Moving Mesh, the leaflet-half slab gets prescribed rigid rotation, and the rest keeps the current Yeoh/sliding settings.
4. The housing slab must cover the whole swept footprint.

**Cost:** 6–8 engineer-days just to rebuild the 118-boundary geometry and re-attach every selection. Flow-on ALE for one opening is then 20–60 h (UNVERIFIED). It is also UNVERIFIED whether automatic remeshing keeps non-conforming pairs in an assembly. **P about 30%. It does not fit in the calendar.** Leave it for a paper revision or PhD work.

##### N5: stationary K screening
Already costed above: about 1.5 h, 0.5 engineer-day. Keep it.

#### 4. Responses to schedule claims
- **Challenger: "only A, possibly B; C–E will not pass".**
  - Agreed for D and E.
  - Disagree for C *if it is redefined*: the film is checked against the N2 body-fitted box, not against the ghost model. The film has under 50k DOF, the one-way coupling is justified (Inventor: hinge torque about 1e-3 of blade torque; Advocate: about 0.2%), and it costs about 2 engineer-days.
  - Their "K = 1e-10 about 13 h" sits inside my 5.7–19.5 h bracket, but it is moot because the K study plateaus (§1).
  - Their "GCI beyond 32 GB" is correct for global refinement. It is solved by coarsening (the current mesh becomes the fine grid) plus GCI inside the box.
- **Advocate: two-way mouth exchange plus pocket blocking.** The flux return to 3D needs geometry cuts plus fixed-point iteration, about 3–4 engineer-days. The physics shows the feedback is negligible, so cut it and keep one-way coupling.
- **Inventor: "30 windows at 0.5°".** Compute is trivial (about 5 h) because each window is short. The cost is per-angle remeshing. Let O6's e decide between 16 and 31 angles.

#### 5. Revised plan: O4 core = N2 + N1 + O6, with O1-C (film vs box) as the novelty and ghost A + N5 as the cross-check

Compute (`r2plan.py`): **27.5 h total, or 55 h with a 2× rerun factor.** About 17 nights × 12 h = 204 h are available on one desktop, so **HPC is not needed.** Engineer time is the limit: about 13 days are needed and about 12.5 are available. N3 is therefore dropped unless there is slack.

| Day | Day work (engineer) | Unattended overnight run |
|---|---|---|
| Tue Nov 10 | T0: read the K = 1e-9 log for steps and s/step. Save the ALE log. FixChi: VA/VB = 208.5. Email the advisor about the closed angle, seat gap and clearance | A: ghost 70/75/80/84 (5.3 h) |
| Wed Nov 11 | O6 pivot fit: eccentricity e and sliding. Phase A analysis (equivalence CI and slope, 4 significant figures) | N5 stationary K screening + K = 1e-10 restart (5.2 h) |
| Thu Nov 12 | Coarsened meshes (ghost and v20). Begin the N2 box geometry | Ghost and v20 coarse GCI runs (1.1 h + reruns) |
| Fri Nov 13 | N2: box, Rotate(th_k), boundary-layer mesh, BC interpolation from global | Box mesh test sweep, 16 angles, mesh only |
| Sat–Sun Nov 14–15 | (light) | N2: 3-mesh GCI at 70/84° (1.8 h) + 16-angle frozen sweep (5.3 h) |
| Mon Nov 16 | Check GCI; mouth flux/Δp extraction | N2 unit Stokes solves for C(θ) (0.2 h) |
| Tue Nov 17 | Moving-wall windows set up (Δθ from e) | N2 moving-wall windows (5 h) |
| Wed Nov 18 | Fit lumped law R, K, C(θ); leave-one-angle-out check | Held-out 77.5° box run |
| Thu Nov 19 | N1: rotate to closed angle, seat gap, quarter model, mesh | N1 Δp continuation + age-of-fluid (≤ 3 h) |
| Fri Nov 20 | N1 post-processing: leakage, gap WSS, washout | N1 2-mesh check |
| Mon Nov 23 | O1-C: Thin-Film Flow on pocket walls (Border p from box), or a Coefficient Form Boundary PDE fallback | Film at 16 angles (0.8 h) |
| Tue Nov 24 | Film vs box: flux and WSS within 10% gate; validity map of eps·Re | Reruns / buffer |
| Wed Nov 25 | Figures; optional N3 start if ahead of schedule | N3 (1.7 h) if started |
| Thu–Fri Nov 26–27 | Thanksgiving / buffer. **Code freeze Nov 27** | Any reruns |
| Sat Nov 28 – Thu Dec 3 | Thesis writing and slides | none |

**Overnight-safe tasks** (scriptable through Java/Study sequences with no interaction): all Phase A runs, N5, the K restart, coarse-mesh GCI runs, the N2 parametric angle sweep and box GCI, the unit Stokes solves, the moving-wall window chain (Batch/Sweep), N1 continuation, the film sweep, and N3. Always check LinRes ≈ 1e-14 the next morning.

**Fallback at each gate:**
- If the box BC interpolation fails, use pressure-only boundary conditions from the 4 existing angles.
- If N1 geometry fails by Nov 20, report N2 at 70° with diastolic Δp applied as a "closed-like" sensitivity, and state it as a limitation.

#### Bottom line (feasibility ranking, R2)
1. **N2 (within O4):** body-fitted 50 µm gap with GCI, fits in 2–8 GB and about 12 h of compute.
2. **N1:** the clinically dominant phase, about 3 days of work; it depends on the advisor giving the closed angle today.
3. **O6:** 1 day, and its eccentricity value sets the N2 angle spacing and the validity of the moving-wall method.
4. **O1 limited to Phase C (film vs N2 box):** novelty for about 2 days of work. Phases D/E are dropped.
5. **O2 + N5:** a cheap blade-scale cross-check, reported as K-sensitivity only.
6. **N3:** cheap but off the critical path.
7. **N4:** 6–8 days of rebuild at about 30% odds, out of reach.
8. **O3:** a kinematics appendix only.
9. **O5:** infeasible.

---

## Round 2 — Inventor

### R2 — Inventor rebuttal: defending N2 (box submodel + lumped law) and N4 (end-gap sliding split)

Numbers come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/inventor/r2.py`, plus `calc.py` from R1.
Geometry not taken from the handoff (leaflet stump size, gap length L, window run time) is marked **UNVERIFIED**.

#### 1. Where the openings agree, and what that does to the debate
Challenger, Numerical, Clinical and Feasibility all rank O4 first or second. Challenger's O7 is the same construction as my N2. Numerical C9 independently proposes "end-face slip with a gap-local mesh that rotates rigidly", which is N4.

So the open questions are not "two-scale or not". They are:
- (i) Is the moving-angle part of N2 accurate and mass-consistent?
- (ii) Is it novel beyond Simon et al.?
- (iii) How does it relate to O1's film?

#### 2. Objections to N2, answered

**(a) Is the full rigid-body wall velocity on a frozen mesh mass-conservative?**
- For a closed rigid surface, ∮(ω×r)·n dA = 0 exactly, by the divergence theorem with div(ω×r) = 0.
- **But in the box submodel the leaflet is an open stump cut by the box faces, so that integral is not zero.**
- I estimated the swept flux of the cut cross-section as Ω·r_cut·t·l_cut. With t = 0.7 mm, l_cut = 4 mm and r_cut = 4 mm (all UNVERIFIED), it is **0.10 mL/s at mean closing speed (8.7 rad/s) and 1.1 mL/s at Ω = 100 rad/s**. The hinge leakage itself is about 0.19 mL/s.
- An inconsistent box boundary condition would therefore swamp the hinge flux. This was a real gap in my R1.
- Fix:
  - Impose a velocity condition from the global model on at most all box faces but one, and an open (pressure) condition on the downstream face, so the displaced volume has somewhere to go.
  - **Gate:** the discrete mass balance |Σ faces ∫u·n + ∫_wall u_w·n| < 0.5% of Q_hinge. COMSOL's P1P1 GLS formulation is not exactly locally conservative, so this must be measured, not assumed.

**(b) Non-concentric arcs, and R = 1.92 mm (Numerical) vs my 2.2 mm.**
- The refined window error, with the mesh frozen at the window midpoint (largest deviation δθ/2), is

  ε_win ≈ 3·e·(δθ/2)/h  +  R·(δθ/2)/L

  - First term: the gap-thickness change from eccentricity e.
  - Second term: translation of a housing-defined mouth. This is conservative; it vanishes for a uniform concentric gap whose mouths are set by the ear itself.

  | e | δθ = 0.25° | 0.5° | 1° |
  |---|---|---|---|
  | 0 um | 0.42% | 0.84% | 1.68% |
  | 50 um | 1.07% | 2.15% | 4.29% |
  | 200 um | 3.04% | 6.07% | 12.2% |

  These values use R = 1.92 mm and L = 1 mm. With R = 2.2 mm, the values change by at most 0.25 percentage points. R also enters the Couette/pressure flux ratio, which moves only from 1.32% to 1.51% (at 80 mmHg), so **the R conflict does not matter for N2**.
- To meet a 3% budget, the number of windows per 15° swing is **11, 22 or 61 for e = 10, 50 or 200 um**. At 0.3–1 h per window (UNVERIFIED), 61 windows cost 18–60 h. That is affordable even in the worst case.
- e is unknown until O6 runs, so **O6 must come before N2's moving windows**.
- Physics in our favour: during the closing swing, Δp is only a few mmHg (Numerical C8), so the window error applies to a small flux. The 80–100 mmHg leakage phase happens at a fixed angle, where frozen geometry has no geometric error at all.

**(c) Closure impact, h²ω/ν ≈ 3.8 (Numerical C8).**
- This kills a quasi-steady film, which drops ∂u/∂t.
- **It does not kill N2.** Each window is a transient Navier–Stokes solve, so the unsteady term is kept. Only the geometry is frozen.
- Over a 0.2 ms impact the angle sweep is 0.37°, 1.1° or 3.4° at Ω = 32, 100 or 300 rad/s. One to fourteen windows at δθ = 0.25° cover it.
- The real limits at impact lie outside the method:
  - θ(t) through the stop comes from the global model's contact law;
  - water hammer and cavitation are compressible effects.
- **Impact therefore stays out of scope in the claims**, as Clinical recommended. N2 reports the closing swing up to the stop, plus the closed state.

**(d) Novelty vs Simon et al. 2010.**
- Simon, Ge, Sotiropoulos & Yoganathan (Ann Biomed Eng 2010) already drove a 3D SJM hinge simulation with measured leaflet motion using an immersed-boundary method. "Two-scale hinge CFD" is **not** new; I conceded this in R1.
- The defensible contributions are:
  - (1) a body-fitted, no-slip, local-GCI result at the **real 50 um** gap. Challenger reports Simon's gap as about 150 um, which I have not verified (needs literature check). An immersed boundary has its own interface layer, so a body-fitted reference is a step forward.
  - (2) an **a-priori window error bound**, together with its test (Section 5).
  - (3) a reusable lumped law, Δp_h = R(θ)(Q − C θ̇) + K(θ)ρQ|Q|/(2A²).
  - (4) a quantitative proof that one-way coupling is valid (friction torque about 1e-4 of blade torque).
  - (5) Section 3 below.
- I found no published lumped hinge law (needs literature check). (1) alone is incremental; (2)–(5) together make a methods paper.

#### 3. How N2 and O1 relate: a model hierarchy, and the combined novelty
Integrate the Reynolds film across the gap width for a uniform gap:

  R_film = 12μL/(w h³),   C_film = (couette term from U_t h/2)

These are the linear terms of N2's law. **N2's 0D law is the film PDE projected to 0D. The film is N2's law resolved in space. The 3D box is the truth for both.** What the film lacks is exactly what K(θ) captures in 0D: losses at the mouths and corners. Finite-width side-wall effects are only about 0.5–0.8% (0.63·h/w).

The combined contribution is a **three-level hierarchy, 3D body-fitted → 2D film → 0D law, with the error measured at each level and at each angle**. N2 is O1's missing validation reference; Challenger and Numerical both demand one. Two outcomes are possible:
- If the film matches 3D within 10% in the uniform region, O1's film can carry diastolic leakage.
- If it does not, the fitted mouth loss K(θ) can be added to the film's boundary conditions as a 3D-calibrated mouth correction, which is still novel.

Either outcome is publishable. **The ghost-leaflet half of O1 is not needed for this.** Numerical C2/C3 and Challenger 2.2 show the ghost cannot hold the closure pressure: 1.2–3.9× gap flux and m/s-scale leakage at K = 1e-9.

#### 4. Clinical's N1 and N3
- **N1 (closed-angle diastolic leakage)** is N2's first and most valuable case, at θ = θ_closed and θ̇ = 0. The geometry is fixed, so the window error is zero and only the mesh GCI applies. I fully endorse making it **deliverable #1**.
  - Prerequisites: the closed-stop angle from the advisor or the CAD (Clinical: UNVERIFIED, about 25–30°), and a Carreau sensitivity for residence time in the pocket.
- **N3 (2D hinge section with prescribed closing)**: I propose to repurpose it as **the falsification test of N2's window error**.
  - In 2D, true ALE with remeshing is cheap and does not suffer the end-face problem.
  - Run the same kinematics in 2D with full ALE and with frozen windows. The difference is the frozen-geometry error measured directly, with no other modelling difference.
  - N3 then also delivers Clinical's closure-squeeze shear story.
  - Its limit: no end gaps (z), so it says nothing about N4.

#### 5. Falsifiable validation protocol for N2 (numeric thresholds)
| # | Test | Pass threshold | Falsified if |
|---|---|---|---|
| V0 | Discrete mass balance in the box (Section 2a) | < 0.5% of Q_hinge at every saved step | ≥ 1% |
| V1 | Box-size independence (two box sizes) and agreement with the full v20 model at 70° and 84° | Q_h, Δp_h, gap WSS (area-mean and 99th percentile) within 3% | > 5% |
| V2 | Local 3-mesh GCI in the gap (4/6/9 elements across h; bulk mesh fixed) | GCI < 5% on Q_h and mean WSS; < 10% on 99th-percentile WSS; asymptotic ratio 0.95–1.05 | otherwise |
| V3 | 2D ALE vs 2D frozen windows (N3), δθ = 1, 0.5, 0.25° | error ≤ the ε_win bound, and log-log slope vs δθ in 0.8–1.2 | error > 2× the bound, or slope outside 0.6–1.4 (error model wrong) |
| V4 | 3D δθ self-convergence over the closing swing | Richardson change < 2% between 0.5° and 0.25° | > 4% |
| V5 | Lumped law, leave-one-angle-out (70/75/80/84/closed) | Q prediction < 5%, max WSS < 10%; fitted inertance changes Q by < 2% | inertance > 2% (quasi-steady assumption falsified) |
| V6 | One-way coupling: hinge torque / blade torque; hinge leakage / forward volume | < 1% | > 2% (two-way coupling needed) |

Every threshold is fixed before the runs are done.

#### 6. N4, defended and downgraded
- In-plane rotation of a plane normal to the pivot axis is pure sliding, so the end-gap shear drops from 0.67–0.77 to 0 per degree.
- Numerical C9 says the inversion comes from **sliding nodes meeting pinned nodes at corners**. N4 moves that corner to the place where the pair face meets the pocket's lateral wall.
  - If the lateral clearance is also about 50 um and not concentric, the same in-plane shear reappears in 2D there.
  - N4 therefore works only if the lateral thin arcs are concentric to about 5 um, so that they can become cylindrical sliding pairs too, or if the lateral clearance is much larger than 50 um.
- Flow Continuity on a non-conforming pair is enforced weakly, so the flux mismatch across the pair must be under 1% of the gap Q.
- **I lower P(N4) from 40% to about 30%.** It is a stretch goal after O6 has measured the concentricity.

#### Updated ranking
1. **N2 (with N1 as its first case and O4 as the global driver)**: body-fitted at 50 um, local GCI is affordable, the error bounds are falsifiable, and the hierarchy gives the novelty.
2. **O6**: supplies e and h(θ), which set N2's window count and the clinically dominant leakage (∝ h³). Run it first.
3. **O1, film part only, validated against N2**: the combined hierarchy novelty. Drop the ghost-leaflet bulk for closure.
4. **N3 (repurposed as V3)**: cheap, directly falsifies the frozen-geometry error, and adds closure shear.
5. **O2**: the blade-scale bulk flow only, with a K-sensitivity check rather than a convergence claim.
6. **N4**: high payoff only if the arcs are concentric; ~30%.
7. **N5**: an accelerator for O2's K screening.
8. **O3 as it stands**: kinematics appendix only.
9. **O5**: infeasible in the time available.
