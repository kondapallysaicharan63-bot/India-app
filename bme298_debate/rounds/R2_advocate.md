# R2 — Advocate rebuttal: O1 revised to "O1-R", a validated Reynolds-film hinge law

Scratch scripts: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/advocate/r2.py`. My R1 numbers come from `calc.py` in the same folder. Labels follow OPTION_MAP_R2.

Summary of this round: most of the Round-1 attacks are correct, and they are aimed at the ghost-leaflet half of O1. The Reynolds-film half survives. It also turns out to be the physics-based version of Inventor's N2 hinge law. So the case I now make is for O1-R, which keeps the film and gives up the ghost-film coupling.

## 1. Challenger: "29.7x Poiseuille flux at K = 1e-9; ~1300 mL/s through the closed leaflet"

**Conceded, with a reconciliation.** I re-solved the 1-D composite slot: a 50 um gap next to a 300 um Brinkman layer. The gap carries 3.63x / 1.70x / 1.20x / 1.06x the no-slip flux and the porous layer adds 26.1x / 2.97x / 0.31x / 0.03x, for K = 1e-9 / 1e-10 / 1e-11 / 1e-12. That gives 29.7x total at K = 1e-9, exactly Challenger's number.

Numerical's "3.9x" is the gap-only part of the same physics (3.63x in my configuration). The two figures are consistent, not in conflict.

Darcy transleaflet velocity at 80 mmHg for a 0.7 mm leaflet: 4.35 m/s (K=1e-9), 0.44 m/s (1e-10), 4.4 cm/s (1e-11). This reproduces Numerical C2.

**Consequences for O1:**
1. The ghost leaflet must never be used to block the pocket, and never at the closed angle. My R1 §4 step 4 offered "penalize as solid" as one option for blocking. **That option is withdrawn.**
2. Wherever the film is active, the pocket is removed geometrically, and the ear-to-housing gap exists only in the film.
3. The leaflet is impermeable in the film by construction: the Reynolds BC on the ear is u·n = wall velocity.
4. Closure and diastole are taken **out of the ghost model altogether** (see §3).

## 2. Numerical

**"The R² > 0.95 gate always passes."** Conceded. I replace it with an equivalence test (TOST): the 90% CI of the bias must lie inside ±5%. With n = 4 effective pairs the CI half-width is 1.18% at sd = 1% and 2.35% at sd = 2%. The test therefore passes only if |bias| < 3.8% (sd 1%) or < 2.6% (sd 2%). I also add the slope CI. Torques are to be reported to 4 significant figures.

**"The K study plateaus, because mesh 0.61 mm >> sqrt(K)."** Conceded. I checked local refinement. A shell around both leaflet faces (about 930 mm², UNVERIFIED area) costs about 31 M tets at 30 um elements and about 557 M at 10 um. Neither fits on the desktop. So the ghost result becomes a **K sensitivity at fixed mesh** with p fixed at 1, reported honestly as such. It is no longer a convergence claim.

**"The GCI does not fit in memory."** Conceded. For the bulk model, the current mesh is the finest grid and two coarsened grids are added. The real-gap GCI moves to the small N2 hinge box, where it is cheap (Numerical C11, Inventor).

## 3. Clinical: "70–85° excludes the high-risk diastole and closure; O1 gives no pocket washout"

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

## 4. Feasibility: "14–18 engineer-days against about 13 available"

**Conceded.**
- Phase E (two-way ghost FSI) is **dropped**. Its missing dχ/dθ Jacobian and added-mass ratio of 1.3–6 do not fit in this timeline.
- Phase D (ghost + film, coupled during motion) is also dropped.
- What remains for O1-R is Phase C: about 3–4 days of film development (Feasibility's estimate), run on geometry and boundary conditions that N1/N2 produce anyway.

## 5. Inventor: "N2, and the hinge barely back-couples, so one-way O4 is justified"

**Agreed.** My R1 gave 0.2% for the hinge torque relative to the blade, and Inventor's number is 1e-4 to 1.1e-3. One-way coupling is justified. This also removes the main reason for a two-way ghost-film coupling.

N2's lumped law is dp_h = R(θ)(Q − Cθ̇) + Kρ Q|Q|/(2A²). The Reynolds film supplies its coefficients from first principles:
- R(θ) = 12μ ∫ ds / (w h(s,θ)³). The uniform case is 5.6e10 Pa·s/m³, i.e. Q = 0.19 mL/s per hinge at 80 mmHg (L = 1 mm and w = 6 mm UNVERIFIED).
- C(θ) comes from the Couette and squeeze terms.
- K is the only empirical term, and Inventor estimates the entrance loss at 0.2–3.8% of dp.

So **N2 fits R from 3D runs, while O1-R predicts R from the CAD clearance map without a 3D solve and validates it against N2.** Agreement within 10% would turn a fitted law into a predictive one that can be reused at any gap width, eccentricity or θ (including the 40/50/60 um sensitivity) without remeshing a 50 um gap. To my knowledge that predictive capability is the novelty (needs literature check; no search this round or last found a Reynolds-film hinge law for a mechanical heart valve).

## 6. Challenger's other points, briefly

- **Mouth pressures taken from the least accurate ghost region:** resolved. They now come from body-fitted N1/N2 runs.
- **Sharp χ validated vs. smoothed χ used:** moot. The ghost leaflet is now used only at fixed angles with sharp χ, as a blade-scale cross-check, with no moving χ.
- **Opposite-signed Tz/Fy errors:** still open. It is now a Phase A item. FixChi has to be run, followed by a traction-based check on a surface outside χ. If |ΔFy| stays above 10% after FixChi, the ghost model is reported as unsuitable for forces. That would not affect O1-R.

## 7. Revised O1 ("O1-R"): scope, phases, gates

**Scope.** Real-gap hinge flux and WSS from a Reynolds film on the CAD clearance map. It is validated against body-fitted 3D (N1 at the closed angle under diastolic dp; N2 frozen-angle boxes in systole) and then used for gap-width and eccentricity sensitivity. The ghost leaflet is optional and limited to a systolic, fixed-angle blade cross-check.

| Phase | Work | Gate |
|---|---|---|
| A′ | FixChi and ghost fixed angles 70/75/80/84 (optional cross-check) | TOST 90% CI inside ±5%; Fy within 10% |
| C1 | Film on CAD h(s); mouth dp from N1 (closed angle, 80–100 mmHg) and from N2 at 2 systolic angles | Film Q and area-mean gap WSS within 10% of body-fitted 3D in the uniform-gap region; eps·Re map reported |
| C2 | Film-predicted R(θ), C(θ) vs. N2-fitted values | Within 10–15%; otherwise the film is reported as a bounded approximation |
| C3 | Gap-width (40/50/60 um) and eccentricity (from O6) sensitivity using the film only | None (an output) |

Dropped: D, E, the K-convergence claim, and the ghost-film coupling. Fallback: if C1 fails, N1/N2 stand alone and the film is reported as a negative result with its validity map. That is still publishable as a methods note.

## 8. Primary deliverable or add-on?

**Add-on.** The primary deliverable should be **O4 + N2 + N1**: body-fitted, monolithic, local GCI, and covering diastole. **O1-R is the novelty layer on top**: the reduced, predictive hinge law. I now consider it the strongest novelty claim the thesis can make, but the thesis must not depend on it. This is a change from my R1 position.

## Bottom line (updated ranking)

1. **O4 + N2 (+ N1)**: an accurate, body-fitted real-gap result that covers diastole, with a GCI that fits on the desktop.
2. **O1-R**: a predictive Reynolds hinge law validated against #1, and the main novelty; cheap because it reuses N1/N2 boundary conditions.
3. **O6**: about 1 day, and it supplies the eccentricity e needed by O1-R, N2 and N4.
4. **O2 (fixed-angle, systole only)**: a blade-scale cross-check; it cannot be used in the gap or at closure (29.7x leak, 4.4 m/s Darcy).
5. **N5**: a cheap K-sensitivity screen, useful only if O2 is kept.
6. **N4**: a high-payoff rescue for O3 (about 40% chance of success); a stretch goal.
7. **O3**: kinematics only. **O5**: future work.
