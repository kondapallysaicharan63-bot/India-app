# R2: Feasibility Engineer rebuttal (new options N1–N5, revised O1, revised plan)

Scripts are in `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/feasibility/`: `r2.py` (sizes and K runtime) and `r2plan.py` (plan hours).

Cost model:
- 45.7 s per time step at 1.58 M DOF. This is measured, from the K = 1e-11 run.
- Time per step scales as N^1.7. The exponent is ASSUMED; anything from 1.5 to 2 moves the box-model estimates by at most a factor of 2.
- Memory: 15–28.7 GB × (N/1.58 M)^1.3–1.5. This is the Numerical analyst's bracket.

## 1. K-runtime conflict (my 6.5 / 19 h vs Numerical's 2.8 / 5.7 h)
Both numbers are extrapolations from one stalled run, and they differ only in the assumed time-step law.

| K | My estimate: dt ≤ 2.8τ | Numerical: power law K^-0.63 | 0.025 s restart window (my basis) |
|---|---|---|---|
| 3e-10 | 6.5 h | 2.8 h | 1.2 h |
| 1e-10 | 19.5 h | 5.7 h | 3.7 h |

- **My estimate (upper bound).** I take the step at 2.8τ, which was observed when the step froze. τ = ρK/μ is the penalty time scale. I assume the Study 2 window runs 0.08 → 0.21 s, which gives 788 steps in 10 h and 45.7 s/step.
- **Numerical's estimate.** They assume the K = 1e-11 run started at t = 0 and kept a constant sim-s/h rate. That implies 3.5 s/step if the step had been 8.5 µs throughout, which is not credible for MUMPS at 1.58 M DOF. So the step froze only late in that run, and the run's rate was not constant.
- **The truth lies between the two.** The K = 1e-9 solver log settles it in 5 minutes: count the steps and read the time per step. Record this as task T0.
- **It does not change the plan.** Numerical's own finding C4 is that the K study plateaus on the 0.61 mm mesh. Only restart windows (≤ 3.7 h) and N5 stationary screening (about 1.5 h for 4 K values) are worth running. All of them fit in one night under either estimate.

## 2. O1 re-assessed: what would local refinement cost?
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

## 3. New options: COMSOL steps, engineer-days, compute

### N2: frozen-angle body-fitted hinge box with moving-wall velocity and a fitted lumped law (core)
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

### N1: closed-angle body-fitted leakage (core, clinically dominant)
**COMSOL steps:**
1. Rotate the v20 leaflets to the closed angle. This angle is unknown, so **ask Dr. Bellofiore today**.
2. Check the rim gap at the housing seat. If the CAD shows contact, a seat gap must be set. That is a model value, so it needs advisor approval.
3. Cut along the x = 0 and z = 0 symmetry planes to a quarter model, about 0.6 M DOF. Steady flow at Re_h ≈ 10 makes symmetry plausible (UNVERIFIED).
4. Stationary solve with Auxiliary Sweep continuation on Δp from 0 to 100 mmHg.
5. Age-of-fluid: add a Transport of Diluted Species with source = 1 to get pocket residence time.

**Cost:** 2.5–3 engineer-days, ≤ 3 h compute. **P about 70%.** The main risk is geometry: seat contact and rotating the leaflet by a large angle.

### N3: 2D hinge cross-section with prescribed closing motion (optional)
**COMSOL steps:** a 2D xy section through the ear. If the ear outline is concentric with the pivot, use a Rotating Domain with a circular identity pair (Flow Continuity) in the mid-gap. Otherwise use 2D Moving Mesh with remeshing, which is cheap in 2D. Put at least 10 cells across the gap and use dt of about 1e-5 over 30 ms.

**Cost:** 2–3 engineer-days, about 2 h compute. **P about 80%.**

**Limit:** it misses the end-face and axial leakage paths. It is a good closure-squeeze appendix but **not in the critical path**.

### N4: mid-plane sliding split in the ALE model (stretch, drop for thesis)
**COMSOL steps:**
1. Rebuild v20 with **Form Assembly** so each end gap splits into two 25 µm slabs: a housing half (fixed) and a leaflet half (rigid rotation).
2. On the mid-plane Identity Pair, apply Laminar Flow **Flow Continuity**. Make the destination side finer.
3. In Moving Mesh, the leaflet-half slab gets prescribed rigid rotation, and the rest keeps the current Yeoh/sliding settings.
4. The housing slab must cover the whole swept footprint.

**Cost:** 6–8 engineer-days just to rebuild the 118-boundary geometry and re-attach every selection. Flow-on ALE for one opening is then 20–60 h (UNVERIFIED). It is also UNVERIFIED whether automatic remeshing keeps non-conforming pairs in an assembly. **P about 30%. It does not fit in the calendar.** Leave it for a paper revision or PhD work.

### N5: stationary K screening
Already costed above: about 1.5 h, 0.5 engineer-day. Keep it.

## 4. Responses to schedule claims
- **Challenger: "only A, possibly B; C–E will not pass".**
  - Agreed for D and E.
  - Disagree for C *if it is redefined*: the film is checked against the N2 body-fitted box, not against the ghost model. The film has under 50k DOF, the one-way coupling is justified (Inventor: hinge torque about 1e-3 of blade torque; Advocate: about 0.2%), and it costs about 2 engineer-days.
  - Their "K = 1e-10 about 13 h" sits inside my 5.7–19.5 h bracket, but it is moot because the K study plateaus (§1).
  - Their "GCI beyond 32 GB" is correct for global refinement. It is solved by coarsening (the current mesh becomes the fine grid) plus GCI inside the box.
- **Advocate: two-way mouth exchange plus pocket blocking.** The flux return to 3D needs geometry cuts plus fixed-point iteration, about 3–4 engineer-days. The physics shows the feedback is negligible, so cut it and keep one-way coupling.
- **Inventor: "30 windows at 0.5°".** Compute is trivial (about 5 h) because each window is short. The cost is per-angle remeshing. Let O6's e decide between 16 and 31 angles.

## 5. Revised plan: O4 core = N2 + N1 + O6, with O1-C (film vs box) as the novelty and ghost A + N5 as the cross-check

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

## Bottom line (feasibility ranking, R2)
1. **N2 (within O4):** body-fitted 50 µm gap with GCI, fits in 2–8 GB and about 12 h of compute.
2. **N1:** the clinically dominant phase, about 3 days of work; it depends on the advisor giving the closed angle today.
3. **O6:** 1 day, and its eccentricity value sets the N2 angle spacing and the validity of the moving-wall method.
4. **O1 limited to Phase C (film vs N2 box):** novelty for about 2 days of work. Phases D/E are dropped.
5. **O2 + N5:** a cheap blade-scale cross-check, reported as K-sensitivity only.
6. **N3:** cheap but off the critical path.
7. **N4:** 6–8 days of rebuild at about 30% odds, out of reach.
8. **O3:** a kinematics appendix only.
9. **O5:** infeasible.
