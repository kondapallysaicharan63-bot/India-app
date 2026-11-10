# R1 — Feasibility Engineer (COMSOL 6.4, 4-core / 28–32 GB desktop, ~3.5 weeks)

Scripts: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/feasibility/sched.py` and `plan.py`.

## 0. What the measured run times tell us
I took these numbers from the handoff and checked them in Python:
- **Ghost cost per step.** The K = 1e-11 run did 788 steps (0.0067 s of simulated time at dt = 8.5e-6 s) in 10 h, which is **45.7 s per time step** at 1.58 M DOF with MUMPS. The K = 1e-9 Study 2 run took 80 min, so it used about 105 steps over 0.13 s (mean dt about 1.24e-3 s). Step size limits the cost. The cost of a single step does not.
- **Penalty time scale.** τ = ρK/μ is 3.0e-4, 9.1e-5, 3.0e-5 and 3.0e-6 s for K = 1e-9, 3e-10, 1e-10 and 1e-11. In the 1e-11 run the step froze at about 2.8τ. With dt ≈ min(1.24e-3, 2.8τ), the full 0.13 s window would take **1.3 h (measured) / 6.5 h / 19.4 h / 194 h**. A 0.025 s restart window around peak flow (t = 0.11 → 0.135 s) would take **0.3 / 1.2 / 3.7 / 37 h**. K = 1e-11 is infeasible on any calendar. The K study has to use restart windows.
- **ALE (O3).** 0.654 s and 9.2° in 25.2 h, flow OFF. That is **0.026 sim-s/h, or 0.36°/h**, well below the 0.06–0.1 sim-s/h figure quoted in the task. Remeshing cost goes with the angle swept, not with simulated time. A real opening (15° in about 50 ms) still needs roughly 15+ remesh restarts, and those were 2–6 h each once out-of-core set in.
- **GCI meshes (refinement ratio r = 1.3, 3D).** Coarse 0.72 M DOF, fine 3.46 M DOF. With MUMPS, memory grows about N^(4/3) and time per step about N² (UNVERIFIED scaling rule), so the fine mesh needs **about 2.9× the memory and 4.8× the time per step**. The ALE run already went out-of-core on 28.7 GB. Assume the fine mesh runs out-of-core on the desktop (×3 slowdown, as seen in ALE intervals 5–6), or use HPC or an iterative solver (GMRES + SCGS/AMG) for it.
- **Added mass.** A flat-plate estimate (chord about 10 mm, span 21.6 mm, rotating about an edge) gives I_added ≈ 5.1e-8 kg m², which is **about 5.3× the leaflet inertia (9.5e-9 kg m²)**. The geometry here is assumed, so this is UNVERIFIED. With a ratio this far above 1, segregated or fixed-point coupling between the leaflet ODE and the fluid will diverge unless it is relaxed or stabilized. This is the core risk of Phase E.

**Main finding:** for the phased plan, compute is not the bottleneck. Engineer time is (section 3).

## 1. Option-by-option implementation and risks

### O1 Hybrid ghost leaflet + Reynolds thin-film hinge
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

### O2 Ghost leaflet alone
Same as O1 without the film; this is the existing v22 model.
- **Risks:** with sqrt(K) ≈ 32 µm against a 50 µm gap, the hinge is *not* resolved. The local mesh needs about 5 elements across the gap, which is unaffordable in a fixed mesh. The gap flux will be wrong, but blade torque within 5–9% is plausible.
- **Cost:** the K study + GCI above, about 45 h total. **Feasible and low-risk, but it says nothing new about the hinge.**

### O3 Body-fitted ALE with automatic remeshing
Moving Mesh (Yeoh), prescribed mesh velocity + Nitsche, Sliding on 18 faces, remeshing at 1.5. These settings are proven.
- **Blocker:** inverted elements on the ear end faces (z = ±10.8 mm) from about 3° onward. The likely fix is to give the 50 µm end gaps an Identity Pair or a separate swept ring domain that rotates rigidly. That is a geometry rebuild, 3–5 days, uncertain.
- **Cost:** with flow ON, expect 2–3× the flow-off rate (UNVERIFIED), which is ≥ 60–120 h for one opening. Two-way FSI adds the Global ODE plus added-mass instability.
- **Not feasible as a primary result.** Keep it as a kinematics/remeshing demonstration (the existing run).

### O4 Two-scale
- **(i) Blade-scale model:** the ghost leaflet (O2) or the pin-joint v21 for the FSI opening history θ(t).
- **(ii) Hinge model:** the real-hinge body-fitted model at fixed angles (already exists, v20_fixedangles), with the 3D pressure from (i) as boundary conditions. Optionally, the Thin-Film Flow model in the same pocket for comparison.
- **Risks:** low, because every piece already runs. The weakness is one-way coupling (hinge → blade feedback is ignored, but that is small since the gap flux is ≪ the main flow; UNVERIFIED magnitude).
- **Cost:** reuses Phase A/B runs plus about 1 h per angle of hinge runs. **Highest probability of a finished, defensible result.**

### O5 Manual overset
- **How it would be built:** two components (a rotating leaflet-body mesh and a fixed background mesh), with fringe interpolation through General Extrusion both ways and hole cutting by χ.
- **Why it fails here:** there is no built-in overset, so the Jacobian across the fringes is missing in both directions, and the overlap must still fit inside the 50 µm gap. DOF is about 2×, so out-of-core is certain.
- **Development time:** ≥ 3 weeks. **Infeasible in this calendar.**

### O6 Pivot optimization / rotation + sliding
- **Geometry check:** sample the ear and pocket surfaces, rotate them about a candidate pivot, and minimize the maximum penetration and the clearance variance. Python (scipy.optimize) on exported STL/mesh points, or a COMSOL parametric sweep with `minop` distance probes.
- **Cost:** about 1 day, 2 h compute.
- **What it gives:** it explains the ear-face inversion, and it tests whether pure rotation is kinematically admissible. If it is not, the ALE needs a translational DOF (prescribed mesh displacement = R(θ)x + s(θ)). **Cheap, novel, do it.**

### O7 (new) Quasi-steady stationary ghost solves for the K and mesh trend
- **What:** stationary solves at the peak-flow boundary conditions (t = 0.13 s) avoid the τ step limit entirely. About 10 Newton iterations × 46 s per solve, roughly 15–30 min each (UNVERIFIED that steady converges at Re near peak).
- **Use:** as a *screening* tool for the K → 0 extrapolation (Richardson in sqrt(K)). Confirm only the final pair with transient restart windows.

## 2. Solver notes (all options)
- **Linear solver:** keep MUMPS and always check LinRes at about 1e-14. Do not use cuDSS: on the T4 it silently gave LinRes = 1.
- **Penalty problems:** fully coupled, Jacobian update "On every iteration" (the ALE fix also applies here).
- **K ≤ 1e-10:** consider excluding the penalty-region velocity from the time-error estimate, or setting a larger absolute tolerance there. Otherwise BDF tracks the τ transient (UNVERIFIED that this relaxes the step limit without hurting the torque; test against the K = 1e-9 run).
- **Fine GCI mesh:** run on HPC, or use iterative GMRES with an AMG/SCGS preconditioner on the desktop.

## 3. Schedule arithmetic

### Compute
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

### Engineer time
This is the real limit.
- A: 1 day, B: 2 days, C: 3–4 days, D: 3–4 days, E: 4–6 days, O6: 1 day.
- **Total 14–18 working days**, against about 13 working days available (Nov 10–27, minus weekends and Thanksgiving Nov 26).
- **E does not fit with margin.** E is therefore a stretch goal. A+B+C+O6+O4 fits (about 8–9 days), and D fits if C passes by about Nov 19.

## 4. Gantt (best feasible plan: A → B → C → O4, with O6 in parallel; D is conditional; E is stretch)

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

## 5. Uncertainties (UNVERIFIED)
- How Thin-Film Flow connects to 3D Laminar Flow in 6.4 (no built-in coupling found), and how its Border node behaves.
- Whether General Extrusion destination-map derivatives reach the Jacobian when the map depends on a Global ODE variable.
- Whether excluding variables from the time-error estimate relaxes the τ step limit safely.
- The body-fitted reference run time and memory at the fine mesh.
- The added-mass ratio, which needs the real leaflet chord and pivot location.

Documentation found: the COMSOL 6.4 "Thin-Film Flow Interface (3D)" page (https://doc.comsol.com/6.4/doc/com.comsol.help.cfd/cfd_ug_fluidflow_thinfilm.12.03.html) and the 6.4 "Common Settings for Nonlocal Couplings" page (https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_definitions.21.104.html), both seen only as search snippets. doc.comsol.com was blocked in this session, so I could not read the full text. The wall/base velocity options I cite come from the 5.3 Fluid-Film Properties page (https://doc.comsol.com/5.3/doc/com.comsol.help.mems/mems_ug_fluidflow.07.07.html) and may differ in 6.4. No journal literature was checked: needs literature check.

## Bottom line (feasibility ranking)
1. **O4 Two-scale.** Every component already runs, and it fits with days to spare. This is the guaranteed publishable core.
2. **O1 Hybrid (phases A–D).** Compute is trivial (about 25 h). Phases C/D are 6–8 days of development and the novel part. Phase E (two-way FSI) is a stretch because of the missing ∂χ/∂θ Jacobian and an added-mass ratio of about 5.
3. **O2 Ghost alone.** Cheap and low-risk (about 45 h for K + GCI), but it cannot resolve the 50 µm gap (sqrt(K) ≈ 32 µm). Treat it as the validation layer for O1/O4, not as the result.
4. **O6 Pivot optimization.** About 1 day and almost no compute. It explains the ALE inversion and adds a kinematic novelty point.
5. **O7 Stationary K screening.** A cheap accelerator for Phase B. Not a result on its own.
6. **O3 Body-fitted ALE.** More than 60–120 h per opening with flow, and the ear-face inversion is unsolved. Report it as a kinematics demonstration only.
7. **O5 Overset.** No built-in support, at least 3 weeks of development, and still needs an overlap inside the 50 µm gap. Infeasible before defense.
