# R2 — Inventor rebuttal: defending N2 (box submodel + lumped law) and N4 (end-gap sliding split)

Numbers come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/inventor/r2.py`, plus `calc.py` from R1.
Geometry not taken from the handoff (leaflet stump size, gap length L, window run time) is marked **UNVERIFIED**.

## 1. Where the openings agree, and what that does to the debate
Challenger, Numerical, Clinical and Feasibility all rank O4 first or second. Challenger's O7 is the same construction as my N2. Numerical C9 independently proposes "end-face slip with a gap-local mesh that rotates rigidly", which is N4.

So the open questions are not "two-scale or not". They are:
- (i) Is the moving-angle part of N2 accurate and mass-consistent?
- (ii) Is it novel beyond Simon et al.?
- (iii) How does it relate to O1's film?

## 2. Objections to N2, answered

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

## 3. How N2 and O1 relate: a model hierarchy, and the combined novelty
Integrate the Reynolds film across the gap width for a uniform gap:

  R_film = 12μL/(w h³),   C_film = (couette term from U_t h/2)

These are the linear terms of N2's law. **N2's 0D law is the film PDE projected to 0D. The film is N2's law resolved in space. The 3D box is the truth for both.** What the film lacks is exactly what K(θ) captures in 0D: losses at the mouths and corners. Finite-width side-wall effects are only about 0.5–0.8% (0.63·h/w).

The combined contribution is a **three-level hierarchy, 3D body-fitted → 2D film → 0D law, with the error measured at each level and at each angle**. N2 is O1's missing validation reference; Challenger and Numerical both demand one. Two outcomes are possible:
- If the film matches 3D within 10% in the uniform region, O1's film can carry diastolic leakage.
- If it does not, the fitted mouth loss K(θ) can be added to the film's boundary conditions as a 3D-calibrated mouth correction, which is still novel.

Either outcome is publishable. **The ghost-leaflet half of O1 is not needed for this.** Numerical C2/C3 and Challenger 2.2 show the ghost cannot hold the closure pressure: 1.2–3.9× gap flux and m/s-scale leakage at K = 1e-9.

## 4. Clinical's N1 and N3
- **N1 (closed-angle diastolic leakage)** is N2's first and most valuable case, at θ = θ_closed and θ̇ = 0. The geometry is fixed, so the window error is zero and only the mesh GCI applies. I fully endorse making it **deliverable #1**.
  - Prerequisites: the closed-stop angle from the advisor or the CAD (Clinical: UNVERIFIED, about 25–30°), and a Carreau sensitivity for residence time in the pocket.
- **N3 (2D hinge section with prescribed closing)**: I propose to repurpose it as **the falsification test of N2's window error**.
  - In 2D, true ALE with remeshing is cheap and does not suffer the end-face problem.
  - Run the same kinematics in 2D with full ALE and with frozen windows. The difference is the frozen-geometry error measured directly, with no other modelling difference.
  - N3 then also delivers Clinical's closure-squeeze shear story.
  - Its limit: no end gaps (z), so it says nothing about N4.

## 5. Falsifiable validation protocol for N2 (numeric thresholds)
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

## 6. N4, defended and downgraded
- In-plane rotation of a plane normal to the pivot axis is pure sliding, so the end-gap shear drops from 0.67–0.77 to 0 per degree.
- Numerical C9 says the inversion comes from **sliding nodes meeting pinned nodes at corners**. N4 moves that corner to the place where the pair face meets the pocket's lateral wall.
  - If the lateral clearance is also about 50 um and not concentric, the same in-plane shear reappears in 2D there.
  - N4 therefore works only if the lateral thin arcs are concentric to about 5 um, so that they can become cylindrical sliding pairs too, or if the lateral clearance is much larger than 50 um.
- Flow Continuity on a non-conforming pair is enforced weakly, so the flux mismatch across the pair must be under 1% of the gap Q.
- **I lower P(N4) from 40% to about 30%.** It is a stretch goal after O6 has measured the concentricity.

## Updated ranking
1. **N2 (with N1 as its first case and O4 as the global driver)**: body-fitted at 50 um, local GCI is affordable, the error bounds are falsifiable, and the hierarchy gives the novelty.
2. **O6**: supplies e and h(θ), which set N2's window count and the clinically dominant leakage (∝ h³). Run it first.
3. **O1, film part only, validated against N2**: the combined hierarchy novelty. Drop the ghost-leaflet bulk for closure.
4. **N3 (repurposed as V3)**: cheap, directly falsifies the frozen-geometry error, and adds closure shear.
5. **O2**: the blade-scale bulk flow only, with a K-sensitivity check rather than a convergence claim.
6. **N4**: high payoff only if the arcs are concentric; ~30%.
7. **N5**: an accelerator for O2's K screening.
8. **O3 as it stands**: kinematics appendix only.
9. **O5**: infeasible in the time available.
