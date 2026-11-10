# R2: Numerical Analyst verdicts on the Round-1 numbers

Source: `bme298_debate/math_checks.py`, sections R2-C1 to R2-C8, appended to the Round-1 checks. The whole file runs cleanly, and its output is in `math_checks_output.txt`. Assumptions are as in R1. I did not check any literature this round; every literature item still needs a check.

## 1. Leakage conflict (Challenger 29.7x vs Numerical 3.9x): no conflict (R2-C1)

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

## 2. K-study run time (R2-C2): I withdraw my R1 estimate

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

## 3. Inventor (R2-C3)

**Frozen-angle error.** The Reynolds film has a pressure-driven flux q ∝ h³. An axis offset e moves the gap arc locally by dh = e·dθ·cos φ, so δQ/Q = (1 + dh/h)³ − 1 ≈ 3e·dθ/h. **PASS**: the table is reproduced (e = 10 µm: 0.26 / 0.52 / 1.05%; e = 25 µm: 0.65 / 1.31 / 2.62%; e = 50 µm: 1.31 / 2.62 / 5.24%, exact 5.33%). The 5% limits of 4.8°, 1.9° and 0.95° are reproduced. Two caveats:
- It is a local upper bound. cos φ changes sign along the arc, so the net flux error over a long arc is smaller.
- The WSS error, which is linear in h, is one third of the flux error.

**Back-coupling.** D = μR²A/h = 2.0e-9 N m s per hinge, against Inventor's 2.2e-9. The friction torque relative to 0.2 mN m is 0.9e-4 for the mean closing motion and 1.0e-3 at peak. With both ears it is 1.8e-4 to 2.0e-3. **PASS.** One caveat: 0.2 mN m is the systolic torque. The closing torque is not known, and if it is larger the ratio only gets smaller.

**Leakage.** With w = 6 mm, L = 1 mm and 80 mmHg, the leak is **0.38 mL/beat**. Clinical's w = 4 mm, L = 2 mm and 100 mmHg give 0.16 mL/beat. **Both PASS.** The 2.4x difference comes only from the unverified w/L ratio, so I would report 0.15–0.4 mL/beat until w and L are measured in the CAD. The gap sensitivity is −49% at 40 µm and +73% at 60 µm (PASS).

**Other Inventor numbers.** All PASS: 0.77 shear per degree, 11.5 at 15°, 576 µm of slide, h²/ν = 0.76 ms, and a 3.8% entrance-loss share.

## 4. Clinical (R2-C4)

**Stress accumulation per pass.** At 100 mmHg and L = 2 mm: τ_w = 167 Pa, U = 0.40 m/s, transit time 5.04 ms, so τ_w·t = **0.84 Pa s. PASS.** About 4.2 passes reach 3.5 Pa s. Two things are new here:
- **τ_w·t = 6μL/h, which is independent of dp.** Linear stress accumulation per pass in a slot depends only on L/h. A narrower gap gives more stress accumulation per pass, which matches Clinical's washout argument.
- 0.84 Pa s is the conservative bound. The flux-weighted mean over streamlines is 3μL/h = **0.42 Pa s, or about 8 passes**. Report both, or use pathlines.

**Giersiepen power-law damage.** HI = 0.133% and platelet index = 0.383%. PASS, with the constants still needing a literature check.

**Carreau.** Clinical's percentages (0.7 / 2.2 / 8 / 20 / 37 / 57%) are reproduced, but they are measured against μ∞ = 3.45 mPa s. Against the model's 3.5 mPa s:
- In the gap the Newtonian model is −0.8%, a slight over-prediction.
- In the pocket at 250 / 100 / 50 s⁻¹ the error is +18 / +35 / +55%.

The conclusion does not change: Newtonian is fine in the gap, and a Carreau case is needed for residence time.

## 5. Inventor N4: split each end gap at its mid-plane with a sliding interface (R2-C5)

**Shear.** **PASS, with a condition.** A rotation about z maps every plane z = const onto itself (numerically |dz| = 0). The end-gap mesh shear therefore drops from 11.5 to 0, and the in-plane slide is carried by the identity pair. The condition is that the rotating slab must be bounded laterally by a **cylinder about the pivot axis**, not by the ear outline. Otherwise its side faces sweep through the fixed slab. In other words, N4 is a thin rotor disc of radius at least 2.3 mm with planar and cylindrical sliding pairs, which is rotating-machinery topology.

**5 µm concentricity tolerance.** PASS as adequate. A cylindrical interface placed inside a lateral gap survives the whole sweep if e·θ + 2σ_r < h/2 − h/6 = 16.7 µm (here e is the pivot-axis offset and σ_r the surface scatter). σ_r = 5 µm then allows e up to about 25 µm. σ_r = 10 µm fails even at e = 0.

**Mesh cost.** **FLAG.** A 25 µm slab with 3 isotropic tet layers needs about 6 M elements per slab, and there are 8 slabs. Swept prisms (3 layers of 75 µm triangles) need about 20 k. **N4 is feasible only with swept or prism meshes.**

## 6. Corrected validation statistics (R2-C6)

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

## 7. Other conflicts (R2-C7)

- **Halo thickness.** Halo = ΔV·t/(2V): 152 / 196 / 240 µm for t = 0.7 / 0.9 / 1.1 mm. Advocate's ~200 µm (t = 0.9 mm) and my 150–240 µm come from the same formula. Either way the halo is 3–5 gap widths, so the hinge is sealed. PASS.
- **Added inertia.** I_a = (9π/128)ρc⁴s gives 5.3x for Feasibility's chord and span, and 1.3–6x across plausible geometry. All three agents agree within the geometry spread. **Measure it directly:** in the body-fitted model, prescribe α with the flow at rest and take I_a = −T/α.
- **Challenger spot checks.** All PASS: Couette 1.3 Pa, gap Womersley 0.076, and GCI memory 2.9x under N^(4/3) scaling.

## 8. Proposed numeric gates for the Judge (R2-C8)

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

## Bottom line (ranking, updated)

1. **N2 (real-hinge box submodel) together with O4 (two-scale):** body-fitted, so the leak does not arise (R2-C1). Its frozen-angle error is bounded and computable (R2-C3), one-way coupling is justified to about 1e-3, and a local gap GCI is affordable.
2. **N1 (closed-angle leakage case):** the only route to the clinically dominant diastolic numbers. Stress accumulation is 0.42–0.84 Pa s per pass, and leakage is 0.15–0.4 mL/beat.
3. **O6 (pivot optimization):** it supplies e, which sets gates G10 and G14 and the leakage sensitivity (h³).
4. **O1 (hybrid ghost + film):** valid only for bulk flow and the film, with the gap taken from the film. Phase A fails until FixChi is run; the K study is feasible in 0.025 s windows.
5. **N4 (mid-plane split), as an O3 rescue:** geometrically exact, but it needs rotor-disc topology and swept meshes.
6. **N5 (stationary K screening) and N3 (2D closure model):** cheap supporting tools.
7. **O2 (ghost alone):** disqualified in the gap at every K that can be run.
8. **O3 as it stands:** kinematics only.
9. **O5 (overset):** too costly in degrees of freedom.
