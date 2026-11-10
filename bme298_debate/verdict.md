# Judge verdict: BME 298 real 50 um hinge-gap method selection

Judge: independent (did not argue). Written Tue Nov 10, 2026.
Inputs read: HANDOFF_BME298_2026-11-10.md, BRIEF_common.md, debate_transcript.md (Rounds 1 and 2), math_checks_output.txt.
Judge-run checks: `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/judge/judge.py`. Results are cited as [J1]–[J5] below.
- [J1] A/B asymmetry recomputed: 1.04 / 0.62 / 1.80 / 2.74 %.
- [J2] TOST power by simulation, n = 4, margin ±5 %, alpha 0.05.
- [J3] Slope-CI behaviour at n = 4.
- [J4] Bland–Altman LoA CIs.
- [J5] Elements left across the gap after global coarsening.

Option labels follow the canonical Round-2 map (O1–O6, N1–N5).

Literature caveat (applies to every paper named in this file). The transcript reports these as found only at abstract or snippet level, because full-text hosts were blocked in the session:
- Simon et al. 2010 (two Ann Biomed Eng papers)
- Jun et al. 2014
- Dumont et al. 2007
- Govindarajan et al.
- Klusak & Quinlan
- Hester, Vasil & Burns
- Angot, Bruneau & Fabrie 1999

Every one needs a full-text literature check before it is cited in the thesis. I add no citations of my own.

---

## 0. Where the debate converged, and where it converged too easily

By Round 2 the team agreed on five points:
- The ghost leaflet cannot be used inside the gap or at closure (R2-C1, C2, C3).
- One-way coupling between hinge and blade is justified (R2-C3: friction torque is 1e-4 to 2e-3 of blade torque).
- The core should be body-fitted: O4 + N2 + N1, with O6 first.
- The Reynolds film becomes a validated reduced model (O1-R), not a coupled solver.
- Phases D and E are dropped.

I agree with the direction. Six points were accepted too readily:

1. **Nobody has checked that the body-fitted reference actually resolves the 50 um gap.** Every agent treated "body-fitted v20" as the truth. The number of elements across the gap in v20_fixedangles is not stated anywhere.
   - If it is 1–2 elements, the "reference" gap flux and WSS are themselves unresolved.
   - Global coarsening by r² = 1.69 leaves 1.8 / 2.4 / 3.0 elements from 3 / 4 / 5 [J5].
   - This is now Gate G0b, checked today.
2. **N2's moving windows are under-specified.**
   - No validated global model gives θ(t) or box boundary data through closure. The torque table covers only 70–84°, and v21 is a stationary pin-joint.
   - The Challenger (R2 §3(iv)) and Clinical (R2 §3) both noted this, but the Feasibility Gantt still schedules the windows.
   - I make them a stretch goal.
3. **N1's P ≈ 70 % (Feasibility) ignores the 55° rotation of the v20 ears into the closed position.** That rotation is exactly where O6 may show interference, and it adds a possible seat-contact model value that needs the advisor. I use about 55 %.
4. **"Box vs full v20" (Inventor V1) is a self-consistency check, not independent validation.** The thesis has no experimental validation at the hinge. The independent anchors are:
   - the analytic slot solution;
   - order-of-magnitude literature ranges (Simon 612 Pa diastolic peak; Dumont 0.81 %);
   - PIV data, if the advisor has any.

   The thesis must say this plainly.
5. **n = 4 statistics are close to theatre.** The model-vs-model error is deterministic and depends on angle. It is not sampling noise.
   - The primary decision rule should be deterministic, per angle, against numerical uncertainty (a V&V-20-style E vs u_val; needs literature check).
   - TOST/BH are secondary, for the committee.
   - The Numerical analyst's G4 (slope 95 % CI inside [0.95, 1.05]) **fails even for an unbiased model** at n = 4: CI [0.956, 1.056] for a ~0 % bias case [J3]. It is as non-discriminating as R², only in the opposite direction.
6. **Is 50 um the physiological clearance, or only the CAD value?** The advisor rule ("model values only if physiologically correct") applies here.
   - Clinical shows leakage ∝ h³: ±10 um gives −49 % / +73 % (R2-C3).
   - The 40/50/60 um sensitivity is therefore mandatory, not optional.

---

## 1. Ranking of all options with probability of success

P is my judgement, informed by the transcript estimates and adjusted for the points above. "Success" is defined per option. All deadlines are on or before code freeze, **Fri Nov 27**.

| Rank | Option | Success means | P(success) | Justification (transcript / math checks) |
|---|---|---|---|---|
| 1 | **O4 Two-scale framework** (blade-scale loads from existing body-fitted fixed-angle runs + real-hinge submodel) | Thesis chapter with body-fitted real-hinge results at 70/75/80/84° with bulk GCI (G5) and a justified one-way coupling (G11), in the thesis by Nov 27 | **80 %** | Every component already runs (Feasibility R1 §1 O4). One-way coupling verified at 1e-4 to 2e-3 (R2-C3). Monolithic, so no added-mass issue (C12). Bulk GCI is feasible only by coarsening (C11). |
| 2 | **N2-lite**: real-hinge box at frozen angles 70/75/80/84°, wall velocity 0 (fixed angle), local gap GCI | V0 mass balance, V1 box-vs-full ≤ 3 %, V2 gap GCI < 5 % (flux, mean WSS) passed by Nov 18 | **70 %** | Fits in 2–8 GB (Feasibility R2 §3). Local GCI is affordable (C11, R2-C8 G5b). Risk: boundary-layer meshing in a sliver gap at some angles (Challenger R2 §3(iii)). |
| 2b | N2 moving windows (frozen-geometry transpiration through a swing) | V4 δθ self-convergence < 2 % plus V0 on the stump-flux-corrected box, by Nov 25 | **30 %** | No validated θ(t) or box BC through closure (Challenger R2 §3(iv)). Stump swept flux 0.1–1.1 mL/s against a hinge Q of ~0.19 mL/s (Inventor R2 §2a). Window count depends on the unmeasured e (R2-C3, G10). |
| 3 | **O6 Pivot fit / rotation+sliding** | e, σ_r and the min-clearance map h(θ) for 70–85° **and** along the swing to the closed angle, from CAD point clouds, by Nov 12 | **85 %** | About 1 day, Python on exported surfaces (Feasibility R1 O6). Sets G10, G14 and the leakage input (∝ h³). Explains the ALE inversion (C9). |
| 4 | **N1 Closed-angle body-fitted diastolic leakage** | Closed-angle case with gap GCI < 5 % on Q and mean WSS, mass imbalance < 1 %, uniform-gap WSS within ±30 % of slot analytic, Carreau sensitivity, by Nov 24 | **55 %** | Clinically dominant phase (Clinical R1 §1, R2 §5). Steady, reuses the v20 workflow. Risks: closed angle UNVERIFIED (~30° inferred), seat contact, ear interference across the 55° swing (point 0.3). |
| 5 | **O1-R** (Reynolds film on CAD h(s), validated against N2/N1; no ghost coupling) | G8 valid and G9 film-vs-3D flux within 10 % and mean WSS within 15 % at ≥ 3 of 5 cases, by Nov 25 | **45 %** | Physics sound in the uniform gap (C8: eps·Re 0.01–0.08 systole, 0.12–0.60 diastole). Mouths and corners are outside validity. Thin-Film ↔ 3D coupling is UNVERIFIED in 6.4 (Feasibility R1 §5). A negative result with a validity map is still reportable (Advocate R2 §7). |
| 6 | **O2 Ghost leaflet alone**, fixed angles, systole only (FixChi, Phase A cross-check) | G1 volume < 1 %, then torque within ±5 % and Fy within ±10 % at all 4 angles by Nov 14. Otherwise reported as a documented negative cross-check | 35 % (pass); 85 % (documented either way) | Current run fails deterministically: −8.69 % torque, +25.6 % Fy (R2-C6a). Halo sealed the hinge (C5, R2-C7). Disqualified in the gap at every runnable K (R2-C1: 3.6x slip flux, 29.7x with Darcy core) and at closure (C2). |
| 7 | **N3 2D hinge section, prescribed closing**, repurposed as the V3 falsification of the frozen-window error | 2D ALE vs frozen windows at δθ = 1/0.5/0.25°, error ≤ bound and slope 0.6–1.4, by Nov 25 | 50 % (if started); in thesis only if slack | Cheap: 2–3 engineer-days, ~2 h compute (Feasibility R2). Misses the butterfly recess and end gaps (Challenger R2 §4). Off the critical path. |
| 8 | **N5 Stationary K screening** | Stationary solves converge at t = 0.13 s boundary conditions for 3 K values, by Nov 13 | 55 % | Steady convergence at peak Re is UNVERIFIED (Feasibility R1 O7). Useful only if O2 is kept. Reports sensitivity, not convergence (C4). |
| 9 | **O3 ALE + remesh** | Kinematics/remeshing appendix: log saved, LinRes ≈ 1e-14 checked, inversion explained by C9 | 90 % (appendix); < 5 % (flow-on FSI result) | Inversion is geometric: 37–40 um/deg slide at the end faces would need 20–45 remeshes (C9). Rate 0.026 sim-s/h (C10). |
| 10 | **N4 End-gap mid-plane sliding split** | Rebuilt v20 assembly runs flow-off ALE to 15° with no inversion and pair flux mismatch < 1 % | 10 % by Nov 27 (~30 % with unlimited time) | Exact for end-gap shear (R2-C5), but needs rotor-disc topology and swept prisms. 6–8 engineer-days to rebuild 118 boundaries (Feasibility R2). Needs σ_r ≤ 5 um (G14). Post-defense / paper revision. |
| 11 | **O1 as originally specified** (ghost + film coupled, Phases A–E) | Phases A–E pass their gates | < 10 % | Ghost leaks at closure (C2). Slip/Darcy contamination of the mouth Δp (Challenger R2 §2c). No dχ/dθ Jacobian with added mass 1.3–6x (C12). Engineer time 14–18 d vs ~13 d (Feasibility R1 §3). Both D and E were withdrawn by the Advocate (R2 §4). |
| 12 | **O5 Manual overset** | Working overset in COMSOL | < 5 % | No built-in support. Overlap still needed inside 50 um. About 2x DOF beyond memory (C11). ≥ 3 weeks of development (Feasibility R1). |

---

## 2. Recommended plan: phases, dates, measurable go/no-go gates

**Core:** O4 framework = N2-lite + N1, with O6 first.
**Novelty layer:** O1-R film plus the 3D → film → 0D hierarchy.
**Optional cross-check:** O2/FixChi, time-capped.
**Stretch:** N2 moving windows, N3.

### Curated gate table (from Numerical R2-C8, with fixes)

| Gate | Metric | Threshold | Change vs R2-C8 |
|---|---|---|---|
| G0 | LinRes (direct solver) on every step of every run used | < 1e-10 (expect ~1e-14) | kept |
| **G0b** (new) | Elements across the 50 um gap in v20_fixedangles and in every box mesh | ≥ 4 on the finest mesh used as reference; ≥ 3 on the coarsest GCI mesh | new: the reference resolution was never checked |
| G1 | Ghost volume vs CAD 208.5 mm³ | < 1 % per leaflet | kept (optional track) |
| G2 | A/B torque symmetry, body-fitted, **4 significant figures** | < 3 % at every angle | kept; current 2.74 % at 84° [J1] is not decidable at 2 sig figs (1.4–4.1 %, C6) |
| G3 | Ghost torque bias (A/B averaged per angle, n = 4) | **Primary:** \|E_i\| < 5 % at every angle. **Secondary:** TOST ±5 %, p < 0.05 | primary changed to deterministic per-angle |
| G3b | Ghost FyA bias | \|E_i\| < 10 % at every angle; TOST ±10 % | kept |
| G4 | Regression slope/intercept with 95 % CI | **Report only, not a gate** | **retired as a gate**: at n = 4 the slope-CI half-width is ~0.04–0.05, so [0.95, 1.05] fails even at ~0 % bias (CI 0.956–1.056) [J3] |
| — | R² | **Report only, never a gate** | R² = 0.9998 for −3 %, −4.3 % and ~0 % bias alike (R2-C6, [J3]); C7: a pure −8.6 % bias gives R² = 1.000 |
| G5 | Bulk GCI_fine on blade torque, current mesh = finest (coarsen r = 1.3 → 0.72 M, 0.33 M DOF) | GCI < 5 %; asymptotic ratio 0.9–1.1 | kept (C11) |
| G5b | Gap GCI in the box (4/6/9 elements across h, bulk fixed) | Q_h, mean WSS < 5 %; 99th-percentile WSS < 10 %; asymptotic ratio 0.9–1.1 | merged with Inventor V2 and Clinical R2 §5 (99th percentile, not max) |
| G6 | Ghost K sensitivity, \|T(1e-10) − T(3e-10)\|/T, 0.025 s restart windows | < 2 %, otherwise labelled "K-limited" | kept; it is a sensitivity, not a convergence claim (C4) |
| G7 | Ghost used for any gap or closure quantity | **Never** (12KT/h³ and the slip excess are both > 5 % at every runnable K) | stated as a prohibition (R2-C1) |
| G8 | Film validity at evaluation points | eps·Re < 0.3 and h²ω/ν < 0.3; points that fail are reported in the validity map and excluded from G9 | kept |
| G9 | Film vs 3D body-fitted, uniform-gap region | Q within 10 %, area-mean WSS within 15 % | kept |
| G10 | Frozen-window error 3e·δθ/h (e from O6) | < 5 %, which sets δθ | kept (stretch only) |
| G11 | Hinge friction torque / blade torque | < 1e-2 | kept, but **already passed analytically** (2e-3, R2-C3). Re-measure once in the box; it is not a decision gate |
| V0 | Box discrete mass balance | < 0.5 % of Q_hinge (fail ≥ 1 %) | from Inventor R2 §5 |
| V1 | Box vs full v20 at 70° and 84° (Q_h, Δp_h, mean and 99th-percentile WSS) | ≤ 3 % pass, > 5 % fail (box enlarged in between) | from Inventor R2 §5 |
| G12, G13 | Phase D mass balance / Phase E stabilization | **retired**: Phases D/E dropped (Advocate R2 §4) | retired |
| G14 | N4 geometry, e·θ + 2σ_r | < 16.7 um | stretch only |
| G15 | Carreau vs Newtonian change in pocket residence time or WSS | > 5 % on residence time puts Carreau in the main results (Clinical uses 10 %; I take the stricter 5 %) | kept |

### Phases

**Phase 0: Triage. Tue Nov 10 – Wed Nov 11.**
Tasks:
- Save the ALE log; check LinRes in every interval.
- Count the steps and s/step in the K = 1e-9 log (T0, Feasibility R2 §1).
- Re-export the body-fitted torques and FyA at 70/75/80/84° to 4 significant figures.
- Measure the elements across the gap (G0b).
- Verify the sign convention of the B-leaflet torque expression (C13, UNVERIFIED).
- Advisor email (§8).
- Overnight: run FixChi + Torque at 70° (optional track, scripts already exist).

Go/no-go on Nov 11:
- G0 on the v20 reference: pass.
- G2 < 3 % at 4 significant figures.
- G0b ≥ 4 elements across the gap.

**Phase 1: O6 pivot fit. Wed Nov 11 – Thu Nov 12.**
Fit the pivot to CAD point clouds of the ears and pockets (scipy). Outputs:
- e
- σ_r
- min clearance along 70→85° and along the swing to the closed angle
- whether pure rotation is admissible (max penetration ≤ 0)

Go/no-go on Nov 12:
- Clearance > 0 along the full swing under pure rotation (or rotation + slide).
- Record e. If e ≤ 50 um, moving windows remain a viable stretch (≤ 22 windows per 15°, Inventor R2 §2b). If e > 50 um, the moving windows are dropped.

**Phase 2: N2-lite box + bulk GCI. Thu Nov 12 – Wed Nov 18.**
Tasks:
- Build the 6×6×4 mm box with Rotate(th_k) and boundary-layer meshing at 70/75/80/84°.
- Box boundary conditions from the full v20 run at the same angle.
- Mouth cut-planes for Q_h and Δp_h.
- Run the box GCI at 70° and 84°.
- Run the bulk GCI by coarsening v20 (and the ghost, if that track is alive).

Go/no-go:
- **Nov 16:** V0 < 0.5 % and V1 ≤ 3 %.
- **Nov 18:** G5 < 5 % with asymptotic ratio 0.9–1.1, and G5b passed.

**Phase 3: N1 closed angle. Thu Nov 19 – Tue Nov 24.** (Geometry prep can start Nov 16 once the advisor gives the angle.)
Tasks:
- Rotate to θ_closed and set the seat gap (needs advisor approval).
- Quarter model if symmetry holds.
- Δp continuation 0 → 100 mmHg.
- Age-of-fluid (Transport of Diluted Species with source = 1).
- Pathline stress accumulation.
- Carreau case.
- Parametric gap 40/50/60 um.

Go/no-go on Nov 24:
- G5b passed; mass imbalance < 1 %.
- Uniform-gap WSS within ±30 % of the analytic slot value (133–333 Pa for L = 1–2 mm, C8).
- Leakage 0.05–1 mL/beat for 4 hinges (analytic band 0.15–0.4, R2-C3).
- Fraction of forward-flow stress accumulation above threshold < the regurgitation fraction (Clinical R2 §5).
- G15 evaluated.

**Phase 4: O1-R film validation + hierarchy. Fri Nov 20 – Wed Nov 25.**
Go condition on Nov 20: Phase 2 has passed.
Tasks:
- Thin-Film Flow (or a Coefficient Form Boundary PDE) on the CAD h(s), with mouth Δp taken from the box or N1. Never from the ghost.
- Compute R_film = 12μ∫ds/(w h³) and compare with the box-fitted R(θ).

Go/no-go on Nov 25: G8 and G9 at ≥ 3 of 5 cases (4 systolic angles + closed).

**Phase 5: Stretch (only if Phases 2–3 pass early). Mon Nov 23 – Wed Nov 25.**
Either the N2 moving windows over one opening (requires e ≤ 50 um; G10; V4 < 2 %) or N3 as V3. Not both.

**Freeze: Fri Nov 27** (Thanksgiving is Nov 26).
**Writing: Sat Nov 28 – Thu Dec 3.**

**Optional ghost track (O2), time-capped.**
- Overnight runs Nov 10–13 only: FixChi, then 4 angles, then N5 and 0.025 s K restarts (1.2 h at 3e-10 and 3.7 h at 1e-10, R2-C2).
- Gates: G1, then G3/G3b, then G6.
- **Hard stop Fri Nov 14.** After that it receives no engineer time and goes into the thesis as whatever result it has.

---

## 3. Fallback at each phase

| Phase | If the gate fails | Fallback |
|---|---|---|
| 0 | G0b < 4 elements across the gap in v20 | v20 is a **blade-scale** reference only. The box (Phase 2) becomes the only gap reference, and V1 is restricted to blade-scale/mouth quantities (Δp_h). |
| 0 | G2 ≥ 3 % at 4 significant figures | Report the A/B asymmetry as a numerical-uncertainty floor and add it to u_val. Check the mesh symmetry and the sign of the B torque expression. |
| 1 | O6 shows interference under pure rotation | Prescribe rotation + slide in N1/N2 geometry. Record it as a finding: it explains the ALE end-face inversion and has kinematic novelty. |
| 1 | No CAD point export by Nov 12 | Use the nominal 50 um, uniform. Drop the moving windows. Widen the gap sensitivity to 40/50/60 um. |
| 2 | V1 > 5 % | Enlarge the box once (e.g. 9×9×6 mm). If it still fails, use the full v20 model with local gap refinement at 70° only. |
| 2 | V0 ≥ 1 % | Velocity BCs on all faces but one pressure outlet (Inventor R2 §2a). If it still fails, run fixed angles only (θ̇ = 0 makes the stump flux zero). |
| 2 | G5/G5b fail | Report GCI values and the observed order honestly, and use the medium/fine difference as u_num. Do not claim mesh convergence. |
| 3 | No closed angle from the advisor by Nov 16 | Use 30° (inferred, Clinical R1 §3) with a ±5° sensitivity, clearly flagged as an assumption. If the geometry still fails by Nov 20, apply diastolic Δp at the nearest angle that can be built and state it as a "closed-like" limitation (Feasibility R2). |
| 4 | G9 fails | Report the film as a bounded approximation with its validity map and the fitted mouth-loss K(θ) correction (Inventor R2 §3). That is still a methods result. |
| Ghost | G1 or G3 fail by Nov 14 | One paragraph in the thesis: the leak-path decomposition (R2-C1) and why the ghost was excluded from the gap and closure. |
| Global | Behind schedule by > 2 engineer-days on Nov 18 | Drop Phase 4 and Phase 5. Thesis = O4 + N2-lite + N1 (+ O6). That is still a complete, body-fitted real-gap result. |

---

## 4. Validation statistics

**Principle.** Ghost vs body-fitted and film vs 3D are model-vs-model comparisons of deterministic numbers. The "scatter" across angles is angle-dependent model error, not sampling noise.

The primary rule is deterministic and per point. Compute E_i = S_i − D_i (as % of D_i) at each case and compare it with:
- the fixed margin, and
- u_val = sqrt(u_S² + u_D²), where u = GCI/1.25 (a V&V-20-style procedure; needs literature check).

If |E_i| < 2·u_val, model error cannot be separated from numerical error, and the thesis must say so.

Statistics are secondary and are reported because the committee expects them.

### Families, n and alpha

| Family | Comparison | n | Test(s) | Margin | alpha |
|---|---|---|---|---|---|
| F1 (optional ghost track) | Ghost vs body-fitted: Tz (A/B averaged per angle), FyA | **n = 4 angles** (70/75/80/84). A and B are mirror images, so n_eff is 4, not 8 (C7) | TOST on % difference; 90 % CI | Tz ±5 %, Fy ±10 % | 0.05 per TOST; BH q = 0.05 across the family (m = 2, or 3 if mouth Δp_h is added per Challenger R2 §2a) |
| F2 | Box vs full v20 (V1) | n = 2 angles (70, 84) × 4 quantities | **No statistics.** Deterministic: ≤ 3 % pass, > 5 % fail | — | — |
| F3 | Film vs 3D (box or N1): Q_h, area-mean WSS, 99th-percentile WSS | n = 5 cases (70/75/80/84 + closed). Add 40/60 um at closed for n = 7 (not independent, which must be stated) | TOST on % difference, Bland–Altman, regression | Q ±10 %, mean WSS ±15 %, p99 WSS ±15 % | 0.05 per TOST; BH q = 0.05 over m = 3 |
| F4 | N1 vs analytic slot and literature ranges | — | Range checks only (±30 % analytic; 0.3–3x of the 612 Pa Simon peak; leakage band) | — | — |

### Bland–Altman
- Compute on % differences (equivalently log-ratios, since the error is proportional and torque varies threefold across angles, C6).
- Report the bias with its 95 % CI and both LoA with their 95 % CIs, using SE(LoA) ≈ sqrt(3 s²/n).
- At n = 4 each LoA is uncertain by ±2.8 % (sd 1 %) or ±5.5 % (sd 2 %) (C7, [J4]).
- So the LoA are **descriptive**, not a gate.

### Regression
- Report slope and intercept with 95 % CIs, for ghost vs body-fitted and film vs 3D.
- **R² is not a gate.** It measures how well the comparison tracks the threefold torque range across angles, not agreement. Any proportional bias still gives R² ≈ 1: 0.9998 under −3 %, −4.3 % and ~0 % bias (R2-C6, [J3]), and 1.000 under a pure −8.6 % bias (C7).
- **The slope CI is not a gate either.** At n = 4 its half-width is ~0.04–0.05, so a [0.95, 1.05] CI criterion fails even for an unbiased model [J3].
- Use the slope point estimate only as a diagnostic of proportional vs additive bias.

### Paired tests
- A paired t-test answers "is bias ≠ 0?". With a deterministic model error the answer is eventually always yes. So it is not a gate (C7).
- Use **TOST**: two one-sided t-tests on the mean % difference, df = n − 1. Equivalence means the 90 % CI lies inside the margin.
- If the committee asks, also report the paired-t p-value of the % difference, included in the same BH family.

### Benjamini–Hochberg
- Apply it within each family only (F1, F3), at q = 0.05, ranking the TOST p-values.
- Critical values are i·0.05/m: for m = 3 they are 0.0167 / 0.033 / 0.05 (R2-C6).
- With m ≤ 3 it changes little. State that openly rather than presenting it as rigour.

### Honest power at n = 4 ([J2], 20,000 simulations, margin ±5 %)
- **sd = 1 %:** power 1.00 / 0.91 / 0.74 / 0.46 at true |bias| 2 / 3 / 3.5 / 4 %. The largest observed bias that can pass is ≈ 3.8 %.
- **sd = 2 %:** power 0.74 at 2 % and 0.46 at 3 %. Pass needs |bias| ≲ 2.65 %.
- A model with a true −4.3 % bias, inside 5 % at 3 of 4 angles, is **not shown** equivalent (p = 0.080, R2-C6).
- Conclusion: n = 4 can show equivalence only for clearly good agreement (≲ 3 %). A "not shown" result is inconclusive, not a failure, and the primary deterministic rule decides.

---

## 5. Top 3 risks and early detection

1. **The body-fitted "truth" does not resolve the gap, or does not converge** (points 0.1 and G0b; C11).
   - Detection, **Nov 10–11:** count the elements across the 50 um gap in v20_fixedangles (mesh statistics on a cut plane through the ear end gap and a lateral gap). Require ≥ 4.
   - Then V1 by **Nov 16** and G5b by **Nov 18**.
   - If V1 fails at > 5 % while the box passes G5b, the full v20 model was the under-resolved one.
2. **N1 geometry is not buildable**: closed angle unknown, ear–pocket interference over the ~55° swing, seat contact.
   - Detection: advisor email **today**.
   - The O6 swing-clearance map by **Nov 12** must show min clearance > 0 at every 1° from 85° down to θ_closed.
   - A successful closed-angle geometry plus mesh build (no zero-thickness gaps, G0b ≥ 3) by **Nov 19**. Otherwise use the Phase 3 fallback on Nov 20.
3. **Engineer-time overrun.** About 13 days are needed against about 12.5 available (Feasibility R2 §5, zero margin). Silent overnight failures add to it (the cuDSS LinRes = 1 precedent).
   - Detection: every morning, check G0 (LinRes) on every overnight run before using it.
   - Milestone check **Mon Nov 16**: V0/V1 done.
   - Milestone check **Wed Nov 18**: G5/G5b done.
   - Slip of > 2 days on Nov 18 triggers the Global fallback (drop Phases 4–5).
   - The ghost track is hard-stopped on Nov 14 so it cannot consume this margin.

Watch list (not top 3):
- Thin-Film ↔ 3D coupling in 6.4 is UNVERIFIED (Feasibility R1 §5). Fallback: Coefficient Form Boundary PDE. Test the interface on a flat slot against the analytic flux by Nov 20.
- Novelty erosion by prior art (§6).

---

## 6. Novelty and publishability (stated conservatively)

**Prior art that bounds every claim (all abstract-level only, all need a full-text literature check):**
- Simon, Ge, Sotiropoulos & Yoganathan, Ann Biomed Eng 2010 (doi 10.1007/s10439-009-9857-0). 3D SJM hinge flow under aortic conditions, immersed-boundary method with measured leaflet motion. Reported peak shear about 192 Pa in systole and 612 Pa in diastole.
- Simon et al., Ann Biomed Eng 2010, 38(11):3295–3310 (doi 10.1007/s10439-010-0086-3). Three hinge designs; gap width matters, recess shape matters more.
- Jun et al. 2014, micro-PIV at three gap widths (PubMed 24976188). Trade-off between shear and washout.
- Dumont et al. 2007. Stress-accumulation fractions.

**Consequence:** "two-scale or prescribed-motion hinge CFD" and "diastolic hinge leakage CFD" are **not novel** in themselves (Inventor R1/R2, Challenger R2 §3(i)). The gap value used by Simon (~150 um per the Challenger and Clinical) is UNVERIFIED.

**Likely publishable** (a methods/verification paper, CMBBE or J Biomech Eng tier, per Challenger R2 §3(i); conditional on the gates passing):
1. **Body-fitted, mesh-verified (local GCI) hinge flow at the real CAD 50 um clearance.** This covers open angles and the closed-angle diastolic leakage, with gap WSS, leakage Q, pocket age-of-fluid and pathline stress accumulation, plus 40/50/60 um and Newtonian/Carreau sensitivity. The novelty is the verification at 50 um, not the topic.
2. **A 3D → Reynolds film → 0D model hierarchy with the error measured at each level and angle** (Inventor R2 §3, Advocate R2 §5). The deliverables are a predictive R(θ) = 12μ∫ds/(w h³) from the CAD clearance map and a fitted mouth-loss K(θ). The team found no published lubrication or lumped hinge law for MHVs, but that absence is **unconfirmed and needs a literature check before any "first" claim**.
3. **A quantitative justification of one-way hinge–blade coupling** (friction torque 1e-4 to 2e-3 of blade torque; hinge leakage ~0.5 % of stroke volume; R2-C3). It belongs as a section of paper 1 or 2, not as a standalone paper.

**Thesis-only (useful, but not a paper on its own):**
- The ALE remeshing/kinematics demonstration and the quantitative explanation of the end-face inversion (C9). Add the O6 rotation+slide finding if one emerges.
- The ghost-leaflet results: the leak-path decomposition near a thin gap (R2-C1: slip 3.6x, Darcy core 12KT/h³, transleaflet leak) and the K sensitivity. A short note is possible, but penalization error analyses already exist (Angot et al. 1999; Hester et al.; a 2024 JFM Brinkman force paper whose authors are unidentified). Needs a literature check.
- N3 2D closure and the frozen-window error test, if done.
- The analytic result that stress accumulation per pass in a slot, τ·t = 6μL/h, does not depend on Δp (R2-C4). It is a neat observation but probably known lubrication algebra; needs a literature check.

**Do not claim:**
- closure dynamics or impact;
- cavitation;
- two-way hinge FSI;
- leaflet-surface WSS from the ghost model;
- any "first" claim before the literature check.

The Clinical counterpoint (Hedayat & Borazjani APS abstract: bulk activation exceeds hinge activation) should be acknowledged.

---

## 7. Corrections to the handoff found by the debate

1. **A/B symmetry is 2.7 % at 84° (1.04 / 0.62 / 1.80 / 2.74 %), not "within 2.2 %".** With 2-significant-figure data the 84° value can only be bounded to 1.4–4.1 % (C6, [J1]). Report 4 significant figures.
2. **ALE rate:**
   - 0.026 sim-s/h on average (0.25–0.75°/h; 0.36°/h overall), at 65.4 % completion. Projected finish Nov 10, 08:40–17:30 (C10, Feasibility R1 §0).
   - 70 + 14 = 84°: the run never reaches the 85° stop (C13).
   - The prescribed 14°/s is about 18x slower than physiological opening, and flow is OFF. It is a kinematics result only (Clinical R1 §2).
3. **The 10.8 mm radius is misused.** 10.8 mm is the z-coordinate of the ear end faces (span direction), not the in-plane pivot-to-tip radius. "~2.6 mm tip displacement at 14°" is not a valid pass/fail check until the tip radius is read from the CAD (C9b).
4. **The halo invalidates every K = 1e-9 hinge-region result, and the −8.6 % / +25.6 % numbers are not validation.**
   - The +44 % volume corresponds to a 150–240 um halo per face, i.e. 3–5 gap widths, which seals the hinge (C5, R2-C7).
   - Phase A currently fails deterministically: −8.69 % (R2-C6a).
   - The opposite-signed Tz/Fy errors are unexplained, and a traction-based check is needed (Challenger R1 §2.3).
5. **The R² > 0.95 gate cannot discriminate.** It passes under any proportional bias (C7, R2-C6). The paired t-test asks the wrong question. The judge adds that a slope-CI-in-[0.95, 1.05] gate fails even for an unbiased model at n = 4 [J3]. Both are replaced (§2, §4).
6. **"~270 Pa at closure" is mislabelled.** It is the Poiseuille term under the closed-valve diastolic load (80–100 mmHg, L ≈ 1–1.5 mm). The closing motion itself contributes only 4–20 Pa Couette (C8). Similarly, "eps·Re ~ 0.24 at closure" should be read as 0.12–0.60 for diastolic leakage, depending on L (C8).
7. **The ghost leaflet cannot hold closure or diastolic pressure at any runnable K.**
   - Darcy leak of 4.4–5.4 m/s at K = 1e-9.
   - About 2 L/min through both leaflets even at K = 1e-11 (C2).
   - Next to a gap it multiplies the slot flux by 1.2–3.9x (slip), or about 30x including the Darcy core (C3, R2-C1).

   The plan's Phases D/E through closure were therefore invalid.
8. **"K study" on a fixed 0.61 mm mesh is a sensitivity, not a convergence study.**
   - Three K values give no Richardson redundancy (C4).
   - K = 1e-10 costs ~19 h for the full window (10.9–29.2 h), not 5.7 h. Use 0.025 s restart windows: 1.2 h at 3e-10, 3.7 h at 1e-10 (R2-C2).
9. **A 3-mesh refinement GCI does not fit in 28–32 GB.** r = 1.3 refinement means 3.5 M DOF and 42–93 GB. The current mesh must be the finest grid, and the gap GCI belongs in a local box (C11). The coarse grids must keep ≥ 3 elements across the gap [J5].
10. **"Mass balance < 1 %" and "stabilization < 2 %" (Phases D/E) had no design behind them** (Challenger R1 §2.6). Both are retired.
11. **The 70–85° study window excludes diastole and closure, the clinically dominant hinge phase** (Clinical R1 §1, §3). The thesis needs N1.
12. **Minor:**
    - TzB difference is −8.74 % (handoff −8.8 %, rounding).
    - The sign convention in the torque expression (x + pin_x)/(y − pin_y) is UNVERIFIED for leaflet B (C13).
    - The closed-stop angle, the waveform source and the physiological validity of the 50 um clearance are all UNVERIFIED and must come from the advisor or CAD.

---

## 8. Single next action today (Tue Nov 10)

**Email Dr. Bellofiore today**, cc Ayah. Ask for:
1. the closed-stop leaflet angle and the seat clearance in the CAD, or approval to impose a minimum seat gap;
2. whether 50 um is the physiological/manufactured hinge clearance and its tolerance (to justify the 40/50/60 um sensitivity);
3. aortic vs mitral position and the source of the Study 2 waveform, plus the diastolic Δp to use (80–100 mmHg);
4. any PIV or experimental hinge data available for validation.

This is the longest-lead dependency on the critical path. N1, the clinically dominant deliverable, cannot start geometry without item 1. Everything else today (ALE log, LinRes, G0b element count, 4-significant-figure re-export, FixChi overnight) is under the student's own control.

---

## 10-line summary
1. Core plan: O4 two-scale = body-fitted real-hinge box (N2-lite, 70/75/80/84°) + closed-angle diastolic leakage (N1), with the O6 pivot fit first.
2. Novelty layer: Reynolds film (O1-R) validated against the box/N1, giving a 3D -> film -> 0D hinge-law hierarchy; ghost-film coupling and Phases D/E are dropped.
3. Ranking: O4 80% > N2-lite 70% > O6 85% > N1 55% > O1-R 45% > O2 35% > N3 50% > N5 55% > O3 appendix > N4 10% > original O1 <10% > O5 <5%.
4. Ghost leaflet: never inside the gap or at closure (slip 3.6x, ~30x with Darcy core, 4.4 m/s transleaflet leak); systolic cross-check only, hard stop Nov 14.
5. New gate G0b today: >= 4 elements across the 50 um gap in v20, or the "body-fitted reference" is not a gap reference.
6. Gates: V1 box-vs-full <= 3% (Nov 16); bulk and gap GCI < 5% (Nov 18); N1 WSS within ±30% of slot analytic (Nov 24); film G9 10%/15% (Nov 25); freeze Nov 27.
7. Statistics: primary rule is per-angle |E| against margin and 2*u_val; TOST + BH secondary; R2 and slope-CI are not gates (both non-discriminating at n = 4).
8. Honest power: n = 4, sd 1% -> TOST shows |bias| < 5% only if bias is <~3.8% (power 0.46 at 4% true bias).
9. Handoff fixes: A/B 2.7% (not 2.2%); ALE 0.026 s/h and ends at 84°; 10.8 mm is the z of the end face, not the tip radius; halo voids hinge results; "270 Pa" is diastolic Poiseuille.
10. Next action today: email Dr. Bellofiore for closed-stop angle, seat gap, 50 um tolerance, waveform/valve position, and PIV data.
