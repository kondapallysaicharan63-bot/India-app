# R2 — Clinical / Physiology Expert: rebuttal

Labels follow OPTION_MAP_R2: my R1 "O7" is now **N1**, my "O8" is now **N3**.
Numbers come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/clinical/r2.py`, plus the R1 `calc.py`.
Literature status is unchanged from R1: I checked abstracts and search snippets only, because the full-text hosts were blocked.

## 1. To the Advocate: is gap WSS from a Reynolds film clinically useful without pocket washout?

**Partly, and only as half of a trade-off.** The best experimental hinge-gap study I found, the SJM micro-PIV study with three gap widths (PubMed 24976188, abstract), reports that the narrowest gap had the **lowest** shear but the **poorest washout**. Its conclusion is that hinge design has to balance the two.

The lubrication scalings show why WSS alone points the wrong way. At fixed Δp:

- Leakage flux Q ∝ h³, wall shear τ_w ∝ h, transit time ∝ 1/h², so stress accumulation per pass (τ·t) ∝ 1/h.
- Going from 50 to 40 µm gives Q ×0.51, τ_w ×0.80 and τ·t per pass **×1.25**.
- Going from 50 to 60 µm gives Q ×1.73, τ_w ×1.20 and τ·t ×0.83.

If you report only film WSS, a narrower hinge looks safer. Exposure time and washout say it is worse. So film WSS has clinical value only when it comes with the film's own transit time (from which it gives an analytic τ·t per pass), and with a pocket residence or washout measure from a 3D body-fitted model. The film can supply the first, cheaply and with high accuracy, in the uniform-gap region. It cannot supply the second.

I accept the Advocate's validity numbers. ε·Re is 0.12 to 0.60 in diastolic leakage, and **in diastole the film does not need the ghost bulk at all**. With the valve closed, the driving Δp is the transvalvular pressure, which can be prescribed from the waveform. That makes a standalone diastolic film model clinically credible for gap WSS and flux, and it avoids the problem in §2.

## 2. To the Challenger and Numerical analyst: ghost leaflets leak at closure

The clinical consequence is decisive. Numerical C2 gives a leak through the ghost leaflet of about 2 L/min at K = 1e-11, a K value that cannot be run. That is about **17 mL over a 0.5 s diastole**. At the K that does run (1e-9) the leak scales to about **1.7 L per diastole**.

Compare these with:

- A true hinge leakage of **0.13 to 0.48 mL/beat** for all 4 hinges (slot band: span 4–6 mm, L = 1–2 mm, Δp = 80–100 mmHg).
- A total in-vitro mechanical-valve regurgitant volume of **5 to 11 mL/beat**. This figure is from a search snippet of a multi-valve in-vitro comparison; the exact citation needs a literature check.

Even at the infeasible K, the false leak is **35 to 130 times the hinge leak** and **larger than the whole-valve regurgitation**. The Challenger's flux-ratio and Numerical's velocity-ratio numbers disagree in magnitude, but both say the same thing clinically: the leak is far above the signal at any runnable K.

So **no diastolic or closure output from a ghost bulk model (O1 bulk, O2) is clinically credible.** That includes leakage jets, regurgitant volume, washout of the pocket in diastole, and mouth pressures fed to a film. Penalized flow dominates exactly the region we care about.

The fix is to restrict ghost use to systole, the opening phase and torque, and to drive every diastolic and closure hinge output from body-fitted geometry (N1/N2) or from a prescribed Δp (film). The thesis text should say this explicitly. A reviewer will catch it otherwise.

## 3. To the Inventor: N2 (frozen-angle body-fitted mesh with transpiration and a lumped law)

- **Leakage phase (about 60% of the cycle): credible.** At the closed angle nothing moves, so transpiration is not needed. N2 reduces to my N1 and is exact apart from mesh error. It delivers gap WSS, leakage Q, jets, pocket residence time and pathline stress accumulation. The prerequisite is the **closed-stop angle**, which is still UNVERIFIED. My estimate is about 30°, inferred from an 85° opening (search-verified) and a 55° traverse (snippet). It must come from the CAD or the advisor.
- **Closure phase: credible for hinge flow, with two clinical caveats.**
  1. Couette (wall-motion) shear is only 1 to 15% of pressure-driven gap shear (4–20 Pa against 133–333 Pa). So closure hinge hemodynamics are governed by the Δp(t) history across the closing valve. N2's credibility at closure therefore depends on where its box boundary conditions come from. If they come from the ghost model, §2 applies and it fails. They must come from the waveform, the v21 pin-joint model, or body-fitted runs.
  2. Closure sweeps 85° → about 30°, but every existing geometry and validation covers 70–85°. N2 needs about 55° of new frozen angles, and the ears must not intersect the pocket along the way. This is O6's job.
- **What N2 does not capture:** leaflet impact, rebound, water-hammer and cavitation. The thesis should state these as out of scope.
- **Clinical implication of negligible back-coupling.** The Inventor's friction-torque ratio of about 1e-3 and leakage of about 0.5% of stroke volume mean the hinge does not change leaflet motion or bulk flow. Hinge thrombogenicity is therefore a **post-processable function of (θ(t), Δp(t), h)**. Gap-width and hinge-shape design studies can then be run on the small box without repeating blade FSI, which is clinically the most useful property in the whole debate.
- **Caveat:** this holds for a clean hinge. The clinical failure mode is thrombus or pannus that *restricts* motion, and then the coupling becomes strong. Mention it as a limitation.
- **The lumped R/K/C law is engineering novelty, not a clinical output.** It is fine if cheap, but it is not part of the clinical minimum.

## 4. To Feasibility: what survives on about 13 engineering days

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

## 5. Minimum clinically defensible deliverable set, with acceptance criteria

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

## Bottom line (updated clinical ranking)

1. **N1 + O4**: the only path to credible dominant-phase (diastolic) hinge outputs. Every item in the minimum set comes from it.
2. **O6**: needed to fix h at the closed angle and to make the closing sweep admissible. Leakage ∝ h³.
3. **N2**: the right path for closure-phase hinge flow, but only with Δp boundary conditions that do not come from the ghost model. It is a stretch.
4. **O1, film part only (diastole-driven standalone)**: high-accuracy gap WSS and transit time. It is clinically incomplete without washout.
5. **N3**: a cheap qualitative closure story.
6. **O2 / O1 ghost bulk**: systolic torque only. Clinically disqualified for diastole and closure because the porous leak is 35–130 times the hinge leak even at an unrunnable K.
7. **O3, N4**: kinematics only. **O5**: infeasible.
