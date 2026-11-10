# R1 — Clinical / Physiology Expert

Scope: which outputs matter clinically for hinge thrombosis and blood damage, which cardiac phase dominates, what each option O1–O6 can credibly deliver, and whether the model inputs are physiologically correct.
Numbers come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/clinical/calc.py` (Python, order-of-magnitude slot and Couette estimates). They are not CFD.

## 1. What the literature says matters (checked this session through search abstracts and snippets only)

I could not open full texts because pmc.ncbi.nlm.nih.gov and link.springer.com are blocked in this session. Each claim below rests on the abstract or snippet returned by search. Before the thesis cites any of them, someone must open the full text.

- **Simon HA, Ge L, Sotiropoulos F, Yoganathan AP**, "Simulation of the three-dimensional hinge flow fields of a bileaflet mechanical heart valve under aortic conditions", *Ann Biomed Eng* 2010 (DOI 10.1007/s10439-009-9857-0, https://link.springer.com/article/10.1007/s10439-009-9857-0). Micro-CT SJM hinge. The abstract reports two strong diastolic leakage jets and peak shear of **1920 dyn/cm² (192 Pa) in systole versus 6115 dyn/cm² (612 Pa) in diastole**, so diastole is about 3.2 times systole. The year in the DOI prefix (2009 vs 2010) still needs checking.
- **Simon HA et al.**, "Numerical investigation of the performance of three hinge designs of bileaflet mechanical heart valves", *Ann Biomed Eng* 2010 (https://experts.umn.edu/en/publications/numerical-investigation-of-the-performance-of-three-hinge-designs/). Gap width matters, and **recess/leaflet shape matters more**, for flow distribution and thromboembolic potential. Elevated shear is linked to hemolysis and platelet activation; stasis and recirculation are linked to poor washout and thrombus.
- **Jun BH, Saikrishnan N, Yoganathan AP et al.**, micro-PIV of SJM valves with three hinge gap widths, *Ann Biomed Eng* 2014 (PubMed 24976188, https://www.pubmed.gov/24976188/; author list and title still need checking). The widest gap gave the highest Reynolds shear stress (about 1000 N/m²) during **leakage flow**. The **narrowest gap gave the poorest washout and the most stasis**. The authors conclude the design must balance shear against washout.
- **Govindarajan V, Udaykumar HS, Chandran KB**, 2D hinge model of the closing phase (https://pmc.ncbi.nlm.nih.gov/articles/PMC2744151; exact title, year and venue still need checking). High shear combined with **long residence time in leakage flow during closure** can activate platelets.
- **Dumont K, Vierendeels J, Kaminsky R et al.**, ATS vs SJM Regent FSI, *J Biomech Eng* 2007 (https://researchconnect.stonybrook.edu/en/publications/comparison-of-the-hemodynamic-and-thrombogenic-performance-of-two/). During forward flow, **no platelets exceeded the 35 dyn·s/cm² stress-accumulation threshold**. During regurgitation, 0.81% (SJM) vs 0.63% (ATS) exceeded it, and the difference was attributed to hinge design.
- Hedayat & Borazjani, APS DFD 2017 abstract (https://archive.aps.org/dfd/2017/d4/5). Bulk-flow platelet activation was several times larger than hinge activation for a 25 mm SJM. This is an abstract, not peer reviewed, but it is a fair counterpoint to examiners who argue "the hinge is everything".
- The Bluestein stress-accumulation and power-law models, and Yun/Dasi lattice-Boltzmann hinge blood-damage studies, appeared in the search results but I did not check them. They need a literature check.

**Clinical conclusion:** four studies (Simon 2010, Jun 2014, Govindarajan, Dumont 2007) agree that hinge risk is dominated by **diastolic leakage plus the closure transient**, not by peak systole. The mechanisms are high gap shear (activation) and poor washout in the pocket (stasis and accumulation). Systolic hinge flow is mostly forward washout.

## 2. Phase-by-phase magnitudes (my estimates)

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

## 3. Physiologic correctness of the inputs

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

## 4. Option-by-option: what each can credibly deliver

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

## 5. Minimum clinically meaningful deliverable set (publishable)

1. Validated fixed-angle torques at 70–84° with mesh GCI. This already exists, and the GCI is still missing.
2. **A closed-angle (O7) leakage case**: hinge gap WSS (max, area-weighted), leakage flow rate per hinge, jet description, and pocket residence time (tracer washout or mean age-of-fluid) with a Newtonian vs Carreau sensitivity.
3. **Particle-based stress accumulation** (linear τ·t plus one power law) along pathlines through the gap, opened vs closed, reported as a distribution and the fraction above 35 dyn·s/cm². Compare it against Dumont 2007's 0.81%, citing it as a published reference to compare against, not as a match.
4. O1 thin-film gap WSS versus 3D (O4/O7) in the uniform-gap region, as the methods-novelty piece.
5. A gap-width sensitivity (for example 40/50/60 µm, motivated by O6 and Jun 2014). This is cheap and clinically framed as the trade-off between shear and washout.

Out of scope, stated explicitly: impact, cavitation, and two-way FSI through closure.

## Bottom line (clinical ranking)

1. **O4, extended with O7 (closed-angle leakage)**: the only path that credibly delivers diastolic gap WSS, pocket washout and SA in the phase that dominates hinge thrombosis.
2. **O1**: an elegant and novel analytic-quality gap WSS, but it cannot give pocket recirculation, leaflet WSS or closure. It is valuable only when validated against item 1.
3. **O6**: cheap, and it fixes the clinically dominant input (h at each angle, especially closed).
4. **O8 (2D closure sub-model)**: low cost, adds the closure-shear story.
5. **O2**: gives torque, but its wall shear and leakage are K-contaminated, so it is weak clinically.
6. **O5**: right in principle for closure, infeasible in 4 weeks.
7. **O3**: no hemodynamics and an invalid mesh at the ears, so it is a feasibility appendix only.
