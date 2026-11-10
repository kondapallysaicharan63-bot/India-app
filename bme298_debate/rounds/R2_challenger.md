# R2 — Challenger rebuttal

Labels follow OPTION_MAP_R2. New numbers come from
`/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/challenger/r2.py`.
Blood: mu = 3.5 mPa s; gap h = 50 um.

## 1. Reconciling 29.7x (mine) with 3.9x (Numerical, C3)

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

## 2. Responses to the Advocate

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

## 3. Attacking my own alternative (N2 = O4 + O6 + frozen-angle body-fitted hinge)

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

## 4. Clinical (N1, N3) and Feasibility

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

## Bottom line (updated ranking)

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

## The one decisive experiment that would change my mind

Run the **FixChi sharp-χ ghost at 70° with local refinement to h_mesh ≤ sqrt(K) around one hinge, at K = 1e-10**.
Compare the **mouth-to-mouth hinge Δp** and FyA with the body-fitted v20 70° run on two meshes. If both agree
within 5%, with the equivalence CI inside ±5%, the ghost can feed the film honestly. I would then rank O1 Phases
A–D above N2, because it covers motion without hundreds of remeshes. If the Δp is off by more than 10%, O1's hinge
WSS is unusable and N2 stands. The run costs about 6 h of compute (Numerical C4 estimate for K = 1e-10) plus one
day of mesh work.
