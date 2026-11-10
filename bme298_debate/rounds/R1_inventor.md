# R1 — Inventor: two new options (O7, O8) for the real 50 um hinge

All numbers below come from `/tmp/claude-0/-home-user-India-app/3e7bbeff-659d-5fae-9fc6-69d16cec2281/scratchpad/inventor/calc.py` (plus one inline check).
Inputs: rho = 1060 kg/m3, mu = 3.5 mPa s, h = 50 um, ear radius from pivot R = 2.2 mm (where the inverted elements appear).
**UNVERIFIED (need CAD/model):** gap path length L_s ~ 1 mm, wetted gap width w ~ 6 mm, end-face disc radius 1.5 mm, closing peak Omega ~ 100 rad/s, opening and closing times of 40 and 30 ms.

## 0. Scale facts that shape the rest

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

## O7 (new): frozen-geometry moving-wall (FGMW) real-hinge submodel + calibrated hinge impedance

### Math
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

### COMSOL sketch (6.4)
1. Build a hinge box submodel (about 6 x 6 x 4 mm around one pocket) from v20.
2. Use a parameter th_k in a Rotate geometry node and remesh per angle with a boundary-layer mesh in the gap. This is the same workflow as v20_fixedangles.
3. Set the leaflet wall to Moving wall with u = -thd*(y - yp), v = thd*(x - xp). theta_dot(t) comes from the blade FSI or the prescribed kinematics.
4. Take box inlet/outlet conditions from the blade-FSI solution via Interpolation functions exported on the box faces.
5. Hand over each window's solution as the next window's initial value. Export the field on a grid, then import it as an interpolation function; this is safer than cross-mesh `initsol` after a geometry rebuild (UNVERIFIED which of the two works in 6.4).
6. Run windows only for opening and closing: about 30 windows each at 0.5 deg. Use plain fixed geometry elsewhere.

There is no ALE, no remeshing, no inverted elements, and no penalization parameter K.

### Novelty and risk
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

## O8 (new): split-gap sliding interface, a rescue for O3 (idea b, reformulated)

### Why not a true rotating frame
In a frame that co-rotates with the leaflet, the NS equations gain rho*(2 Omega x u + Omega x (Omega x r) + Omega_dot x r).
- Coriolis/viscous is 2*rho*Omega*l^2/mu: 0.01 (gap, closing mean) to 0.15 (gap, peak), but **1.3 to 15 at pocket scale (0.5 mm)**. The terms cannot be neglected.
- The butterfly pocket walls are not axisymmetric, so they would move in that frame and still need ALE.
- Rejected. Rigid-rotation ALE (sliding mesh) needs no fictitious forces at all.

### The idea
The inversion happens because one 50 um mesh layer must absorb 576 um of in-plane slide. Split each end gap at its mid-plane, z = +/-10.825 mm (geometry to be confirmed), with a planar sliding identity pair:
- The housing-side half-layer (25 um) is a separate geometry object with a fixed mesh.
- On the leaflet side, the pair face has prescribed in-plane displacement equal to rigid rotation over the ear footprint and free in-plane sliding elsewhere.
- Laminar Flow "Flow Continuity" on the pair carries velocity and traction across the non-conforming meshes. This is the same mechanism COMSOL uses for rotating machinery. COMSOL docs found: Rotating Machinery, Laminar Flow (https://doc.comsol.com/6.4/doc/com.comsol.help.cfd/cfd_ug_fluidflow_single.06.055.html) and the blog "Deformed Mesh Interfaces: Rotations and Linear Translations" (https://www.comsol.com/blogs/deformed-mesh-interfaces-rotations-and-linear-translations/). Both note the identity pair, Form Assembly, and "destination side finer".

A plane normal to the rotation axis is a surface of revolution, so the slide is pure tangential motion with **zero distortion: the end-gap shear drops from 11.5 to 0**. The remaining distortion moves to the lateral gaps. There, a non-concentric arc gives a normal strain of e*theta/h = 0.05 / 0.26 / 1.05 at 15 deg for e = 10 / 50 / 200 um. Yeoh smoothing handles this for e <= 50 um but not for e = 200 um. The O6 fit again decides.

Concentricity check before any build: sample the pocket-wall and ear-side points in the CAD, fit a circle about the pivot axis, and report the radial scatter sigma_r. Requirement: sigma_r << h, about 5 um, along the thin-gap arcs. If this passes, the lateral thin arcs can also become cylindrical sliding pairs.

### Risk
- The v20 single blood domain must be re-split into an assembly.
- The 25 um slab needs at least 3 elements across its thickness.
- The non-conforming pair mesh must be finer on the destination side.
- **The swept-footprint problem:** the fixed slab must cover the whole area the ear sweeps.
- **P ~ 40%**, at about 1.5 weeks of geometry and mesh work. It does produce the first fully body-fitted moving real-hinge FSI, which is the most impressive result if it works.

## Idea (c) alone and other discards
- (c) without (a) is just O4's prescribed-motion part, so I merged it into O7, where its derived error bound is the contribution.
- A literal rotating reference frame is rejected (Section O8).
- Full two-way coupling of the hinge impedance into the blade FSI is unnecessary (friction torque ratio 1e-3). It is optional for completeness.

## Bottom line (ranking)
1. **O4 + O7**: highest probability of an accurate real-hinge result. The one-way coupling is now quantitatively justified, and FGMW plus the calibrated impedance law add publishable novelty.
2. **O6**: cheap, and it supplies the eccentricity e that sets the error of both O7 and O8. Do it first.
3. **O1**: most novel, but the film is invalid at mouths and corners, and moving sharp chi gives noisy forces and has no Jacobian. Keep it as a stretch goal, and validate it against O7.
4. **O2**: needed as the bulk-flow baseline and for the K/GCI study, but on its own it does not resolve the 50 um gap (sqrt(K) ~ 32 um).
5. **O8**: a targeted fix for O3's end-face inversion; high payoff, ~40% probability.
6. **O3 as is**: proven kinematics, but the mesh at the ear ends is invalid, so its flow there cannot be used.
7. **O5**: no built-in support in COMSOL and about 2x the DOF; leave it for follow-up work.

### Sources checked this session
- Simon HA, Ge L, Sotiropoulos F, Yoganathan AP. "Simulation of the three-dimensional hinge flow fields of a bileaflet mechanical heart valve under aortic conditions." Ann Biomed Eng 2010;38(3):841-853. https://pmc.ncbi.nlm.nih.gov/articles/PMC2841271 (abstract-level facts only).
- "High-resolution measurements of leakage flow through bileaflet mechanical heart valve hinges" (Klusak & Quinlan; year and venue not verified). https://researchrepository.universityofgalway.ie/entities/publication/4dc93b67-aa9f-4bc9-9f78-d7100adac0fc. It reports that Simon et al. estimated a hinge Re of about 150. That definition differs from my h-based Re of about 10, and the difference needs a literature check.
- The COMSOL 6.4 documentation and blog listed under O8.
- Transpiration BC background: NASA report https://ntrs.nasa.gov/api/citations/19970016811/downloads/19970016811.pdf (only the search summary was read; needs literature check before citing it in the thesis).
