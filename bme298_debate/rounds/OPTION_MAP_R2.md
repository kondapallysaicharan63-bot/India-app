# Canonical option labels for Round 2 (orchestrator)
Round-1 agents used "O7"/"O8" for different ideas. From Round 2 on, use these labels:
- O1 Hybrid ghost leaflet + Reynolds thin-film hinge
- O2 Ghost leaflet alone
- O3 Body-fitted ALE + automatic remeshing
- O4 Two-scale (blade FSI + separate real-hinge fixed/prescribed model)
- O5 Manual overset mesh
- O6 Pivot optimization / rotation+sliding kinematics
- N1 Closed-angle body-fitted diastolic-leakage case (Clinical's "O7")
- N2 Real-hinge box submodel: frozen-angle body-fitted meshes with full rigid-body wall velocity + fitted lumped hinge
       pressure-leakage law (Inventor's "O7"; merges Challenger's "O7" = O4+O6 + body-fitted hinge submodel)
- N3 2D hinge cross-section with prescribed closing motion (Clinical's "O8")
- N4 ALE end-gap mid-plane sliding split / Flow Continuity pair (Inventor's "O8")
- N5 Stationary peak-flow solves for cheap K screening (Feasibility's "O7")

Known numeric conflicts to resolve in Round 2:
- K-study run time: Feasibility (0.13 s window: 6.5 h @3e-10, 19 h @1e-10) vs Numerical (~2.8 h @3e-10, ~5.7 h @1e-10, K^-0.63).
- Penalized-leaflet leakage: Challenger (29.7x Poiseuille flux at K=1e-9, ~1300 mL/s at closure) vs Numerical (3.9x gap flux; 4.4-5.4 m/s through leaflet at 80-100 mmHg).
- Halo thickness: Advocate ~200 um (0.9 mm leaflet assumed) vs Numerical 150-240 um.
- Added inertia ratio: 1.3-6x (Numerical), 1.5-6.3x (Challenger), ~5.3x (Feasibility).
