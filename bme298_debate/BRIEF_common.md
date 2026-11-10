# Common brief for all debate agents (BME 298 hinge-gap FSI)

Read first: /home/user/India-app/bme298_debate/HANDOFF_BME298_2026-11-10.md (the thesis status).

GOAL being debated: which method(s) will most reliably produce an ACCURATE, NOVEL, PUBLISHABLE result on the REAL
50 um hinge gap of a bileaflet mechanical heart valve within ~4 weeks (defense ~Dec 4-8, 2026), on a 4-core desktop
with 28-32 GB RAM (HPC access expected soon, not guaranteed), COMSOL 6.4. Precision and novelty matter more than
simplicity, but the thesis must NOT depend on an unproven method. Advisor rule: numerical settings may change freely;
model values only if physiologically correct.

OPTIONS THAT MUST BE EVALUATED (by every agent, at least briefly):
 O1 Hybrid ghost leaflet (volume penalization, Brinkman) + Reynolds thin-film hinge (current choice)
 O2 Ghost leaflet alone (K study / local refinement)
 O3 Body-fitted ALE with automatic remeshing (reached ~9 deg+, inverted elements at ear end faces z=+/-10.8 mm)
 O4 Two-scale approach (blade FSI with simplified hinge + separate real-hinge fixed/prescribed-motion model)
 O5 Overset mesh built manually in COMSOL
 O6 Pivot optimization / rotation+sliding kinematics
 O7+ any new option proposed by the team

RULES
- Every quantitative claim must be computed (Python 3 with numpy/scipy is installed; run with `python3`) or derived
  step by step. Mark anything you did not verify as UNVERIFIED.
- Literature: you may try web search (load the WebSearch tool via ToolSearch "select:WebSearch" if you have it).
  Cite ONLY papers you actually found and checked in this session (give title, authors, year, venue, and URL/DOI).
  Otherwise write "needs literature check". NEVER invent citations.
- Be honest: if something will probably not work, say so and why.
- Write your contribution as Markdown to the file path given in your task. Keep it focused: ~800-1500 words.
  Use headings. End with a short "Bottom line" section with your option ranking (O1..O7) and one-sentence reasons.
- Do not edit other agents' files. Do not touch git.
