# EEG Meditation Connectivity

Dynamic EEG functional connectivity analysis of meditation and mind-wandering,
using phase-locking value (PLV) estimated from a time–frequency representation.

> **Volunteer research.** This is an independent, volunteer research project.
> It is not peer reviewed, it is not affiliated with the dataset authors, and it
> is not intended for clinical use. Results are exploratory.

## Overview

The project studies how large-scale functional connectivity between EEG
channels changes over time during meditation, and whether connectivity patterns
differ between focused meditation and mind-wandering episodes.

**Data:** [OpenNeuro ds001787](https://openneuro.org/datasets/ds001787), an EEG
meditation study (64-channel BioSemi recordings, `.bdf`) in which participants
meditated and were periodically probed to report their mental state. The raw
data is **not** included in this repository; download it from OpenNeuro.

## Pipeline

| # | Stage | What it does |
|---|-------|--------------|
| 1 | **Preprocessing** | Import BioSemi `.bdf` files, assign channel locations, re-reference, filter, and resample. |
| 2 | **Artifact removal** | Remove bad channels and high-amplitude bursts with EEGLAB `clean_rawdata` (ASR), then run ICA and reject artifactual components (eye, muscle, line noise). |
| 3 | **Sliding-window PLV** | Compute instantaneous phase from the **Rihaczek time–frequency distribution** and estimate PLV for every channel pair in overlapping sliding windows, giving a dynamic connectivity time series. |
| 4 | **Network measures** | Derive **node strength** (summed connectivity per channel), **connection strength** (per channel pair), and **hemisphere maps** (within- and between-hemisphere connectivity) and plot them on the scalp. |
| 5 | **State split** | Use the experience-sampling probe responses to label the epochs before each probe as **meditation** or **mind-wandering**, then compare connectivity between the two states. |

```
.bdf raw EEG ─► preprocessing ─► clean_rawdata + ICA ─► Rihaczek TFD ─► sliding-window PLV
                                                                              │
                     meditation vs mind-wandering ◄── node / connection / hemisphere maps
```

## Tools

- **MATLAB** (Signal Processing Toolbox recommended)
- **[EEGLAB](https://sccn.ucsd.edu/eeglab/)** with the `clean_rawdata` and
  `BIOSIG` (for `.bdf` import) plugins
- ICA through EEGLAB (`runica`, or ICLabel for component classification)

## Repository layout

```
eeg-meditation-connectivity/
├── src/          MATLAB scripts (.m) for every pipeline stage
├── README.md
├── LICENSE
└── .gitignore    excludes raw data, .mat/.set files, and results folders
```

Data and outputs (`data/`, `results/`, `results_clean/`, `results_state/`) are
generated locally and are intentionally not version-controlled.

## How to run

1. Install MATLAB and EEGLAB, then add the `clean_rawdata` and `BIOSIG` plugins
   through *File → Manage EEGLAB extensions*.
2. Download ds001787 from OpenNeuro and put the subject folders in `data/` at
   the repository root (or edit the data path at the top of the scripts).
3. Start MATLAB in the repository root and add the code to the path:
   ```matlab
   eeglab nogui;            % initialise EEGLAB
   addpath(genpath('src'));
   ```
4. Run the scripts in pipeline order: preprocessing → artifact removal →
   PLV → network measures → state split. Outputs are written to `results*/`.

## License

Released under the [MIT License](LICENSE).
