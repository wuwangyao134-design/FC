# HMD-NSGA-II Reproducibility Repository

This repository contains the MATLAB implementation and reproducibility
materials for the revised HMD-NSGA-II study on multi-objective resource
orchestration in fog-enabled industrial IoT networks.

The reviewed release is organized under
[`R2_Reproducibility`](R2_Reproducibility/README.md). The original files at
the repository root are retained temporarily for backward compatibility.

## What is included

- HMD-NSGA-II and five comparison algorithms;
- the adapted MODDPG baseline and its saved training diagnostics;
- S1--S8 benchmark and statistical-analysis scripts;
- the formal S5 ablation experiment;
- globally certified 2x2 and 3x2 micro-instance validation code and results;
- checksums and documentation for the formal 30-run result archive.

## Quick start

1. Install MATLAB R2024b or later. The evolutionary experiments use the
   Statistics and Machine Learning Toolbox; MODDPG additionally uses the Deep
   Learning Toolbox.
2. Download the release asset `FC_R2_Results_S1-S8.zip` and extract
   `S1_Formal.mat` through `S8_Formal.mat` into
   `R2_Reproducibility/data/formal_results/`.
3. In MATLAB, change directory to
   `R2_Reproducibility/Compare_nsga_slot`.
4. Run `Summarize_S1_S8_Metrics` to reproduce the S1--S8 statistical table,
   or edit `scenario_id` in `M4.m` and run it to repeat a benchmark scenario.

The BARON executable and MATLAB--BARON interface are not redistributed.
They are required only for rerunning the certified micro-instance studies;
the corresponding formal outputs are included in the repository.

See [REPRODUCIBILITY.md](R2_Reproducibility/docs/REPRODUCIBILITY.md) for the
complete workflow and [DATA_DICTIONARY.md](R2_Reproducibility/docs/DATA_DICTIONARY.md)
for the saved-variable definitions.

The prepared release-asset names and archive checksums are listed in
[RELEASE_ASSETS.md](R2_Reproducibility/docs/RELEASE_ASSETS.md).
