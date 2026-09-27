# HMD-NSGA-II: MATLAB Code and Reproducibility Package

This repository contains the MATLAB implementation and reproducibility
materials for HMD-NSGA-II, a multi-objective resource-orchestration framework
for fog-enabled industrial IoT networks. The package accompanies the revised
study and separates formal manuscript results, executable source code, global
validation records, and historical implementation-audit evidence.

The authoritative revision package is located in
[`R2_Reproducibility/`](R2_Reproducibility/README.md).

## Reproducibility at a glance

Readers are encouraged to begin with the saved-result verification before
launching the computationally intensive experiments.

| Level | Purpose | Data or software required | Entry point |
|---|---|---|---|
| Quick verification | Reproduce the S1--S8 numerical summary without rerunning the optimizers | Formal S1--S8 Release asset | `Summarize_S1_S8_Metrics.m` |
| Single-scenario rerun | Verify the complete optimization workflow from a fresh initialization | MATLAB toolboxes listed below | `M4.m` with one `scenario_id` |
| Ablation rerun | Evaluate the complete method and three single-module variants | MATLAB | `Ablation_main.m` |
| Certified validation | Recompute the globally certified 2x2 and 3x2 micro-instances | MATLAB, BARON, and the MATLAB--BARON interface | `Run_Micro_Exact_Validation*.m` |

## Repository contents

```text
R2_Reproducibility/
|-- Compare_nsga_slot/       HMD-NSGA-II, baselines, S1--S8 driver, and metrics
|-- S_3test/                 S5 ablation experiment
|-- data/
|   |-- formal_results/      Extraction target for the formal S1--S8 data
|   `-- global_validation/   Included BARON certificates and gap summaries
|-- audit/MODDPG/            Corrected training diagnostics and audit summaries
`-- docs/                    Workflow, data dictionary, release notes, and checksums
```

Publication-specific plotting and typesetting scripts are intentionally not
included. The package reproduces the optimization outputs and reported
statistics; figure formatting is treated as a separate presentation step.

## Requirements

The code was tested with MATLAB R2024b on Windows. The following components are
used:

- Statistics and Machine Learning Toolbox for the evolutionary experiments;
- Deep Learning Toolbox for the adapted MODDPG baseline;
- BARON and the MATLAB--BARON interface only for recomputing the certified
  micro-instance solutions.

BARON is third-party software and is not redistributed. It is not needed to
inspect the included certified results or to reproduce the S1--S8 summary.

## 1. Quick verification of the reported S1--S8 results

Download
[`FC_R2_Results_S1-S8.zip`](https://github.com/wuwangyao134-design/FC/releases/latest/download/FC_R2_Results_S1-S8.zip)
from the latest GitHub Release. Extract the following files directly into
`R2_Reproducibility/data/formal_results/`:

```text
S1_Formal.mat  S2_Formal.mat  S3_Formal.mat  S4_Formal.mat
S5_Formal.mat  S6_Formal.mat  S7_Formal.mat  S8_Formal.mat
```

The archive also contains `manifest.csv`. Its SHA-256 records are duplicated
in [`data/manifest.csv`](R2_Reproducibility/data/manifest.csv) so that the
formal inputs can be verified independently.

From the repository root, run:

```matlab
repo_root = pwd;
cd(fullfile(repo_root, 'R2_Reproducibility', 'Compare_nsga_slot'));
Summarize_S1_S8_Metrics
```

A successful run identifies all eight formal files in the Command Window and
writes the following outputs to
`R2_Reproducibility/results/generated/S1_S8_Metric_Summary/`:

- `S1_S8_Metric_Summary.csv`;
- `S1_S8_Metric_Summary.mat`;
- `S1_S8_LaTeX_Table.tex`.

This is the recommended first check because it validates the archived data,
scenario identities, statistical definitions, and manuscript-table values
without repeating the full optimization campaign.

## 2. Rerun one benchmark scenario

Open [`M4.m`](R2_Reproducibility/Compare_nsga_slot/M4.m), set
`scenario_id` to an integer from 1 to 8, and run the script from its containing
directory:

```matlab
cd(fullfile(repo_root, 'R2_Reproducibility', 'Compare_nsga_slot'));
M4
```

Each formal scenario uses 30 independent runs, 10 scheduling slots, a
population size of 100, and 200 generations. MODDPG is trained for 20,000
episodes using the fixed configuration in `M4.m`. A timestamped directory is
created below `Compare_nsga_slot/ExperimentResults_Output/`.

For an initial execution check, use S1. A complete formal rerun is
computationally expensive, especially for the larger scenarios, and should be
performed one scenario at a time. Within each independent run, shadow fading
is fixed across scheduling slots and independently resampled between runs.

## 3. Rerun the ablation experiment

```matlab
cd(fullfile(repo_root, 'R2_Reproducibility', 'S_3test'));
Ablation_main
```

The S5 experiment evaluates HMD-NSGA-II, the w/o ISMM, w/o LAMS, and w/o HO
variants, and standard NSGA-II under matched network realizations and random
seeds. See the package README for the status of the archived ablation MAT file.

## 4. Inspect or rerun the certified validation

Formal solver outputs and gap summaries are already included under
`R2_Reproducibility/data/global_validation/`. To recompute them after installing
BARON and adding the MATLAB--BARON interface to the MATLAB path, run:

```matlab
cd(fullfile(repo_root, 'R2_Reproducibility', 'Compare_nsga_slot', ...
    'MicroExactValidation'));
Run_Micro_Exact_Validation('formal')

cd(fullfile(repo_root, 'R2_Reproducibility', 'Compare_nsga_slot', ...
    'MicroExactValidation_3x2'));
Run_Micro_Exact_Validation_3x2('formal')
```

The `quick` mode is an installation smoke test only; manuscript values are
obtained from `formal` mode.

## Metric and naming conventions

- `MyNSGA_II` in saved MATLAB structures denotes HMD-NSGA-II.
- `DRL_Baseline` denotes the adapted MODDPG implementation.
- `MOEA_D` denotes MOEA/D.
- IGD and HV are summarized over runs that produce a valid feasible
  non-dominated set.
- FFR is calculated over all 30 runs.
- Failed runs contribute zero to the reported number of feasible
  non-dominated solutions.

The historical MODDPG audit material is distributed separately as
[`FC_R2_MODDPG_Audit.zip`](https://github.com/wuwangyao134-design/FC/releases/latest/download/FC_R2_MODDPG_Audit.zip).
It documents the corrected implementation and retained pre-correction S2
evidence, but it is not a source of manuscript benchmark metrics.

## Further documentation

- [Complete reproducibility workflow](R2_Reproducibility/docs/REPRODUCIBILITY.md)
- [Saved-variable data dictionary](R2_Reproducibility/docs/DATA_DICTIONARY.md)
- [Release assets and SHA-256 checksums](R2_Reproducibility/docs/RELEASE_ASSETS.md)
- [R2 change log](R2_Reproducibility/docs/CHANGELOG_R2.md)

When reporting a reproducibility issue, please include the MATLAB release,
operating system, selected scenario, and the complete MATLAB error message.
