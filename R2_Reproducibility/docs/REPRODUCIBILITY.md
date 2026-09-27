# Reproducibility workflow

## 1. Environment

- MATLAB R2024b or later;
- Statistics and Machine Learning Toolbox;
- Deep Learning Toolbox for MODDPG;
- BARON plus the MATLAB--BARON interface only when rerunning the certified
  micro-instance calculations.

BARON is third-party software and is not included in this repository.

## 2. Restore the formal data

Download
[`FC_R2_Results_S1-S8.zip`](https://github.com/wuwangyao134-design/FC/releases/latest/download/FC_R2_Results_S1-S8.zip)
from the latest release. Extract its eight MAT files directly into
`data/formal_results/` and compare their SHA-256 values with
`data/manifest.csv`.

## 3. Reproduce the S1--S8 statistics

In MATLAB:

```matlab
cd('R2_Reproducibility/Compare_nsga_slot')
Summarize_S1_S8_Metrics
```

Generated CSV, MAT, and LaTeX summaries are written to
`R2_Reproducibility/results/generated/S1_S8_Metric_Summary/`.

A successful run prints the resolved S1--S8 source paths before reporting the
per-scenario metrics. This saved-result workflow is the recommended first
verification because it does not rerun the optimizers.

IGD and HV are summarized only over runs that produce a valid feasible
non-dominated set. FFR is calculated over all 30 runs. Failed runs contribute
zero to the number of non-dominated solutions.

## 4. Repeat a benchmark scenario

Open `Compare_nsga_slot/M4.m`, set `scenario_id` to an integer from 1 to 8,
and run the script from that directory. The output is written beneath a new
timestamped `ExperimentResults_Output` directory. Each scenario is run
separately so that failures or interruptions do not invalidate other cases.
S1 is recommended for an initial end-to-end execution check. The formal setup
uses 30 independent runs, 10 slots, 100 individuals, 200 generations, and
20,000 MODDPG training episodes, so large scenarios require substantial
computation.

## 5. Repeat the ablation experiment

```matlab
cd('R2_Reproducibility/S_3test')
Ablation_main
```

The program applies matched environmental realizations and random seeds to
the complete method, the three single-module ablations, and standard NSGA-II.

## 6. Repeat the certified validation

After installing BARON and adding its MATLAB interface to the MATLAB path:

```matlab
cd('R2_Reproducibility/Compare_nsga_slot/MicroExactValidation')
Run_Micro_Exact_Validation('formal')

cd('../MicroExactValidation_3x2')
Run_Micro_Exact_Validation_3x2('formal')
```

Use `quick` only as an installation smoke test. Values reported in the paper
come from `formal` mode. The included result records allow the reported gaps
and solver bounds to be inspected without redistributing BARON.
