# R2 Reproducibility Package

This directory is the authoritative package for the revised experiments.
The code and records are separated according to their evidential role so that
formal manuscript data cannot be confused with historical audit material.

## Directory map

- `Compare_nsga_slot/`: S1--S8 benchmark, HMD-NSGA-II, comparison algorithms,
  metric calculation, and certified micro-instance programs.
- `S_3test/`: formal S5 ablation experiment.
- `data/formal_results/`: extraction target for the S1--S8 release archive.
- `data/global_validation/`: certified BARON result records used in the paper.
- `audit/MODDPG/`: numerical training diagnostics and the
  corrected-versus-pre-correction audit summary. The pre-correction raw input
  is not used in the manuscript.
- `docs/`: reproducibility procedure, data dictionary, and revision log.

## Main entry points

| Purpose | MATLAB entry point |
|---|---|
| Run one S1--S8 scenario | `Compare_nsga_slot/M4.m` |
| Summarize formal S1--S8 results | `Compare_nsga_slot/Summarize_S1_S8_Metrics.m` |
| Run formal S5 ablation | `S_3test/Ablation_main.m` |
| Run certified 2x2 validation | `Compare_nsga_slot/MicroExactValidation/Run_Micro_Exact_Validation.m` |
| Run certified 3x2 validation | `Compare_nsga_slot/MicroExactValidation_3x2/Run_Micro_Exact_Validation_3x2.m` |

The formal benchmark uses 30 independent runs, 10 scheduling slots, a
population size of 100, and 200 generations. MODDPG is trained for 20,000
episodes under the fixed settings recorded in `M4.m`.

The ablation program is included, but its archived MAT file is intentionally
withheld from this release pending resolution of a legacy HMD-Full/w/o-HO
field-label inconsistency. It is not part of the formal S1--S8 archive.

## Data convention

`S1_Formal.mat` through `S8_Formal.mat` are the only S1--S8 files designated
as formal revised-manuscript results. Their SHA-256 checksums are recorded in
`data/manifest.csv`. Files in the MODDPG audit directory are supporting
implementation records and must not be pooled with the formal benchmark.
