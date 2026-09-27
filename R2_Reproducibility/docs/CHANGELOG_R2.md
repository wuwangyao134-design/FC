# R2 change log

- Replaced the earlier locally bounded exact comparison with globally
  certified 2x2 and 3x2 micro-instance validation.
- Added the adapted MODDPG implementation, fixed configuration, saved
  training diagnostics, and a separated pre-correction audit record.
- Corrected the dB/linear SNR conversion and bandwidth-dependent receiver
  noise calculation in the evaluation model.
- Recomputed S1--S8 with 30 independent runs and consistent normalized metric
  definitions, including explicit feasible-front accounting.
- Added the S8 empty-front handling required for valid post-processing.
- Repeated the formal S5 ablation experiment with matched environments and
  random seeds.
- Distinguished task completion time from algorithm wall-clock runtime in
  both saved statistics and documentation.
- Replaced machine-specific active paths in the released analysis scripts
  with paths relative to this package.
- Removed legacy root-level duplicates and paper-specific plotting utilities
  from the public reproducibility package.
- Withheld the legacy ablation MAT archive pending resolution of its
  HMD-Full/w/o-HO field-label inconsistency; the reproducible ablation driver
  remains included.
