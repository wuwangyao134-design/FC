BARON micro-instance exact validation
=====================================

Files
-----
Run_Micro_Exact_Validation.m
baron_unpruned_micro_solver.m

Place this folder directly inside Compare_nsga_slot so the main program can
find EvaluateParticle.m, HMD_NSGA_II.m, and REMODDPG_Baseline.m.

First run the smoke test from MATLAB:

    cd('R2_Reproducibility/Compare_nsga_slot/MicroExactValidation')
    validation = Run_Micro_Exact_Validation('quick');

Only after the quick run completes with certified lower/upper bounds and no
physics mismatch, run the formal experiment:

    validation = Run_Micro_Exact_Validation('formal');

The formal mode uses the paper's HMD-NSGA-II population/generation settings,
20,000 MODDPG training episodes, 30 independent runs, and three canonical
certified preferences (energy-only, balanced, and delay-only).
The quick mode is only a functional smoke test and must not be reported.

The certified instance contains I=2 MTDs and M=2 FNs in a 3 m x 3 m
industrial work cell.  Its compact region makes d>=5 m provably impossible,
so the corresponding empty path-loss branch is removed by exact bound
propagation rather than heuristic pruning.

Both FNs use the same 5 GHz-equivalent cycle rate in this correctness-only
instance.  With M=I=2 this gives a formal dominance result: assigning one MTD
to each FN, placing each FN at its MTD (the model clips d to 1 m), and removing
queueing weakly improves both objectives for any original feasible solution.
BARON therefore globally solves the remaining two continuous bandwidth
variables without empirical pruning.
