# Data dictionary

## Formal S1--S8 MAT files

Each `S*_Formal.mat` file contains one scenario and the following principal
variables.

| Variable | Meaning |
|---|---|
| `all_scenario_results` | Run-by-slot metrics for each evaluated algorithm |
| `current_scenario_all_run_final_fronts` | Raw objective vectors indexed by algorithm, slot, and run |
| `all_reference_fronts` | Run-specific empirical non-dominated reference sets |
| `all_normalized_reference_fronts` | Reference sets in the shared normalized objective space |
| `all_metric_normalization_bounds` | Ideal, worst, range, and HV reference points used for metric computation |
| `all_REMODDP_training_info` | MODDPG reward, critic-loss, feasibility, and timing histories |
| `statistical_results_to_save` | Terminal-slot descriptive statistics |
| `alg_names_for_results` | Saved algorithm-field ordering |
| `scenario_catalog` | S1--S8 system-size and region definitions |
| `scenario_id` | Scenario represented by the file |
| `nSlots` | Number of successive scheduling slots |
| `num_stat_runs` | Number of independent runs |

Within `all_scenario_results`, the main fields are `IGD`, `HV`, `Spacing`,
`Spread`, `Runtime`, `NumFeasibleSolutions`, and `Tmax`. Each field contains a
cell whose numeric matrix is indexed by independent run and scheduling slot.

`Runtime` is the wall-clock time of one algorithm invocation for one slot.
`Tmax` is the maximum task completion time represented by the generated
schedule. These quantities describe different stages of operation and should
not be interpreted interchangeably.

## Certified validation records

The variable `validation` stores the problem configuration, BARON lower and
upper bounds, heuristic fronts, preference-wise gaps, runtimes, and the
summary table. `AllScalarizationsCertified` must be true for a record used in
the paper.

## Ablation output

Running `Ablation_main.m` produces a MAT file that stores matched-run IGD,
HV, Spacing, runtime,
archives, seeds, normalization bounds, and reference fronts for the complete
framework, three ablated variants, and standard NSGA-II.
