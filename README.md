## 📂 Folder Structure: `Compare_nsga_slot`

This directory serves as the **core experimental framework** for the research. It contains the primary simulation engine, benchmark algorithms, and evaluation metrics used to generate the comparative data (such as the Radar Charts and Statistical Distributions) presented in the paper.

### 1. Main Execution & Scripting
* **`M3.m`**: The primary entry point for the system experiments. This script initializes the IIoT environment, configures scenario parameters (S1–S6), and manages the multi-slot optimization process across the entire simulation.
* **`mean_slot.m` / `mean_slot2.m`**: Scripts responsible for managing temporal dynamics across consecutive time slots, specifically handling average performance data and inter-slot transitions.
* **`radar.m`**: The visualization script used to generate the **Radar Charts**, illustrating the multi-dimensional trade-off balance between IGD, HV, and Spacing.
* **`plot5.m` / `plotS.m`**: Specialized plotting utilities for generating Pareto front comparisons and statistical distribution graphs.

### 2. Optimization Algorithm Library
This folder implements the proposed framework alongside several baseline multi-objective evolutionary algorithms (MOEAs) used for performance benchmarking:
* **`MyNSGA_II.m`**: The implementation of the proposed **HMD-NSGA-II** algorithm, featuring the Inter-Slot Memory Mechanism (ISMM) and Adaptive Mutation Strategy (AMS).
* **`NSGA_III.m`**: A reference point-based algorithm used to evaluate performance in high-dimensional objective spaces.
* **`MOEAD.m`**: A decomposition-based MOEA that serves as a baseline for convergence and diversity testing.
* **`MOPSO.m` / `MOPO.m`**: Implementations of Multi-Objective Particle Swarm Optimization variants used to assess improvements over traditional heuristics.

### 3. Metric Evaluation & Helper Functions
These functions quantify the quality of the identified Pareto-optimal sets according to the metrics defined in the study:
* **`calculateHV.m`**: Computes the **Hypervolume (HV)**, measuring the volume of the objective space covered by the solution set.
* **`calculateIGD.m` / `IGD.m`**: Calculates the **Inverted Generational Distance (IGD)** to assess both convergence and diversity.
* **`calculateSpread.m`**: Measures the distribution uniformity of solutions, also referred to as the **Spacing value**.
* **`CalculateTmax.m`**: Evaluates the system's real-time processing capability by calculating the maximum completion time within a single time slot.
* **`FindNonDominated.m` / `FindAllFronts.m`**: Fundamental sorting utilities that identify non-dominated solutions and organize the population into Pareto ranks.
* **`EvaluateParticle.m`**: The fitness evaluation engine that interfaces with the **MO-MINLP model** to calculate system latency ($G_1$) and energy consumption ($G_2$).
