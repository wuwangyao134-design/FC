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

## 📂 Folder Structure: `S_3test`

This directory is primarily used for **Ablation Studies** and **Robustness Testing** to evaluate the individual contributions of each module within our proposed framework. It features our final optimized method: **OUSNSGA-II**.

### 1. Proposed Method & Evolution
* **`OUSNSGA_II.m`**: The main implementation of our final proposed method (**OUSNSGA-II**). This script incorporates the complete suite of enhancements, including the Inter-Slot Memory Mechanism (ISMM) and Adaptive Mutation Strategy (AMS), for optimal resource allocation in IIoT fog networks.
* **`IMyNSGA_II.m`**: An **early-stage preliminary version** of our method. This file represents the baseline development phase of the OUSNSGA-II framework before the final optimizations and refinements were integrated.

### 2. Ablation Study Framework
* **`Ablation_main.m`**: The master script for conducting ablation experiments. It systematically enables or disables specific modules (ISMM, AMS, HO) to quantitatively verify the performance gains elicited by each component.
* **`Plot_Ablation_Trend.m`**: A visualization utility used to generate the **Ablation Performance Curves** (such as IGD and HV trajectories) as presented in the experimental analysis.
* **`EvaluateParticle.m`**: A specialized fitness evaluation function tailored for ablation scenarios, ensuring accurate mapping of system latency and energy consumption objectives.

### 3. Baseline & Comparative Algorithms
* **`AMOPSO.m`**: Implementation of the Adaptive Multi-Objective Particle Swarm Optimization algorithm, used as a non-genetic baseline for comparison.
* **`DNSGA_II.m`**: A discrete version of the standard NSGA-II used to test the effectiveness and necessity of our hybrid encoding approach.
* **`MOEAD.m`**: A decomposition-based MOEA baseline specifically configured for testing against our framework under high-dimensional decision spaces.

### 4. Utility Scripts
* **`FindAllFronts.m`**: Performs non-dominated sorting to categorize the population into various Pareto ranks.
* **`renewable.m`**: A utility script that manages environment refreshing and task load states during multi-slot simulation runs to maintain system stability.
* 
## 📂 Experimental datasets S1-S6


🚀 How to Run the Experiments
Follow these steps to reproduce the results and performance metrics presented in the paper:

Step 1: Download and Preparation

Clone the Repository: Download the complete project folder from GitHub to your local machine.


Environment Setup: Ensure you have MATLAB R2024b or a later version installed.

Add to Path: Open MATLAB and add the project folders (Compare_nsga_slot and S_3test) to your MATLAB path to ensure all helper functions are accessible.

Step 2: Running the System Experiments

Execute Main Script: Run the M3.m script located in the Compare_nsga_slot folder.


Optimization Process: This script will initialize the IIoT environment and execute our proposed OUSNSGA-II algorithm alongside other baseline algorithms (NSGA-III, MOEA/D, etc.).




Results Generation: After the simulation completes, the script will automatically calculate the IGD, HV, and Spacing metrics and generate the performance comparison plots.





Step 3: Running Ablation Studies (Optional)

Verify Components: To verify the effectiveness of the ISMM and AMS modules, navigate to the S_3test folder and run Ablation_main.m.


Trend Analysis: Use Plot_Ablation_Trend.m to visualize how the OUSNSGA-II method outperforms its preliminary versions and other variants
