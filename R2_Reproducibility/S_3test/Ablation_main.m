%% HMD-NSGA-II formal ablation experiment (Scenario S5)
% [Ablation rerun fix-01 | 2026-09-25]
% Reasons for this revision:
%   1) align IGD/HV/Spacing with the normalized S1--S8 evaluation;
%   2) ensure that each ablation disables exactly one proposed mechanism;
%   3) use matched optimization seeds across configurations;
%   4) save raw fronts and per-run checkpoints for reproducibility.

clc;
clear;
close all;

%% 0. Paths and formal experiment settings
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end
project_root = fileparts(script_dir);
comparison_dir = fullfile(project_root, 'Compare_nsga_slot');
addpath(script_dir);
addpath(comparison_dir, '-end');

required_functions = {'EvaluateParticle1', 'OUSNSGA_II', ...
    'FindAllFronts', 'FindNonDominatedSolutions', ...
    'calculateIGD', 'calculateHV'};
for k = 1:numel(required_functions)
    assert(~isempty(which(required_functions{k})), ...
        'Ablation:MissingFunction', ...
        'Required function not found: %s', required_functions{k});
end

nSlots = 10;
num_stat_runs = 30;
shadow_std_dev = struct('LoS', 3.0, 'NLoS', 8.29);

% Scenario S5: I=80, M=4, 100 m x 100 m.
scenario_id = 5;
scenario_config = [80, 4, 100, 100];

alg_names = {'HMD_Full', 'HMD_no_ISMM', 'HMD_no_LAMS', ...
    'HMD_no_Hybrid', 'Standard_NSGAII'};
metric_names = {'IGD', 'HV', 'Spacing', 'Runtime'};

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
output_root = fullfile(script_dir, 'AblationResults_Output');
output_folder = fullfile(output_root, ['S5_', timestamp]);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
checkpoint_file = fullfile(output_folder, 'Ablation_S5_CHECKPOINT.mat');
result_file = fullfile(output_folder, ...
    ['Ablation_S5_Results_', timestamp, '.mat']);

fprintf('Formal S5 ablation experiment\n');
fprintf('Output folder: %s\n', output_folder);
fprintf('Runs=%d, slots=%d, population=100, generations=200\n', ...
    num_stat_runs, nSlots);

%% 1. Problem and algorithm parameters
problem = struct();
problem.objFunc = @EvaluateParticle1;
problem.Tslot = 5;
problem.systemTotalBandwidth = 225e6;
problem.nObj = 2;
problem.nTerminals = scenario_config(1);
problem.nFogNodes = scenario_config(2);
problem.area = [0, scenario_config(3); 0, scenario_config(4)];

params_base = struct('N', 100, 'T_max', 200, ...
    'pc', 0.9, 'pm', 0.05, 'mu', 20, 'mum', 20);
params_base.pm_min = 0.01;
params_base.pm_max = 0.10;
params_base.mum_min = 5;
params_base.mum_max = 20;
params_base.pm_cont_coeff = 0.8;
params_base.pm_disc_coeff = 1.2;
params_base.mum_cont_coeff = 1.1;
params_base.mum_disc_coeff = 0.9;

all_results = struct();
for a_idx = 1:numel(alg_names)
    name = alg_names{a_idx};
    for m_idx = 1:numel(metric_names)
        all_results.(name).(metric_names{m_idx}) = ...
            NaN(num_stat_runs, nSlots);
    end
end

% Raw valid objective fronts: algorithm x slot x run.
current_run_archives = cell(numel(alg_names), nSlots, num_stat_runs);
% Per-run unions used to construct matched empirical reference fronts.
all_objs_for_pf_star = cell(num_stat_runs, nSlots);
environment_records = cell(num_stat_runs, 1);
optimization_seeds = zeros(num_stat_runs, nSlots);
completed_runs = 0;

%% 2. Optimization runs
for run_idx = 1:num_stat_runs
    environment_seed = 5000 + run_idx;
    rng(environment_seed, 'twister');

    % One Monte Carlo environment per run, shared by every configuration
    % and held fixed over the ten slots. This isolates cross-slot reuse.
    problem.terminalProperties.positions = [ ...
        problem.area(1,2) * rand(problem.nTerminals, 1), ...
        problem.area(2,2) * rand(problem.nTerminals, 1)];
    problem.terminalProperties.Pt_dbm = ...
        linspace(10, 15, problem.nTerminals);
    problem.terminalProperties.fc = ...
        linspace(2.4e9, 5.8e9, problem.nTerminals);
    problem.fogNodeProperties.cpu_cycle_rate = ...
        linspace(2e9, 5e9, problem.nFogNodes);

    x_c = problem.area(1,2) * rand(problem.nFogNodes, 1);
    y_c = problem.area(2,2) * rand(problem.nFogNodes, 1);
    problem.initial_fog_positions_matrix = [x_c, y_c];
    problem.initial_fog_deployment_flat = ...
        reshape(problem.initial_fog_positions_matrix', 1, []);
    problem.terminalProperties.task_sizes = ...
        0.1e6 + 0.9e6 * rand(1, problem.nTerminals);
    problem.bounds.bandwidth = [ ...
        0.2e6 * ones(1, problem.nTerminals); ...
        10e6 * ones(1, problem.nTerminals)];
    problem.fixed_shadow_LoS_val = ...
        shadow_std_dev.LoS * randn(1, problem.nTerminals);
    problem.fixed_shadow_NLoS_val = ...
        shadow_std_dev.NLoS * randn(1, problem.nTerminals);

    environment_records{run_idx} = struct( ...
        'EnvironmentSeed', environment_seed, ...
        'TerminalPositions', problem.terminalProperties.positions, ...
        'TransmitPower_dBm', problem.terminalProperties.Pt_dbm, ...
        'CarrierFrequency_Hz', problem.terminalProperties.fc, ...
        'TaskSizes_bits', problem.terminalProperties.task_sizes, ...
        'FogCpuRates', problem.fogNodeProperties.cpu_cycle_rate, ...
        'InitialFogPositions', problem.initial_fog_positions_matrix, ...
        'ShadowLoS_dB', problem.fixed_shadow_LoS_val, ...
        'ShadowNLoS_dB', problem.fixed_shadow_NLoS_val);

    LastSlotArchives = cell(1, numel(alg_names));

    for t = 1:nSlots
        fprintf('\nRun %d/%d, slot %d/%d\n', ...
            run_idx, num_stat_runs, t, nSlots);
        combined_objs_this_slot = zeros(0, problem.nObj);

        % All configurations start the slot from the same random stream.
        % At t=1, HMD-Full and w/o ISMM are therefore exactly matched.
        optimization_seed = 100000 * run_idx + 100 * t;
        optimization_seeds(run_idx, t) = optimization_seed;

        for a_idx = 1:numel(alg_names)
            name = alg_names{a_idx};
            rng(optimization_seed, 'twister');

            p = params_base;
            mem = LastSlotArchives{a_idx};

            switch name
                case 'HMD_Full'
                    p.memory_ratio = 0.1;
                    p.adaptive_enabled = true;
                    p.hybrid_enabled = true;

                case 'HMD_no_ISMM'
                    p.memory_ratio = 0;
                    p.adaptive_enabled = true;
                    p.hybrid_enabled = true;
                    mem = [];

                case 'HMD_no_LAMS'
                    % Disable both adaptive scheduling and layer-specific
                    % scaling, leaving fixed uniform mutation parameters.
                    p.memory_ratio = 0.1;
                    p.adaptive_enabled = false;
                    p.hybrid_enabled = true;
                    p.pm_cont_coeff = 1;
                    p.pm_disc_coeff = 1;
                    p.mum_cont_coeff = 1;
                    p.mum_disc_coeff = 1;

                case 'HMD_no_Hybrid'
                    p.memory_ratio = 0.1;
                    p.adaptive_enabled = true;
                    p.hybrid_enabled = false;

                case 'Standard_NSGAII'
                    p.memory_ratio = 0;
                    p.adaptive_enabled = false;
                    p.hybrid_enabled = false;
                    p.pm_cont_coeff = 1;
                    p.pm_disc_coeff = 1;
                    p.mum_cont_coeff = 1;
                    p.mum_disc_coeff = 1;
                    mem = [];

                otherwise
                    error('Ablation:UnknownConfiguration', ...
                        'Unknown configuration: %s', name);
            end

            tic;
            Pop = OUSNSGA_II(problem, p, mem);
            all_results.(name).Runtime(run_idx, t) = toc;

            Archive = getFirstFront(FindAllFronts(Pop));
            Archive = filterValidArchive(Archive);
            LastSlotArchives{a_idx} = Archive;

            objs = getObjectivesMatrix(Archive, problem.nObj);
            current_run_archives{a_idx, t, run_idx} = objs;
            if ~isempty(objs)
                combined_objs_this_slot = ...
                    [combined_objs_this_slot; objs]; %#ok<AGROW>
            end
        end

        all_objs_for_pf_star{run_idx, t} = combined_objs_this_slot;
    end

    completed_runs = run_idx;
    save(checkpoint_file, 'completed_runs', 'all_results', ...
        'current_run_archives', 'all_objs_for_pf_star', ...
        'environment_records', 'optimization_seeds', 'alg_names', ...
        'metric_names', 'params_base', 'scenario_id', ...
        'scenario_config', 'nSlots', 'num_stat_runs', '-v7.3');
    fprintf('Checkpoint saved after run %d: %s\n', ...
        run_idx, checkpoint_file);
end

%% 3. Metric calculation in the shared normalized objective space
fprintf('\nAll optimization runs completed. Calculating metrics...\n');

normalization_bounds = repmat(struct( ...
    'IdealPoint', [], 'WorstPoint', [], 'ObjectiveRange', [], ...
    'HVReferencePoint', [1.1, 1.1]), nSlots, 1);
reference_fronts = cell(num_stat_runs, nSlots);
normalized_reference_fronts = cell(num_stat_runs, nSlots);

for t = 1:nSlots
    all_objs_this_slot = zeros(0, problem.nObj);
    for r = 1:num_stat_runs
        objs = all_objs_for_pf_star{r, t};
        objs = filterValidObjectives(objs, problem.nObj);
        if ~isempty(objs)
            all_objs_this_slot = [all_objs_this_slot; objs]; %#ok<AGROW>
        end
    end

    if isempty(all_objs_this_slot)
        error('Ablation:NoValidObjectives', ...
            'No valid objective vectors were found at slot %d.', t);
    end

    ideal_point = min(all_objs_this_slot, [], 1);
    worst_point = max(all_objs_this_slot, [], 1);
    objective_range = worst_point - ideal_point;
    objective_range(objective_range < 1e-12) = 1;
    hv_reference_point = [1.1, 1.1];

    normalization_bounds(t).IdealPoint = ideal_point;
    normalization_bounds(t).WorstPoint = worst_point;
    normalization_bounds(t).ObjectiveRange = objective_range;
    normalization_bounds(t).HVReferencePoint = hv_reference_point;

    fprintf('Slot %d: ideal=[%.6g %.6g], worst=[%.6g %.6g]\n', ...
        t, ideal_point(1), ideal_point(2), ...
        worst_point(1), worst_point(2));

    for r = 1:num_stat_runs
        raw_union = filterValidObjectives( ...
            all_objs_for_pf_star{r, t}, problem.nObj);
        if isempty(raw_union)
            raw_reference = zeros(0, problem.nObj);
            normalized_reference = zeros(0, problem.nObj);
        else
            ref_idx = FindNonDominatedSolutions(raw_union);
            raw_reference = raw_union(ref_idx, :);
            [~, order] = sort(raw_reference(:, 1));
            raw_reference = raw_reference(order, :);
            normalized_reference = ...
                (raw_reference - ideal_point) ./ objective_range;
        end

        reference_fronts{r, t} = raw_reference;
        normalized_reference_fronts{r, t} = normalized_reference;

        for a_idx = 1:numel(alg_names)
            name = alg_names{a_idx};
            raw_archive = filterValidObjectives( ...
                current_run_archives{a_idx, t, r}, problem.nObj);

            if isempty(raw_archive) || isempty(normalized_reference)
                all_results.(name).IGD(r, t) = NaN;
                all_results.(name).HV(r, t) = NaN;
                all_results.(name).Spacing(r, t) = NaN;
                continue;
            end

            normalized_archive = ...
                (raw_archive - ideal_point) ./ objective_range;
            all_results.(name).IGD(r, t) = ...
                calculateIGD(normalized_archive, normalized_reference);
            all_results.(name).HV(r, t) = ...
                calculateHV(normalized_archive, hv_reference_point);
            all_results.(name).Spacing(r, t) = ...
                calculateNormalizedSpacing(normalized_archive);
        end
    end
end

%% 4. Terminal-slot summary and final save
fprintf('\n%s\n', repmat('=', 1, 105));
fprintf('%-22s | %-20s | %-20s | %-20s\n', ...
    'Configuration', 'IGD (Mean +/- Std)', ...
    'HV (Mean +/- Std)', 'Spacing (Mean +/- Std)');
fprintf('%s\n', repmat('-', 1, 105));

for a_idx = 1:numel(alg_names)
    name = alg_names{a_idx};
    fIGD = all_results.(name).IGD(:, end);
    fHV = all_results.(name).HV(:, end);
    fSP = all_results.(name).Spacing(:, end);
    fprintf('%-22s | %.6f +/- %.6f | %.6f +/- %.6f | %.6f +/- %.6f\n', ...
        name, mean(fIGD, 'omitnan'), std(fIGD, 0, 'omitnan'), ...
        mean(fHV, 'omitnan'), std(fHV, 0, 'omitnan'), ...
        mean(fSP, 'omitnan'), std(fSP, 0, 'omitnan'));
end
fprintf('%s\n', repmat('=', 1, 105));

save(result_file, 'all_results', 'current_run_archives', ...
    'all_objs_for_pf_star', 'reference_fronts', ...
    'normalized_reference_fronts', 'normalization_bounds', ...
    'environment_records', 'optimization_seeds', 'alg_names', ...
    'metric_names', 'params_base', 'scenario_id', 'scenario_config', ...
    'nSlots', 'num_stat_runs', 'completed_runs', 'output_folder', ...
    'result_file', '-v7.3');

fprintf('\nFinal result saved: %s\n', result_file);
fprintf('\nFormal S5 ablation experiment completed successfully.\n');

%% Local helper functions
function first_front_archive = getFirstFront(fronts_cell_array)
    if ~isempty(fronts_cell_array) && ~isempty(fronts_cell_array{1})
        first_front_archive = fronts_cell_array{1};
    else
        first_front_archive = [];
    end
end

function archive = filterValidArchive(archive)
    if isempty(archive)
        return;
    end
    objectives = vertcat(archive.Objectives);
    valid_mask = ~any(objectives >= 1e9 | isnan(objectives) | ...
        isinf(objectives), 2);
    archive = archive(valid_mask);
end

function obj_matrix = getObjectivesMatrix(archive_struct, n_obj)
    if isempty(archive_struct)
        obj_matrix = zeros(0, n_obj);
    else
        obj_matrix = vertcat(archive_struct.Objectives);
        obj_matrix = filterValidObjectives(obj_matrix, n_obj);
    end
end

function obj_matrix = filterValidObjectives(obj_matrix, n_obj)
    if isempty(obj_matrix)
        obj_matrix = zeros(0, n_obj);
        return;
    end
    if size(obj_matrix, 2) ~= n_obj
        error('Ablation:ObjectiveDimensionMismatch', ...
            'Expected %d objectives, but received %d.', ...
            n_obj, size(obj_matrix, 2));
    end
    valid_mask = ~any(obj_matrix >= 1e9 | isnan(obj_matrix) | ...
        isinf(obj_matrix), 2);
    obj_matrix = obj_matrix(valid_mask, :);
end

function spacing = calculateNormalizedSpacing(points)
    n_points = size(points, 1);
    if n_points < 2
        spacing = NaN;
        return;
    end
    distance_matrix = pdist2(points, points);
    distance_matrix(1:n_points+1:end) = inf;
    nearest_distances = min(distance_matrix, [], 2);
    spacing = sqrt(sum((nearest_distances - ...
        mean(nearest_distances)).^2) / (n_points - 1));
end
