function validation = Run_Micro_Exact_Validation_3x2(run_mode)
%RUN_MICRO_EXACT_VALIDATION_3X2 Validate a contended 3-MTD/2-FN instance.
%
% First smoke test:
%   Run_Micro_Exact_Validation_3x2('quick')
%
% Formal reviewer experiment:
%   Run_Micro_Exact_Validation_3x2('formal')
%
% The exact model uses the same EvaluateParticle physical parameters,
% including the +5 dB receiver noise figure.  A single fixed stochastic
% channel realization is shared by BARON and every algorithmic run.
%
% [Exact-validation addition 03 | 2026-09-24]
% Reason: add an independent, more demanding micro instance in which three
% MTDs share two heterogeneous FNs.  Every association, FCFS arrival order,
% and queue branch is enumerated; no heuristic ranking or Top-K pruning is
% used in the BARON reference calculation.

    if nargin < 1
        run_mode = 'quick';
    end
    run_mode = lower(string(run_mode));
    if ~ismember(run_mode, ["quick", "formal"])
        error('MicroValidation:InvalidMode', ...
              'run_mode must be ''quick'' or ''formal''.');
    end

    clc;
    close all;
    start_time = datetime('now');
    this_folder = fileparts(mfilename('fullpath'));
    algorithm_folder = fileparts(this_folder);
    addpath(this_folder);
    addpath(algorithm_folder);

    assert_dependency('baron');
    assert_dependency('baron_unpruned_micro_solver_3x2');
    assert_dependency('EvaluateParticle');
    assert_dependency('HMD_NSGA_II');
    assert_dependency('REMODDPG_Baseline');

    cfg = experiment_configuration(run_mode);
    fprintf('\n3x2 contended micro exact-validation mode: %s\n', upper(run_mode));
    fprintf('Independent heuristic runs: %d\n', cfg.NumRuns);
    fprintf('Exact delay weights: %s\n', mat2str(cfg.DelayWeights', 3));

    %% 1. Create one fixed micro instance using the original structure.
    problem_base = struct();
    problem_base.objFunc = @EvaluateParticle;
    problem_base.Tslot = 5;
    problem_base.systemTotalBandwidth = 225e6;
    problem_base.nObj = 2;

    rng(cfg.InstanceSeed, 'twister');
    problem = create_fog_problem( ...
        {cfg.NumTerminals, cfg.NumFogNodes, cfg.RegionWidth, cfg.RegionHeight}, ...
        problem_base);

    % Retain heterogeneous FN rates. Since I=3>M=2, at least two MTDs must
    % share an FN, so the exact experiment preserves association, FCFS
    % ordering, and queue-contention effects.

    rng(cfg.ShadowSeed, 'twister');
    problem.fixed_shadow_LoS_val = ...
        cfg.ShadowStdLoS * randn(1, problem.nTerminals);
    problem.fixed_shadow_NLoS_val = ...
        cfg.ShadowStdNLoS * randn(1, problem.nTerminals);

    %% 2. Certified, fully enumerated BARON reference solutions.
    exact = baron_unpruned_micro_solver_3x2( ...
        problem, cfg.DelayWeights, cfg.BARONOptions);
    if ~exact.AllScalarizationsCertified
        error('MicroValidation:UncertifiedExactResult', ...
              'At least one scalarized reference problem lacks a certificate.');
    end

    %% 3. Original HMD-NSGA-II and adapted MODDPG implementations.
    hmd_params = hmd_parameter_settings(cfg);
    drl_params = remoddpg_parameter_settings(cfg);

    hmd_fronts = cell(cfg.NumRuns, 1);
    drl_fronts = cell(cfg.NumRuns, 1);
    hmd_runtime = nan(cfg.NumRuns, 1);
    drl_training_time = nan(cfg.NumRuns, 1);
    drl_inference_time = nan(cfg.NumRuns, 1);

    for run_idx = 1:cfg.NumRuns
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Algorithmic run %d/%d\n', run_idx, cfg.NumRuns);
        fprintf('------------------------------------------------------------\n');

        rng(cfg.HMDSeedBase + run_idx, 'twister');
        timer_hmd = tic;
        pop_hmd = HMD_NSGA_II(problem, hmd_params, []);
        hmd_runtime(run_idx) = toc(timer_hmd);
        hmd_fronts{run_idx} = extract_feasible_front(pop_hmd, problem);

        rng(cfg.DRLSeedBase + run_idx, 'twister');
        [pop_drl, ~, drl_info] = REMODDPG_Baseline( ...
            problem, drl_params, []);
        drl_training_time(run_idx) = drl_info.TrainingTime;
        drl_inference_time(run_idx) = drl_info.InferenceTime;
        drl_fronts{run_idx} = extract_feasible_front(pop_drl, problem);
    end

    %% 4. Preference-wise normalized optimality gaps.
    hmd_gaps = front_optimality_gaps(hmd_fronts, exact);
    drl_gaps = front_optimality_gaps(drl_fronts, exact);

    hmd_run_mean = mean(hmd_gaps, 2, 'omitnan');
    drl_run_mean = mean(drl_gaps, 2, 'omitnan');

    algorithm = ["HMD-NSGA-II"; "MODDPG (adapted)"];
    mean_gap_percent = [mean(hmd_run_mean, 'omitnan'); ...
                        mean(drl_run_mean, 'omitnan')];
    std_gap_percent = [std(hmd_run_mean, 'omitnan'); ...
                       std(drl_run_mean, 'omitnan')];
    max_gap_percent = [max(hmd_gaps, [], 'all', 'omitnan'); ...
                       max(drl_gaps, [], 'all', 'omitnan')];
    valid_run_rate_percent = [100*mean(~cellfun(@isempty, hmd_fronts)); ...
                              100*mean(~cellfun(@isempty, drl_fronts))];
    mean_online_runtime_s = [mean(hmd_runtime, 'omitnan'); ...
                             mean(drl_inference_time, 'omitnan')];

    summary = table(algorithm, mean_gap_percent, std_gap_percent, ...
        max_gap_percent, valid_run_rate_percent, mean_online_runtime_s, ...
        'VariableNames', {'Algorithm','MeanGapPercent','StdGapPercent', ...
        'MaxGapPercent','ValidRunRatePercent','MeanOnlineRuntimeSeconds'});

    fprintf('\n==================== GAP SUMMARY ====================\n');
    disp(summary);
    fprintf(['Gap definition: 100 times the additive difference between ', ...
             'the best heuristic and certified optimum in the common ', ...
             'normalized weighted objective space.\n']);

    %% 5. Save auditable data, a CSV summary, and a representative plot.
    output_folder = fullfile(algorithm_folder, ...
        'MicroExactValidation_3x2_Output', ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    validation = struct();
    validation.RunMode = char(run_mode);
    validation.Configuration = cfg;
    validation.Problem = problem;
    validation.Exact = exact;
    validation.HMDFronts = hmd_fronts;
    validation.MODDPGFronts = drl_fronts;
    validation.HMDGapPercent = hmd_gaps;
    validation.MODDPGGapPercent = drl_gaps;
    validation.HMDRuntime = hmd_runtime;
    validation.MODDPGTrainingTime = drl_training_time;
    validation.MODDPGInferenceTime = drl_inference_time;
    validation.Summary = summary;
    validation.StartTime = start_time;
    validation.EndTime = datetime('now');

    mat_path = fullfile(output_folder, ...
        sprintf('Micro_Exact_Validation_3x2_%s.mat', upper(char(run_mode))));
    csv_path = fullfile(output_folder, 'Micro_Exact_Gap_Summary_3x2.csv');
    save(mat_path, 'validation', '-v7.3');
    writetable(summary, csv_path);
    make_validation_plot(exact, hmd_fronts, drl_fronts, ...
        hmd_run_mean, drl_run_mean, output_folder);

    fprintf('\nSaved validation data: %s\n', mat_path);
    fprintf('Saved gap summary:     %s\n', csv_path);
end

function cfg = experiment_configuration(run_mode)
    cfg = struct();
    cfg.NumTerminals = 3;
    cfg.NumFogNodes = 2;
    % The compact work cell has diagonal < 1 m. Consequently, the model's
    % d=max(raw distance,1 m) rule makes every deployment point physically
    % equivalent, while association, heterogeneous service rates, queueing,
    % and continuous bandwidth allocation remain active. This exact bound
    % propagation is not heuristic pruning.
    cfg.RegionWidth = 0.5;
    cfg.RegionHeight = 0.5;
    cfg.InstanceSeed = 9401;
    cfg.ShadowSeed = 9402;
    cfg.ShadowStdLoS = 3.0;
    cfg.ShadowStdNLoS = 8.29;
    cfg.HMDSeedBase = 9500;
    cfg.DRLSeedBase = 9600;

    cfg.BARONOptions = struct( ...
        'MaxTimePerCase', 120, ...
        'AbsoluteTolerance', 1e-7, ...
        'RelativeTolerance', 1e-7, ...
        'PrintLevel', 0, ...
        'ProgressEvery', 15, ...
        'StopOnIncomplete', false, ...
        'ScratchDirectory', fullfile(tempdir, 'baron_micro_3x2_scratch'));

    if run_mode == "quick"
        % Smoke test only; these reduced algorithm budgets must not be used
        % for the numerical values reported in the response letter.
        cfg.NumRuns = 1;
        % The contended 3x2 instance certifies the two objective endpoints.
        % The balanced certified preference is already supplied by the
        % complementary 2x2 experiment.
        cfg.DelayWeights = [0; 1];
        cfg.PopulationSize = 30;
        cfg.Generations = 30;
        cfg.DRLEpisodes = 1000;
        cfg.BARONOptions.MaxTimePerCase = 20;
        cfg.BARONOptions.ProgressEvery = 12;
    else
        cfg.BARONOptions.MaxTimePerCase = 300;
        cfg.NumRuns = 30;
        % [Exact-validation modification 02 | 2026-09-24]
        % [Exact-validation modification 05 | 2026-09-24]
        % On this contended instance, report only the two globally certified
        % endpoint preferences. The complementary 2x2 validation provides
        % the certified balanced preference w=0.5.
        cfg.DelayWeights = [0; 1];
        cfg.PopulationSize = 100;
        cfg.Generations = 200;
        cfg.DRLEpisodes = 20000;
    end
end

function params = hmd_parameter_settings(cfg)
    params = struct('N', cfg.PopulationSize, ...
                    'T_max', cfg.Generations, ...
                    'pc', 0.9, 'pm', 0.05, 'mu', 20, 'mum', 20, ...
                    'pm_cont_coeff', 0.8, 'pm_disc_coeff', 1.2, ...
                    'mum_cont_coeff', 1.1, 'mum_disc_coeff', 0.9);
    params.adaptive_enabled = true;
    params.hybrid_enabled = true;
    params.memory_ratio = 0.1;
    params.pm_min = 0.01;
    params.pm_max = 0.10;
    params.mum_min = 5;
    params.mum_max = 20;
end

function params = remoddpg_parameter_settings(cfg)
    params = struct();
    params.N = cfg.PopulationSize;
    params.drl_train_episodes = cfg.DRLEpisodes;
    params.drl_batch_size = 64;
    params.drl_hidden = [128, 64];
    params.drl_actor_lr = 1e-4;
    params.drl_critic_lr = 3e-4;
    params.drl_preference_bins = 21;
    params.drl_elites_per_bin = 4;
    params.drl_inference_candidates_per_preference = 8;
    params.drl_inference_noise = [0.02, 0.05, 0.10];
end

function problem = create_fog_problem(scenario, base)
    problem = base;
    problem.nTerminals = scenario{1};
    problem.nFogNodes = scenario{2};
    problem.area = [0 scenario{3}; 0 scenario{4}];

    term_x = problem.area(1,1) + diff(problem.area(1,:)) .* ...
        rand(problem.nTerminals, 1);
    term_y = problem.area(2,1) + diff(problem.area(2,:)) .* ...
        rand(problem.nTerminals, 1);
    problem.terminalProperties.positions = [term_x, term_y];
    problem.terminalProperties.Pt_dbm = ...
        linspace(10, 15, problem.nTerminals);
    problem.terminalProperties.fc = ...
        linspace(2.4e9, 5.8e9, problem.nTerminals);

    problem.fogNodeProperties.cpu_cycle_rate = ...
        linspace(2e9, 5e9, problem.nFogNodes);
    fog_x = problem.area(1,1) + diff(problem.area(1,:)) .* ...
        rand(problem.nFogNodes, 1);
    fog_y = problem.area(2,1) + diff(problem.area(2,:)) .* ...
        rand(problem.nFogNodes, 1);
    problem.initial_fog_positions_matrix = [fog_x, fog_y];
    problem.initial_fog_deployment_flat = ...
        reshape(problem.initial_fog_positions_matrix', 1, []);

    problem.terminalProperties.task_sizes = ...
        0.1e6 + (1e6-0.1e6) .* rand(1, problem.nTerminals);
    problem.bounds.bandwidth = [0.2e6*ones(1, problem.nTerminals); ...
                                10e6*ones(1, problem.nTerminals)];
end

function front = extract_feasible_front(population, problem)
    if isempty(population)
        front = zeros(0, 2);
        return;
    end
    evaluated = EvaluateParticle(population, problem);
    valid = evaluated.IsFeasible(:) & ...
            all(isfinite(evaluated.RawObjectives), 2);
    objectives = evaluated.RawObjectives(valid, :);
    if isempty(objectives)
        front = zeros(0, 2);
        return;
    end
    objectives = unique(objectives, 'rows', 'stable');
    front = nondominated_rows(objectives);
    front = sortrows(front, 1);
end

function gaps = front_optimality_gaps(fronts, exact)
    n_runs = numel(fronts);
    n_weights = numel(exact.Weights);
    gaps = nan(n_runs, n_weights);

    for r = 1:n_runs
        if isempty(fronts{r})
            continue;
        end
        normalized = (fronts{r} - exact.IdealPoint) ./ ...
                     exact.ObjectiveRange;
        for k = 1:n_weights
            weight = [exact.Weights(k), 1-exact.Weights(k)];
            heuristic_best = min(normalized * weight(:));
            additive_gap = heuristic_best - exact.ExactScalarScores(k);
            % Negative values within solver tolerance are numerical noise.
            gaps(r, k) = 100 * max(0, additive_gap);
        end
    end
end

function make_validation_plot(exact, hmd_fronts, drl_fronts, ...
        hmd_run_mean, drl_run_mean, output_folder)
    hmd_index = representative_index(hmd_run_mean, hmd_fronts);
    drl_index = representative_index(drl_run_mean, drl_fronts);

    figure_handle = figure('Color', 'w', 'Position', [100 100 760 540]);
    hold on;
    plot(exact.ExactFront(:,1), exact.ExactFront(:,2), ...
    'ko', ...
    'LineStyle', 'none', ...
    'LineWidth', 1.6, ...
    'MarkerSize', 7, ...
    'MarkerFaceColor', 'w', ...
    'DisplayName', 'Certified exact endpoints');
    if ~isnan(hmd_index)
        plot(hmd_fronts{hmd_index}(:,1), hmd_fronts{hmd_index}(:,2), ...
            'r-s', 'LineWidth', 1.2, 'MarkerSize', 5, ...
            'DisplayName', 'HMD-NSGA-II');
    end
    if ~isnan(drl_index)
        plot(drl_fronts{drl_index}(:,1), drl_fronts{drl_index}(:,2), ...
            'b-^', 'LineWidth', 1.2, 'MarkerSize', 5, ...
            'DisplayName', 'MODDPG');
    end
    grid on;
    box on;
    xlabel('Average completion time (s)');
    ylabel('Energy consumption (J)');
    title('Unpruned 3-MTD/2-FN Global Validation');
    legend('Location', 'best');
    hold off;

    exportgraphics(figure_handle, ...
        fullfile(output_folder, 'Micro_Exact_Front_Comparison_3x2.pdf'), ...
        'ContentType', 'vector');
    savefig(figure_handle, ...
        fullfile(output_folder, 'Micro_Exact_Front_Comparison_3x2.fig'));
end

function index = representative_index(run_values, fronts)
    valid = find(isfinite(run_values) & ~cellfun(@isempty, fronts));
    if isempty(valid)
        index = NaN;
        return;
    end
    target = median(run_values(valid));
    [~, local] = min(abs(run_values(valid)-target));
    index = valid(local);
end

function nd = nondominated_rows(values)
    n = size(values, 1);
    keep = true(n, 1);
    for i = 1:n
        for j = 1:n
            if i ~= j && all(values(j,:) <= values(i,:)) && ...
                    any(values(j,:) < values(i,:))
                keep(i) = false;
                break;
            end
        end
    end
    nd = values(keep, :);
end

function assert_dependency(function_name)
    if isempty(which(function_name))
        error('MicroValidation:MissingDependency', ...
            'Required function %s is not on the MATLAB path.', function_name);
    end
end
