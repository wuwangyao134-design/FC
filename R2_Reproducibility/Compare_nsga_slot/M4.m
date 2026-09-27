% =========================================================================
% NSGA-II、MyNSGA-II、MOPO、INSGA-II、NSGA-III 主执行脚本 (main_compare.m)
% 描述: 对比原始 NSGA-II、带有记忆机制及自适应参数的 MyNSGA-II、
%       MOPO (2026)、INSGA-II (2022) 和 NSGA-III。
% =========================================================================
clc;
clear;
close all;

%% 0. 全局参数设定
nSlots = 10; 
num_stat_runs = 30; % 建议在正式跑数据时改为 30
shadow_std_dev.LoS = 3.0;  
shadow_std_dev.NLoS = 8.29; 

% --- 场景目录：[ID, I, M, width, height] ---
scenario_catalog = {
    1,  20,  2, 100, 100;
    2,  20,  2, 200, 200;
    3,  45,  3, 100, 100;
    4,  45,  3, 200, 200;
    5,  80,  4, 100, 100;
    6,  80,  4, 200, 200;
    7, 150,  8, 300, 300;
    8, 300, 15, 500, 500
};

% 正式实验时只需修改这一行，取值范围为 1--8。
scenario_id = 1;

catalog_ids = cell2mat(scenario_catalog(:, 1));
catalog_row = find(catalog_ids == scenario_id, 1);
if isempty(catalog_row)
    error('M4:InvalidScenarioID', 'scenario_id 必须是 1--8 之间的整数。');
end
experimental_scenarios = cell(1, 1);
experimental_scenarios{1} = scenario_catalog(catalog_row, 2:5);
num_scenarios = 1;

% --- 更新算法名称列表 ---
alg_names_for_results = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III','DRL_Baseline'};
metric_names = {'IGD', 'HV', 'Spacing', 'Spread', 'Runtime', 'NumFeasibleSolutions','Tmax'}; 

all_scenario_results = struct();
for a_name = alg_names_for_results
    for m_name = metric_names
        all_scenario_results.(a_name{1}).(m_name{1}) = cell(num_scenarios, 1);
        for s = 1:num_scenarios
            all_scenario_results.(a_name{1}).(m_name{1}){s} = zeros(num_stat_runs, nSlots);
        end
    end
end

% 输出文件夹配置
output_base_folder = pwd; 
output_main_folder = fullfile(output_base_folder, 'ExperimentResults_Output');
output_timestamp_folder = fullfile(output_main_folder, string(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
if ~exist(output_timestamp_folder, 'dir'), mkdir(output_timestamp_folder); end

% [新增] 创建一个全局变量，用于存储所有场景的专属参考前沿，供最后画图使用
all_reference_fronts = cell(num_scenarios, 1);

% 指标统一在归一化目标空间中计算。每个场景、每个时隙使用
% 全部算法和全部独立运行共享的 ideal/worst 边界。
all_normalized_reference_fronts = cell(num_scenarios, 1);
all_metric_normalization_bounds = cell(num_scenarios, 1);

% [新增] 专门用来存放每次运行的“被剪枝的精确解”，留着给审稿人看画图用
all_pruned_exact_pfs = cell(num_stat_runs, 1);

% 保存每个场景、每次独立运行的 REMODDPG 训练诊断信息。
all_REMODDP_training_info = cell(num_scenarios, num_stat_runs);

%% 5. 场景循环
problem_base.objFunc = @EvaluateParticle;
problem_base.Tslot = 5;
problem_base.systemTotalBandwidth = 225e6;
problem_base.nObj = 2; 

for s_idx = 1:num_scenarios
    % 1. 从列表中提取当前场景的配置参数
    current_scenario_config = experimental_scenarios{s_idx};
    % 2. 用场景编号设定种子，确保物理终端位置固定
    rng(scenario_id * 100); 
    problem = CreateFogProblem(current_scenario_config, problem_base); 
    
    % --- [核心：基础参数] ---
    base_params = struct('N', 100, 'T_max', 200, 'pc', 0.9, 'pm', 0.05, 'mu', 20, 'mum', 20, ...
                         'pm_cont_coeff', 0.8, 'pm_disc_coeff', 1.2, ...
                         'mum_cont_coeff', 1.1, 'mum_disc_coeff', 0.9);
        
    params_my_nsga = base_params; 
    params_my_nsga.adaptive_enabled = true; 
    params_my_nsga.hybrid_enabled   = true; 
    params_my_nsga.memory_ratio     = 0.1;  
    params_my_nsga.pm_min = 0.01; 
    params_my_nsga.pm_max = 0.1;
    params_my_nsga.mum_min = 5; 
    params_my_nsga.mum_max = 20;
        
    params_nsga = base_params; 
    
    params_moead = base_params; 
    params_moead.T = 20;             
    params_moead.delta = 0.9;        
    
    params_insga = base_params; 
    params_nsga3 = base_params; 
    params_nsga3.p = 99;
    
    params_drl = base_params;
    params_drl.drl_train_episodes = 20000;
    params_drl.drl_batch_size = 64;
    params_drl.drl_hidden = [128, 64];
    % [REMODDPG修改-01 | 2026-09-22]
    % 原因：修复S1/S2中推理前沿多样性不足；所有场景固定使用同一组参数，
    % 避免为了改善单个场景结果而进行逐场景调参。
    params_drl.drl_preference_bins = 21;
    params_drl.drl_elites_per_bin = 4;
    params_drl.drl_inference_candidates_per_preference = 8;
    params_drl.drl_inference_noise = [0.02, 0.05, 0.10];
    
    current_scenario_all_run_final_fronts = cell(length(alg_names_for_results), nSlots, num_stat_runs);
    all_combined_solutions_for_global_pf = cell(num_stat_runs, nSlots);
   
    %% 6. 统计运行循环
    for run_idx = 1:num_stat_runs
        seed_value = (scenario_id - 1) * 1000 + run_idx;
        rng(seed_value);
        
        % 1.阴影衰落刷新 (蒙特卡洛环境)
        problem.fixed_shadow_LoS_val = shadow_std_dev.LoS * randn(1, problem.nTerminals);
        problem.fixed_shadow_NLoS_val = shadow_std_dev.NLoS * randn(1, problem.nTerminals);
        
        REMODDP_agent = [];
        REMODDP_training_info = [];
        LastSlotArchive_MyNSGA = []; 
        
        % % --- 🌟 修复并恢复精确求解器 ---
        % fprintf('\n【Run %d】动态阴影已刷新 → 正在跑传统的剪枝精确求解器 (耗时较长)...\n', run_idx);
        % [Exact_PF_this_run, time_exact] = exact_enumeration_solver(problem);
        % 
        % % 将当次环境跑出来的精确前沿单独存起来！(不参与IGD计算，纯留档)
        % all_pruned_exact_pfs{run_idx} = Exact_PF_this_run;

        %% 7. 时隙循环
        for t_slot = 1:nSlots
            fprintf('场景%d-运行%d: 执行时隙 %d/%d\n', scenario_id, run_idx, t_slot, nSlots);
            
            % 1. OUSNSGA_II
            tic;
            Pop_MyNSGA = HMD_NSGA_II(problem, params_my_nsga, LastSlotArchive_MyNSGA);
            all_scenario_results.MyNSGA_II.Runtime{s_idx}(run_idx, t_slot) = toc;
            Archive_MyNSGA = getFirstFront(FindAllFronts(Pop_MyNSGA));
            LastSlotArchive_MyNSGA = Archive_MyNSGA;
            t_MyNSGA = [Pop_MyNSGA.Tmax]; t_MyNSGA = t_MyNSGA(t_MyNSGA < 1e8 & ~isnan(t_MyNSGA));
            all_scenario_results.MyNSGA_II.Tmax{s_idx}(run_idx, t_slot) = mean(t_MyNSGA, 'omitnan');
            
            % 2. NSGA-II
            tic;
            Pop_NSGA2 = DNSGA_II(problem, params_nsga, []);
            all_scenario_results.NSGA_II.Runtime{s_idx}(run_idx, t_slot) = toc;
            Archive_NSGA2 = getFirstFront(FindAllFronts(Pop_NSGA2));
            t_NSGA2 = [Pop_NSGA2.Tmax]; t_NSGA2 = t_NSGA2(t_NSGA2 < 1e8 & ~isnan(t_NSGA2));
            all_scenario_results.NSGA_II.Tmax{s_idx}(run_idx, t_slot) = mean(t_NSGA2, 'omitnan');
            
            % 3. MOEA/D
            tic;
            Pop_MOEAD = MOEAD(problem, params_moead, []);
            all_scenario_results.MOEA_D.Runtime{s_idx}(run_idx, t_slot) = toc;
            Archive_MOEAD = Pop_MOEAD;
            t_MOEAD = [Pop_MOEAD.Tmax]; t_MOEAD = t_MOEAD(t_MOEAD < 1e8 & ~isnan(t_MOEAD));
            all_scenario_results.MOEA_D.Tmax{s_idx}(run_idx, t_slot) = mean(t_MOEAD, 'omitnan');
            
            % 4. INSGA-II
            tic;
            Pop_INSGA = INSGA_II(problem, params_insga, []);
            all_scenario_results.INSGA_II.Runtime{s_idx}(run_idx, t_slot) = toc;
            Archive_INSGA = getFirstFront(FindAllFronts(Pop_INSGA));
            t_INSGA = [Pop_INSGA.Tmax]; t_INSGA = t_INSGA(t_INSGA < 1e8 & ~isnan(t_INSGA));
            all_scenario_results.INSGA_II.Tmax{s_idx}(run_idx, t_slot) = mean(t_INSGA, 'omitnan');
            
            % 5. NSGA-III
            tic;
            Pop_NSGA3 = NSGA_III(problem, params_nsga3, []);
            all_scenario_results.NSGA_III.Runtime{s_idx}(run_idx, t_slot) = toc;
            Archive_NSGA3 = getFirstFront(FindAllFronts(Pop_NSGA3));
            t_NSGA3 = [Pop_NSGA3.Tmax]; t_NSGA3 = t_NSGA3(t_NSGA3 < 1e8 & ~isnan(t_NSGA3));
            all_scenario_results.NSGA_III.Tmax{s_idx}(run_idx, t_slot) = mean(t_NSGA3, 'omitnan');
            
            % 6. REMODDPG
            [Pop_DRL, REMODDP_agent, remoddp_info] = REMODDPG_Baseline(problem, params_drl, REMODDP_agent);
            all_scenario_results.DRL_Baseline.Runtime{s_idx}(run_idx, t_slot) = remoddp_info.InferenceTime;
            if remoddp_info.WasTrained
                REMODDP_training_info = remoddp_info;
                all_REMODDP_training_info{s_idx, run_idx} = remoddp_info;
            end
            Archive_DRL = getFirstFront(FindAllFronts(Pop_DRL));
            t_DRL = [Pop_DRL.Tmax]; t_DRL = t_DRL(t_DRL < 1e8 & ~isnan(t_DRL));
            all_scenario_results.DRL_Baseline.Tmax{s_idx}(run_idx, t_slot) = mean(t_DRL, 'omitnan');
            
            % --- 收集前沿数据用于指标计算 ---
            current_scenario_all_run_final_fronts{1, t_slot, run_idx} = getObjectivesMatrix(Archive_MyNSGA);
            current_scenario_all_run_final_fronts{2, t_slot, run_idx} = getObjectivesMatrix(Archive_NSGA2);
            current_scenario_all_run_final_fronts{3, t_slot, run_idx} = getObjectivesMatrix(Archive_MOEAD);
            current_scenario_all_run_final_fronts{4, t_slot, run_idx} = getObjectivesMatrix(Archive_INSGA);
            current_scenario_all_run_final_fronts{5, t_slot, run_idx} = getObjectivesMatrix(Archive_NSGA3); 
            current_scenario_all_run_final_fronts{6, t_slot, run_idx} = getObjectivesMatrix(Archive_DRL); 
            
            % --- 收集【本次运行、本时隙】所有算法产生的解，用于构建专属PF* ---
            temp_objs = {getObjectivesMatrix(Archive_MyNSGA), getObjectivesMatrix(Archive_NSGA2), ...
                         getObjectivesMatrix(Archive_MOEAD), getObjectivesMatrix(Archive_INSGA), ...
                         getObjectivesMatrix(Archive_NSGA3), getObjectivesMatrix(Archive_DRL)};
            % 注意：这里不再是全局汇总，而是严格对应 run_idx 存储
            all_combined_solutions_for_global_pf{run_idx, t_slot} = vertcat(temp_objs{~cellfun('isempty', temp_objs)});
        end 
    end

    %% 8. 后处理：在场景-时隙共享的归一化空间中计算指标
    fprintf('\n===== 场景 %d 运行完毕，正在统一归一化并计算指标... =====\n', scenario_id);

    % 原始参考前沿供物理量画图；归一化参考前沿供指标复现。
    scenario_reference_fronts = cell(num_stat_runs, nSlots);
    scenario_normalized_reference_fronts = cell(num_stat_runs, nSlots);
    scenario_normalization_bounds = repmat(struct( ...
        'IdealPoint', [], 'WorstPoint', [], 'ObjectiveRange', [], ...
        'HVReferencePoint', [1.1, 1.1]), nSlots, 1);

    for t_slot = 1:nSlots
        % 1. 汇总当前场景、当前时隙下所有算法和所有独立运行的有效解。
        all_objs_this_scenario_slot = [];
        for r_idx = 1:num_stat_runs
            objs = all_combined_solutions_for_global_pf{r_idx, t_slot};
            if ~isempty(objs)
                valid_mask = ~any(objs >= 1e9 | isnan(objs) | isinf(objs), 2);
                all_objs_this_scenario_slot = [all_objs_this_scenario_slot; objs(valid_mask, :)]; %#ok<AGROW>
            end
        end

        if isempty(all_objs_this_scenario_slot)
            error('M4:NoValidObjectives', ...
                '场景 %d 时隙 %d 没有可用于指标计算的有效目标值。', scenario_id, t_slot);
        end

        % 2. 所有 30 轮共享同一组 ideal/worst 边界。
        ideal_point = min(all_objs_this_scenario_slot, [], 1);
        worst_point = max(all_objs_this_scenario_slot, [], 1);
        objective_range = worst_point - ideal_point;
        objective_range(objective_range < 1e-12) = 1;
        hv_reference_point = [1.1, 1.1];
        num_objectives = numel(ideal_point);

        scenario_normalization_bounds(t_slot).IdealPoint = ideal_point;
        scenario_normalization_bounds(t_slot).WorstPoint = worst_point;
        scenario_normalization_bounds(t_slot).ObjectiveRange = objective_range;
        scenario_normalization_bounds(t_slot).HVReferencePoint = hv_reference_point;

        fprintf('  时隙 %d: ideal=[%.6g %.6g], worst=[%.6g %.6g]\n', ...
            t_slot, ideal_point(1), ideal_point(2), worst_point(1), worst_point(2));

        for r_idx = 1:num_stat_runs
            % 3. 本轮的经验参考前沿仍由本轮全部算法的并集构建。
            all_objs_this_run = all_combined_solutions_for_global_pf{r_idx, t_slot};
            if ~isempty(all_objs_this_run)
                valid_mask = ~any(all_objs_this_run >= 1e9 | isnan(all_objs_this_run) | isinf(all_objs_this_run), 2);
                all_objs_this_run = all_objs_this_run(valid_mask, :);
            end

            raw_reference_front = [];
            if ~isempty(all_objs_this_run)
                ref_idx = FindNonDominatedSolutions(all_objs_this_run);
                raw_reference_front = all_objs_this_run(ref_idx, :);
                [~, sort_idx] = sort(raw_reference_front(:, 1));
                raw_reference_front = raw_reference_front(sort_idx, :);
            end
            % [后处理修复-01 | 2026-09-24]
            % 极端场景中，某次运行可能没有任何算法返回有效前沿。此时普通
            % 空数组为0x0，不能与1x2归一化边界做隐式扩展运算。将其显式
            % 保存为0xM目标矩阵，使该轮IGD/HV保持NaN而不中断全部后处理。
            if isempty(raw_reference_front)
                raw_reference_front = zeros(0, num_objectives);
                normalized_reference_front = zeros(0, num_objectives);
            else
                if size(raw_reference_front, 2) ~= num_objectives
                    error('M4:ReferenceObjectiveDimensionMismatch', ...
                        ['场景 %d 时隙 %d 运行 %d 的参考前沿包含%d个目标，' ...
                         '但归一化边界包含%d个目标。'], ...
                        scenario_id, t_slot, r_idx, ...
                        size(raw_reference_front, 2), num_objectives);
                end
                normalized_reference_front = ...
                    (raw_reference_front - ideal_point) ./ objective_range;
            end

            scenario_reference_fronts{r_idx, t_slot} = raw_reference_front;
            scenario_normalized_reference_fronts{r_idx, t_slot} = normalized_reference_front;

            % 4. IGD/HV/Spacing/Spread 全部在同一归一化空间中计算。
            for alg_idx = 1:length(alg_names_for_results)
                current_alg_name = alg_names_for_results{alg_idx};
                raw_archive = current_scenario_all_run_final_fronts{alg_idx, t_slot, r_idx};
                if ~isempty(raw_archive)
                    valid_mask = ~any(raw_archive >= 1e9 | isnan(raw_archive) | isinf(raw_archive), 2);
                    raw_archive = raw_archive(valid_mask, :);
                end
                % [后处理修复-01 | 2026-09-24]
                % 对单个算法的空前沿采用相同的0xM表示，避免S8等场景下
                % 0x0空数组与1xM边界相减产生“数组大小不兼容”。
                if isempty(raw_archive)
                    raw_archive = zeros(0, num_objectives);
                    normalized_archive = zeros(0, num_objectives);
                else
                    if size(raw_archive, 2) ~= num_objectives
                        error('M4:ArchiveObjectiveDimensionMismatch', ...
                            ['场景 %d 时隙 %d 运行 %d 算法%s的前沿包含%d个目标，' ...
                             '但归一化边界包含%d个目标。'], ...
                            scenario_id, t_slot, r_idx, current_alg_name, ...
                            size(raw_archive, 2), num_objectives);
                    end
                    normalized_archive = ...
                        (raw_archive - ideal_point) ./ objective_range;
                end

                runtime_for_this_run = ...
                    all_scenario_results.(current_alg_name).Runtime{s_idx}(r_idx, t_slot);
                temp_archive_for_metrics = ...
                    createArchiveFromObjectives(normalized_archive, runtime_for_this_run);
                metrics_current_run = CalculateMetricsOnly(temp_archive_for_metrics, ...
                    normalized_reference_front, hv_reference_point);

                all_scenario_results.(current_alg_name).IGD{s_idx}(r_idx, t_slot) = metrics_current_run.IGD;
                all_scenario_results.(current_alg_name).HV{s_idx}(r_idx, t_slot) = metrics_current_run.HV;
                all_scenario_results.(current_alg_name).Spacing{s_idx}(r_idx, t_slot) = metrics_current_run.Spacing;
                all_scenario_results.(current_alg_name).Spread{s_idx}(r_idx, t_slot) = metrics_current_run.Spread;
                all_scenario_results.(current_alg_name).NumFeasibleSolutions{s_idx}(r_idx, t_slot) = metrics_current_run.NumFeasibleSolutions;
            end
        end
    end

    all_reference_fronts{s_idx} = scenario_reference_fronts;
    all_normalized_reference_fronts{s_idx} = scenario_normalized_reference_fronts;
    all_metric_normalization_bounds{s_idx} = scenario_normalization_bounds;

    fprintf('===== 场景 %d 归一化指标计算完成 =====\n', scenario_id);
    clear all_combined_solutions_for_global_pf scenario_reference_fronts ...
        scenario_normalized_reference_fronts scenario_normalization_bounds;
end % 场景循环结束

%% 9. 最终结果展示和保存
fprintf('\n============== 所有场景执行完毕，正在输出最终结果 ===============\n');
statistical_results_to_save = struct(); 
s_idx_to_display_plot = num_scenarios;
t_slot_to_display_plot = nSlots;
current_scenario_display_info = experimental_scenarios{s_idx_to_display_plot};
current_scenario_display_name = sprintf('S%d_I%d_M%d_R2%dx%d', ...
                                        scenario_id, current_scenario_display_info{1}, ...
                                        current_scenario_display_info{2}, ...
                                        current_scenario_display_info{3}, ...
                                        current_scenario_display_info{4});
statistical_results_to_save.(current_scenario_display_name) = struct();

fprintf('\n--- 场景 %d (I=%d, M=%d, R2=%dx%d) 的平均性能 (时隙 %d) ---\n', ...
        scenario_id, current_scenario_display_info{1}, ...
        current_scenario_display_info{2}, current_scenario_display_info{3}, ...
        current_scenario_display_info{4}, t_slot_to_display_plot);
fprintf('------------------------------------------------------------------\n');
slot_field_name_for_save = sprintf('Slot_%d', t_slot_to_display_plot);
statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save) = struct();

for alg_name = alg_names_for_results
    current_alg_name = alg_name{1};
    
    avg_igd = mean(all_scenario_results.(current_alg_name).IGD{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_igd = std(all_scenario_results.(current_alg_name).IGD{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_hv = mean(all_scenario_results.(current_alg_name).HV{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_hv = std(all_scenario_results.(current_alg_name).HV{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_spacing = mean(all_scenario_results.(current_alg_name).Spacing{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_spacing = std(all_scenario_results.(current_alg_name).Spacing{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_spread = mean(all_scenario_results.(current_alg_name).Spread{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_spread = std(all_scenario_results.(current_alg_name).Spread{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_runtime = mean(all_scenario_results.(current_alg_name).Runtime{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_runtime = std(all_scenario_results.(current_alg_name).Runtime{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    avg_tmax = mean(all_scenario_results.(current_alg_name).Tmax{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    std_tmax = std(all_scenario_results.(current_alg_name).Tmax{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    avg_num_feasible = mean(all_scenario_results.(current_alg_name).NumFeasibleSolutions{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    fprintf('    %s:\n', current_alg_name);
    fprintf('      IGD: %.4f (%.4f)\n', avg_igd, std_igd);
    fprintf('      HV: %.4f (%.4f)\n', avg_hv, std_hv);
    fprintf('      Spacing: %.4f (%.4f)\n', avg_spacing, std_spacing);
    fprintf('      Spread: %.4f (%.4f)\n', avg_spread, std_spread);
    fprintf('      Runtime: %.4f (%.4f) s\n', avg_runtime, std_runtime);
    fprintf('      Avg Nondominated Solutions: %.1f\n', avg_num_feasible);
    fprintf('      Min Tmax (Makespan): %.4f (%.4f)\n', avg_tmax, std_tmax);
    
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).IGD_avg = avg_igd;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).IGD_std = std_igd;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).HV_avg = avg_hv;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).HV_std = std_hv;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spacing_avg = avg_spacing;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spacing_std = std_spacing;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spread_avg = avg_spread;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Spread_std = std_spread;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Runtime_avg = avg_runtime;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Runtime_std = std_runtime;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).NumFeasibleSolutions_avg = avg_num_feasible;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Tmax_avg = avg_tmax;
    statistical_results_to_save.(current_scenario_display_name).(slot_field_name_for_save).(current_alg_name).Tmax_std = std_tmax;
end

fprintf('  ----------------------------------------\n');
output_filename = fullfile(output_timestamp_folder, sprintf('Experiment_Statistical_Results_%s_Slot%d_%s.mat', current_scenario_display_name, t_slot_to_display_plot, string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));

save(output_filename, ...
    'all_scenario_results', ...                
    'current_scenario_all_run_final_fronts', ... 
    'all_reference_fronts', ...                % [已修改] 替换为所有场景、所有运行的独立参考前沿矩阵合集
    'all_normalized_reference_fronts', ...
    'all_metric_normalization_bounds', ...
    'statistical_results_to_save', ...         
    'alg_names_for_results', ...
    'metric_names', ...
    'all_pruned_exact_pfs', ...
    'all_REMODDP_training_info', ...
    'params_drl', ...
    'scenario_id', ...
    'scenario_catalog', ...
    'experimental_scenarios', ...
    'nSlots', ...
    'num_stat_runs');

fprintf('\n所有统计结果已保存到文件: %s\n', output_filename);
fprintf('\n============== 实验流程正常结束 ===============\n');

%% 辅助函数
function first_front_archive = getFirstFront(fronts_cell_array)
    if ~isempty(fronts_cell_array) && ~isempty(fronts_cell_array{1})
        first_front_archive = fronts_cell_array{1};
    else
        first_front_archive = [];
    end
end

function obj_matrix = getObjectivesMatrix(archive_struct)
    if ~isempty(archive_struct)
        obj_matrix = vertcat(archive_struct.Objectives);
        valid_mask = ~any(obj_matrix >= 1e9 | isnan(obj_matrix) | isinf(obj_matrix), 2);
        obj_matrix = obj_matrix(valid_mask, :);
    else
        obj_matrix = [];
    end
end

function archive_struct_out = createArchiveFromObjectives(obj_matrix, runtime_val)
    if isempty(obj_matrix)
        archive_struct_out = [];
        return;
    end
    archive_struct_out = repmat(struct('Objectives', [], 'RunTime', NaN), size(obj_matrix, 1), 1);
    for i = 1:size(obj_matrix, 1)
        archive_struct_out(i).Objectives = obj_matrix(i,:);
        archive_struct_out(i).RunTime = runtime_val;
    end
end

function problem = CreateFogProblem(scenario_config, problem_base)
    problem = problem_base;
    problem.nTerminals = scenario_config{1};
    problem.nFogNodes = scenario_config{2};
    problem.area = [0 scenario_config{3}; 0 scenario_config{4}];
    
    term_x_coords = problem.area(1,1) + (problem.area(1,2) - problem.area(1,1)) .* rand(problem.nTerminals, 1);
    term_y_coords = problem.area(2,1) + (problem.area(2,2) - problem.area(2,1)) .* rand(problem.nTerminals, 1);
    problem.terminalProperties.positions = [term_x_coords, term_y_coords];
    problem.terminalProperties.Pt_dbm = linspace(10, 15, problem.nTerminals); 
    problem.terminalProperties.fc = linspace(2.4e9, 5.8e9, problem.nTerminals); 
    
    problem.fogNodeProperties.cpu_cycle_rate = linspace(2e9, 5e9, problem.nFogNodes); 
    x_coords = problem.area(1,1) + (problem.area(1,2) - problem.area(1,1)) .* rand(problem.nFogNodes, 1);
    y_coords = problem.area(2,1) + (problem.area(2,2) - problem.area(2,1)) .* rand(problem.nFogNodes, 1);
    problem.initial_fog_positions_matrix = [x_coords, y_coords];
    problem.initial_fog_deployment_flat = reshape(problem.initial_fog_positions_matrix', 1, []);
    
    problem.terminalProperties.task_sizes = 0.1e6 + (1e6 - 0.1e6) .* rand(1, problem.nTerminals);
    problem.bounds.bandwidth = [ones(1, problem.nTerminals)*0.2e6; ones(1, problem.nTerminals)*10e6];
end
