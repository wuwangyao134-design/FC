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
num_stat_runs = 30; % 统计运行次数
shadow_std_dev.LoS = 3.0;  
shadow_std_dev.NLoS = 8.29; 

% --- 定义实验场景参数 ---
experimental_scenarios = {
    % {20, 2, 100, 100};  % S1
    % {20, 2, 200, 200};  % S2
    % {45, 3, 100, 100};  % S3
    % {45, 3, 200, 200};  % S4
    % {80, 4, 100, 100};  % S5
     {80, 4, 200, 200};  % S6
};
num_scenarios = length(experimental_scenarios); 

% --- [修改] 更新算法名称列表 ---
alg_names_for_results = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III'};
metric_names = {'IGD', 'HV', 'Spacing', 'Spread', 'Runtime', 'NumFeasibleSolutions'}; 

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

%% 5. 场景循环
problem_base.objFunc = @EvaluateParticle;
problem_base.Tslot = 5;
problem_base.systemTotalBandwidth = 225e6;
problem_base.nObj = 2; 

for s_idx = 1:num_scenarios
    problem = problem_base; 
    current_scenario_config = experimental_scenarios{s_idx};
    problem.nTerminals = current_scenario_config{1};
    problem.nFogNodes = current_scenario_config{2};
    problem.area = [0 current_scenario_config{3}; 0 current_scenario_config{4}];
            
    % --- [核心：基础参数] ---
    % 定义一个通用的基础结构体，供所有对比算法继承
    base_params = struct('N', 100, 'T_max', 200, 'pc', 0.9, 'pm', 0.05, 'mu', 20, 'mum', 20, ...
                         'pm_cont_coeff', 0.8, 'pm_disc_coeff', 1.2, ...
                         'mum_cont_coeff', 1.1, 'mum_disc_coeff', 0.9);
        
    % --- [1. MyNSGA-II (OUSNSGA_II) 配置] ---
    params_my_nsga = base_params; 
    params_my_nsga.adaptive_enabled = true; 
    params_my_nsga.hybrid_enabled   = true; 
    params_my_nsga.memory_ratio     = 0.1;  
    params_my_nsga.pm_min = 0.01; 
    params_my_nsga.pm_max = 0.1;
    params_my_nsga.mum_min = 5; 
    params_my_nsga.mum_max = 20;
        
    % --- [2. 其他对比算法继承基础参数] ---
    params_nsga = base_params; % 原始 NSGA-II
    
    params_moead = base_params; % MOEA/D
    params_moead.T = 20;             
    params_moead.delta = 0.9;        
    
    params_insga = base_params; % INSGA-II
    
    params_nsga3 = base_params; % NSGA-III
    params_nsga3.p = 99;

    current_scenario_all_run_final_fronts = cell(length(alg_names_for_results), nSlots, num_stat_runs);
    all_combined_solutions_for_global_pf = cell(num_stat_runs, nSlots);

    %% 6. 统计运行循环
    for run_idx = 1:num_stat_runs
        seed_value = (s_idx - 1) * 1000 + run_idx;
        rng(seed_value);
        
        % 问题实例初始化
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
        problem.fixed_shadow_LoS_val = shadow_std_dev.LoS * randn(1, problem.nTerminals);
        problem.fixed_shadow_NLoS_val = shadow_std_dev.NLoS * randn(1, problem.nTerminals);
        
        LastSlotArchive_MyNSGA = []; 
        
        %% 7. 时隙循环
        for t_slot = 1:nSlots
            fprintf('场景%d-运行%d: 执行时隙 %d/%d\n', s_idx, run_idx, t_slot, nSlots);
            
            % 1. MyNSGA-II
            tic;
            Pop_MyNSGA = OUSNSGA_II(problem, params_my_nsga, LastSlotArchive_MyNSGA);
            time_my_nsga = toc;
            all_scenario_results.MyNSGA_II.Runtime{s_idx}(run_idx, t_slot) = time_my_nsga;
            Archive_MyNSGA = getFirstFront(FindAllFronts(Pop_MyNSGA));
            LastSlotArchive_MyNSGA = Archive_MyNSGA;
            
            % 2. NSGA-II
            tic;
            Pop_NSGA2 = DNSGA_II(problem, params_nsga, []);
            time_nsga2 = toc;
            all_scenario_results.NSGA_II.Runtime{s_idx}(run_idx, t_slot) = time_nsga2;
            Archive_NSGA2 = getFirstFront(FindAllFronts(Pop_NSGA2));
            
            % 3. [修改] 运行 MOEAD 
            tic;
            Pop_MOEAD = MOEAD(problem, params_moead, []);
            time_moead = toc;
            all_scenario_results.MOEA_D.Runtime{s_idx}(run_idx, t_slot) = time_moead;
            % MOEA/D 产生的种群即为各个权重下的最优解集合
            Archive_MOEAD = Pop_MOEAD;
            
            % 4. [修改] 运行 INSGA-II 
            tic;
            Pop_INSGA = INSGA_II(problem, params_insga, []);
            time_insga = toc;
            all_scenario_results.INSGA_II.Runtime{s_idx}(run_idx, t_slot) = time_insga;
            Archive_INSGA = getFirstFront(FindAllFronts(Pop_INSGA));
            
            % 5. 运行 NSGA-III
            tic;
            Pop_NSGA3 = NSGA_III(problem, params_nsga3, []);
            time_nsga3 = toc;
            all_scenario_results.NSGA_III.Runtime{s_idx}(run_idx, t_slot) = time_nsga3;
            Archive_NSGA3 = getFirstFront(FindAllFronts(Pop_NSGA3));
            
            % --- 收集前沿数据用于指标计算 ---
            current_scenario_all_run_final_fronts{1, t_slot, run_idx} = getObjectivesMatrix(Archive_MyNSGA);
            current_scenario_all_run_final_fronts{2, t_slot, run_idx} = getObjectivesMatrix(Archive_NSGA2);
            current_scenario_all_run_final_fronts{3, t_slot, run_idx} = getObjectivesMatrix(Archive_MOEAD);
            current_scenario_all_run_final_fronts{4, t_slot, run_idx} = getObjectivesMatrix(Archive_INSGA);
            current_scenario_all_run_final_fronts{5, t_slot, run_idx} = getObjectivesMatrix(Archive_NSGA3); 
            
            % 收集所有解构建全局 PF*
            temp_objs = {getObjectivesMatrix(Archive_MyNSGA), getObjectivesMatrix(Archive_NSGA2), ...
             getObjectivesMatrix(Archive_MOEAD), getObjectivesMatrix(Archive_INSGA), ...
             getObjectivesMatrix(Archive_NSGA3)};
            all_combined_solutions_for_global_pf{run_idx, t_slot} = vertcat(temp_objs{~cellfun('isempty', temp_objs)});
            
        end 
    end
    %% 8. 后处理：构建 PF* 并计算指标
    fprintf('\n===== 场景 %d 运行完毕，正在计算指标... =====\n', s_idx);
    
    for t_slot = 1:nSlots
        fprintf('  处理时隙 %d 的指标计算...\n', t_slot);
    
        % 1. 构建近似参考前沿 (PF*)
        temp_pf_collector = all_combined_solutions_for_global_pf(:, t_slot);
        all_objs_for_pf_star = vertcat(temp_pf_collector{:});
        
        if ~isempty(all_objs_for_pf_star)
            % 过滤无效解
            valid_objs_mask = ~any(all_objs_for_pf_star >= 1e9 | isnan(all_objs_for_pf_star) | isinf(all_objs_for_pf_star), 2);
            all_objs_for_pf_star = all_objs_for_pf_star(valid_objs_mask, :);
        end
        
        % 计算 IGD 参考前沿
        igd_reference_front_obj = [];
        if size(all_objs_for_pf_star, 1) > 1
            igd_ref_front_idx = FindNonDominatedSolutions(all_objs_for_pf_star); 
            igd_reference_front_obj = all_objs_for_pf_star(igd_ref_front_idx, :);
            [~, sort_idx_igd] = sort(igd_reference_front_obj(:,1));
            igd_reference_front_obj = igd_reference_front_obj(sort_idx_igd, :);
        end
        
        % 计算 HV 参考点 (1.1倍最大值)
        hv_reference_point = [];
        if ~isempty(all_objs_for_pf_star)
            hv_reference_point = max(all_objs_for_pf_star, [], 1) * 1.1; 
        end
        if isempty(hv_reference_point), hv_reference_point = [0.8, 0.35]; end

        % --- 2. 核心修改：遍历更新后的算法列表计算指标 ---
        for alg_idx = 1:length(alg_names_for_results)
            current_alg_name = alg_names_for_results{alg_idx}; % 例如 'MOPO', 'INSGA_II'
            
            for r_idx = 1:num_stat_runs
                % 从保存的缓存中提取当前算法、当前时隙、当前运行的前沿矩阵
                current_run_archive_obj_matrix = current_scenario_all_run_final_fronts{alg_idx, t_slot, r_idx}; 
                
                % 提取对应算法的运行时间
                runtime_for_this_run = all_scenario_results.(current_alg_name).Runtime{s_idx}(r_idx, t_slot);
                
                % 封装为临时存档结构以兼容指标计算函数
                temp_archive_for_metrics = createArchiveFromObjectives(current_run_archive_obj_matrix, runtime_for_this_run);
                
                % 调用计算函数
                metrics_current_run = CalculateMetricsOnly(temp_archive_for_metrics, ...
                    igd_reference_front_obj, hv_reference_point);
                
                % 将结果填回 all_scenario_results 结构体
                all_scenario_results.(current_alg_name).IGD{s_idx}(r_idx, t_slot) = metrics_current_run.IGD;
                all_scenario_results.(current_alg_name).HV{s_idx}(r_idx, t_slot) = metrics_current_run.HV;
                all_scenario_results.(current_alg_name).Spacing{s_idx}(r_idx, t_slot) = metrics_current_run.Spacing;
                all_scenario_results.(current_alg_name).Spread{s_idx}(r_idx, t_slot) = metrics_current_run.Spread;
                all_scenario_results.(current_alg_name).NumFeasibleSolutions{s_idx}(r_idx, t_slot) = metrics_current_run.NumFeasibleSolutions;
            end
        end
    end
    fprintf('\n===== 场景 %d 指标计算完成 =====\n', s_idx);
    
    clear all_combined_solutions_for_global_pf; 
end % 场景循环结束

%% 9. 最终结果展示和保存
fprintf('\n============== 所有场景执行完毕，正在输出最终结果 ===============\n');
statistical_results_to_save = struct(); 
s_idx_to_display_plot = num_scenarios;
t_slot_to_display_plot = nSlots;
current_scenario_display_info = experimental_scenarios{s_idx_to_display_plot};
current_scenario_display_name = sprintf('S%d_I%d_M%d_R2%dx%d', ...
                                        s_idx_to_display_plot, current_scenario_display_info{1}, ...
                                        current_scenario_display_info{2}, ...
                                        current_scenario_display_info{3}, ...
                                        current_scenario_display_info{4});
statistical_results_to_save.(current_scenario_display_name) = struct();
fprintf('\n--- 场景 %d (I=%d, M=%d, R2=%dx%d) 的平均性能 (时隙 %d) ---\n', ...
        s_idx_to_display_plot, current_scenario_display_info{1}, ...
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
    avg_num_feasible = mean(all_scenario_results.(current_alg_name).NumFeasibleSolutions{s_idx_to_display_plot}(:, t_slot_to_display_plot), 'omitnan');
    
    fprintf('    %s:\n', current_alg_name);
    fprintf('      IGD: %.4f (%.4f)\n', avg_igd, std_igd);
    fprintf('      HV: %.4f (%.4f)\n', avg_hv, std_hv);
    fprintf('      Spacing: %.4f (%.4f)\n', avg_spacing, std_spacing);
    fprintf('      Spread: %.4f (%.4f)\n', avg_spread, std_spread);
    fprintf('      Runtime: %.4f (%.4f) s\n', avg_runtime, std_runtime);
    fprintf('      Avg Feasible Solutions: %.1f\n', avg_num_feasible);

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
end
fprintf('  ----------------------------------------\n');

output_filename = fullfile(output_timestamp_folder, sprintf('Experiment_Statistical_Results_%s_Slot%d_%s.mat', current_scenario_display_name, t_slot_to_display_plot, string(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
save(output_filename, ...
    'all_scenario_results', ...                % 包含 30次运行 x 10时隙 的原始指标数据
    'current_scenario_all_run_final_fronts', ... % 包含所有运行的前沿坐标数据 (画散点图用)
    'igd_reference_front_obj', ...             % 包含参考前沿数据
    'statistical_results_to_save', ...         % 包含汇总统计表
    'alg_names_for_results', ...
    'metric_names', ...
    'experimental_scenarios', ...
    'nSlots', ...
    'num_stat_runs');
fprintf('\n所有统计结果已保存到文件: %s\n', output_filename);
fprintf('\n============== 实验流程正常结束 ===============\n');

%% 辅助函数
% --- 辅助函数：提取第一非支配前沿 ---
function first_front_archive = getFirstFront(fronts_cell_array)
    if ~isempty(fronts_cell_array) && ~isempty(fronts_cell_array{1})
        first_front_archive = fronts_cell_array{1};
    else
        first_front_archive = [];
    end
end

% --- 辅助函数：从存档结构体中提取目标值矩阵 ---
function obj_matrix = getObjectivesMatrix(archive_struct)
    if ~isempty(archive_struct)
        obj_matrix = vertcat(archive_struct.Objectives);
        % --- [已修改] 使用更健壮的、通用的过滤方法 ---
        valid_mask = ~any(obj_matrix >= 1e9 | isnan(obj_matrix) | isinf(obj_matrix), 2);
        obj_matrix = obj_matrix(valid_mask, :);
    else
        obj_matrix = [];
    end
end

% --- 辅助函数：创建Archive结构体以兼容CalculateMetricsOnly的输入 ---
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