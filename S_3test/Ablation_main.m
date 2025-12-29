%% =========================================================================
% HMD-NSGAII 消融实验主执行脚本 - 严谨科研版
% =========================================================================
clc; clear; close all;

%% 0. 全局参数设定
nSlots = 10; 
num_stat_runs = 30; % 正式实验请务必保持 30 次以满足统计学显著性
shadow_std_dev = struct('LoS', 3.0, 'NLoS', 8.29); 

% 设定测试场景 S3 (I=45, M=3)
experimental_scenarios = {{80, 4, 100, 100}};
alg_names = {'HMD_Full', 'HMD_no_ISMM', 'HMD_no_AMS', 'HMD_no_Hybrid', 'Standard_NSGAII'};
metric_names = {'IGD', 'HV', 'Spacing', 'NumFeasibleSolutions', 'Runtime'}; 

% 初始化结果存储
all_results = struct();
for a = alg_names 
    for m = metric_names
        all_results.(a{1}).(m{1}) = zeros(num_stat_runs, nSlots);
    end 
end

%% 1. 场景与统计循环
problem_base.objFunc = @EvaluateParticle;
problem_base.Tslot = 5;
problem_base.systemTotalBandwidth = 225e6;
problem_base.nObj = 2; 

for s_idx = 1:length(experimental_scenarios)
    config = experimental_scenarios{s_idx};
    problem = problem_base;
    problem.nTerminals = config{1};
    problem.nFogNodes = config{2};
    problem.area = [0 config{3}; 0 config{4}];
    
    % 算法基础参数
    params_base = struct('N', 100, 'T_max', 200, 'pc', 0.9, 'pm', 0.05, 'mu', 20, 'mum', 20);
    params_base.pm_min = 0.01; params_base.pm_max = 0.1;
    params_base.mum_min = 5; params_base.mum_max = 10;

    current_run_archives = cell(length(alg_names), nSlots, num_stat_runs);
    all_objs_for_pf_star = cell(num_stat_runs, nSlots);

    for run_idx = 1:num_stat_runs
        rng(run_idx + 100); 
        % --- 环境初始化 (保持各变体在同一环境下竞争) ---
        problem.terminalProperties.positions = [problem.area(1,2)*rand(problem.nTerminals,1), problem.area(2,2)*rand(problem.nTerminals,1)];
        problem.terminalProperties.Pt_dbm = linspace(10, 15, problem.nTerminals); 
        problem.terminalProperties.fc = linspace(2.4e9, 5.8e9, problem.nTerminals); 
        problem.fogNodeProperties.cpu_cycle_rate = linspace(2e9, 5e9, problem.nFogNodes); 
        x_c = problem.area(1,2)*rand(problem.nFogNodes, 1); y_c = problem.area(2,2)*rand(problem.nFogNodes, 1);
        problem.initial_fog_positions_matrix = [x_c, y_c];
        problem.initial_fog_deployment_flat = reshape(problem.initial_fog_positions_matrix', 1, []);
        problem.terminalProperties.task_sizes = (0.1e6 + 0.9e6*rand(1, problem.nTerminals));
        problem.bounds.bandwidth = [ones(1,problem.nTerminals)*0.2e6; ones(1,problem.nTerminals)*10e6];
        problem.fixed_shadow_LoS_val = shadow_std_dev.LoS * randn(1, problem.nTerminals);
        problem.fixed_shadow_NLoS_val = shadow_std_dev.NLoS * randn(1, problem.nTerminals);
        
        LastSlotArchives = cell(1, length(alg_names)); 

        for t = 1:nSlots
            fprintf('Run %d, Slot %d\n', run_idx, t);
            combined_objs_this_slot = []; 

            for a_idx = 1:length(alg_names)
                name = alg_names{a_idx};
                p = params_base; mem = LastSlotArchives{a_idx};
                % 消融实验开关
                switch name
                    case 'HMD_Full',      p.memory_ratio=0.1; p.adaptive_enabled=true;  p.hybrid_enabled=true;
                    case 'HMD_no_ISMM',   p.memory_ratio=0;   p.adaptive_enabled=true;  p.hybrid_enabled=true; mem=[];
                    case 'HMD_no_AMS',    p.memory_ratio=0.1; p.adaptive_enabled=false; p.hybrid_enabled=true;
                    case 'HMD_no_Hybrid', p.memory_ratio=0.1; p.adaptive_enabled=true;  p.hybrid_enabled=false;
                    case 'Standard_NSGAII', p.memory_ratio=0; p.adaptive_enabled=false; p.hybrid_enabled=false; mem=[];
                end
                
                tic;
                Pop = OUSNSGA_II(problem, p, mem); 
                runtime = toc;
                
                Archive = getFirstFront(FindAllFronts(Pop));
                LastSlotArchives{a_idx} = Archive;
                
                objs = getObjectivesMatrix(Archive);
                current_run_archives{a_idx, t, run_idx} = objs;
                combined_objs_this_slot = [combined_objs_this_slot; objs];
                all_results.(name).Runtime(run_idx, t) = runtime;
            end
            all_objs_for_pf_star{run_idx, t} = combined_objs_this_slot;
        end
    end

    %% 2. 指标重计算阶段 (使用独立函数)
    fprintf('\n正在进行指标重计算与归一化分析...\n');
    for t = 1:nSlots
        % 聚合所有算法的解来构建该时隙的 PF*
        objs_all = vertcat(all_objs_for_pf_star{:, t});
        valid_mask = ~any(objs_all >= 1e9 | isnan(objs_all), 2);
        objs_all = objs_all(valid_mask, :);
        
        if ~isempty(objs_all)
            pf_idx = FindNonDominatedSolutions(objs_all);
            pf_star = objs_all(pf_idx, :);
            % 动态参考点：防止目标量级不同导致 HV/IGD 偏移
            hv_ref = max(objs_all, [], 1) * 1.1; 
        else
            pf_star = []; hv_ref = [1.1, 1.1];
        end
        
        for a_idx = 1:length(alg_names)
            name = alg_names{a_idx};
            for r = 1:num_stat_runs
                objs = current_run_archives{a_idx, t, r};
                if isempty(objs)
                    all_results.(name).IGD(r,t)=1; all_results.(name).HV(r,t)=0; all_results.(name).Spacing(r,t)=1;
                    continue;
                end
                
                % 调用你的独立函数：calculateHV
                all_results.(name).HV(r, t) = calculateHV(objs, hv_ref);
                
                % 调用你的独立函数：calculateIGD (传入归一化数据以确保公平)
                norm_objs = objs ./ hv_ref;
                norm_pf = pf_star ./ hv_ref;
                all_results.(name).IGD(r, t) = calculateIGD(norm_objs, norm_pf);
                
                % Spacing 计算 (衡量分布均匀性)
                if size(objs, 1) > 1
                    d_mat = pdist2(objs, objs); d_mat(d_mat==0) = inf;
                    di = min(d_mat, [], 2);
                    all_results.(name).Spacing(r, t) = sqrt(sum((mean(di)-di).^2)/(size(objs,1)-1));
                else
                    all_results.(name).Spacing(r, t) = 1;
                end
                all_results.(name).NumFeasibleSolutions(r, t) = size(objs, 1);
            end
        end
    end
end

%% 3. 统计结果展示
fprintf('\n%s\n', repmat('=', 1, 120));
fprintf('%-20s | %-18s | %-18s | %-18s | %-12s\n', ...
    'Variant', 'IGD (Mean±Std)', 'HV (Mean±Std)', 'Spacing (Mean±Std)', 'Feasibility');
fprintf('%s\n', repmat('-', 1, 120));
for i = 1:length(alg_names)
    name = alg_names{i};
    % 以最后一个时隙 (Slot 10) 作为最终性能对比
    fIGD = all_results.(name).IGD(:, end); fHV = all_results.(name).HV(:, end);
    fSP = all_results.(name).Spacing(:, end); fFeas = all_results.(name).NumFeasibleSolutions(:, end);
    % --- 找到这一行进行修改 ---
    fprintf('%-20s | %.4f ± %.4f | %.4f ± %.4f | %.4f ± %.4f | %.2f%%\n', ...
        name, mean(fIGD), std(fIGD), mean(fHV), std(fHV), mean(fSP), std(fSP), ...
        (mean(fFeas)/100)*100); % 直接除以 100 (种群大小 N)
end
fprintf('%s\n', repmat('=', 1, 120));

%% ========================== 脚本专用辅助函数 ==========================

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