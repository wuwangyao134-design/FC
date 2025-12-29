% % =========================================================================
% % Tmax 分析与提取脚本 (Tmax_Analysis.m)
% % 描述: 
% %   本脚本用于处理已保存的单次运行实验结果。
% %   它会自动加载每个场景、每个算法的 .mat 文件，
% %   从 Pareto 最优解集中找到拐点 (Knee Point)，
% %   并提取该拐点解所对应的 Tmax 值。
% %   最后，将所有结果汇总到一个表格中。
% %
% % 使用前请确认:
% %   1. `base_path` 变量已设置为您存放 .mat 文件的正确路径。
% %   2. `scenario_configs` 中的文件名与您磁盘上的实际文件名相符。
% % =========================================================================
% 
% clc;
% clear;
% close all;
% 
% %% 1. 配置区域 (请根据您的实际情况修改)
% 
% % --- 数据文件的根目录 ---
% base_path = 'E:\guthub-matlab\第二篇论文数据'; 
% 
% % --- 场景配置 (更新为更稳健的格式) ---
% % 格式: { '用于显示的场景名', '您磁盘上实际的文件名.mat' }
% scenario_configs = {
%     'S1_20_2_100x100',   '30-5-20-2-100-100.mat';
%     'S2_20_2_200x200',   '30-5-20-200-200.mat';   % 已根据您的描述修正文件名
%     'S3_45_3_100x100',   '30-5-45-100-100.mat';   % 已根据您的描述修正文件名
%     'S4_45_3_200x200',   '30-5-45-200-200.mat';   % 已根据您的描述修正文件名
%     'S5_80_4_100x100',   '30-5-80-100-100.mat';   % 已根据您的描述修正文件名
%     'S6_80_4_200x200',   '30-5-80-200-200.mat'    % 已根据您的描述修正文件名
% };
% 
% % --- 算法配置 ---
% % 这个结构体将图例名称和 .mat 文件中变量名的特征字符串关联起来
% alg_configs = {
%     'Ours',      'FinalArchive_MyNSG...';
%     'NSGA-II',   'FinalArchive_NSGA_Sl...';
%     'A-MOPSO',   'FinalArchive_MOPSO_...';
%     'ARG',       'FinalArchive_Baseline...'
% };
% 
% %% 2. 初始化结果存储
% num_scenarios = size(scenario_configs, 1);
% num_algs = size(alg_configs, 1);
% 
% % --- 从配置中直接获取场景名称 (更新) ---
% scenario_names = scenario_configs(:, 1); % 直接取第一列作为行名
% 
% % 使用 table 格式存储结果，更清晰
% alg_legend_names = alg_configs(:, 1)';
% results_table = array2table(NaN(num_scenarios, num_algs), ...
%     'RowNames', scenario_names, 'VariableNames', alg_legend_names);
% 
% %% 3. 主循环：遍历所有场景和算法
% 
% fprintf('开始处理数据...\n');
% 
% % --- 遍历每个场景 ---
% for i = 1:num_scenarios
%     scenario_name = scenario_configs{i, 1}; % 从配置中获取场景名
%     file_name = scenario_configs{i, 2};     % 从配置中获取文件名
%     file_path = fullfile(base_path, file_name); % 拼接完整路径
% 
%     fprintf('\n--- 正在处理场景: %s ---\n', scenario_name);
% 
%     if ~exist(file_path, 'file')
%         warning('文件未找到，跳过场景: %s', file_path);
%         continue;
%     end
% 
%     % 加载该场景的 .mat 文件
%     % 文件中的变量会被加载到一个 struct 中，避免污染当前工作区
%     data = load(file_path);
% 
%     % --- 遍历每个算法 ---
%     for j = 1:num_algs
%         alg_legend_name = alg_configs{j, 1};
%         alg_var_pattern = alg_configs{j, 2};
% 
%         % 从加载的数据中找到对应的 FinalArchive 变量
%         % 使用 startsWith 来模糊匹配变量名
%         all_vars = fieldnames(data);
%         target_var_name = '';
%         for k = 1:length(all_vars)
%             if startsWith(all_vars{k}, strrep(alg_var_pattern, '...', ''))
%                 target_var_name = all_vars{k};
%                 break;
%             end
%         end
% 
%         if isempty(target_var_name)
%             fprintf('  在文件 %s 中未找到算法 [%s] 的数据，跳过。\n', file_name, alg_legend_name);
%             continue;
%         end
% 
%         final_archive = data.(target_var_name);
% 
%         if isempty(final_archive)
%             fprintf('  算法 [%s] 的解集为空，跳过。\n', alg_legend_name);
%             continue;
%         end
% 
%         % 核心步骤：找到拐点并提取 Tmax
%         try
%             knee_point_solution = findKneePoint(final_archive);
%             tmax_value = knee_point_solution.Tmax;
% 
%             % 将结果存入表格
%             results_table{scenario_name, alg_legend_name} = tmax_value;
%             fprintf('  算法 [%-10s]: 找到拐点 Tmax = %.4f\n', alg_legend_name, tmax_value);
% 
%         catch ME
%             warning('在处理算法 [%s] 时出错: %s', alg_legend_name, ME.message);
%         end
%     end
% end
% 
% %% 4. 显示最终结果
% 
% fprintf('\n\n================================================================\n');
% fprintf('                Tmax 性能比较 (拐点解)\n');
% fprintf('================================================================\n');
% disp(results_table);
% 
% 
% %% 5. 辅助函数：寻找拐点 (Knee Point)
% 
% function knee_point_solution = findKneePoint(archive)
%     % 描述:
%     %   从一个 Pareto 最优解集 (archive) 中找到拐点解。
%     %   方法：将目标函数值归一化后，找到距离理想点 [0, 0] 最近的点。
%     %
%     % 输入:
%     %   archive - 结构体数组，每个元素包含 .Objectives 字段
%     %
%     % 输出:
%     %   knee_point_solution - 拐点对应的那个完整的解结构体
% 
%     if isempty(archive) || ~isfield(archive, 'Objectives')
%         error('输入的解集 (archive) 为空或缺少 Objectives 字段。');
%     end
% 
%     % 步骤 1: 提取所有目标函数值
%     objectives = vertcat(archive.Objectives);
% 
%     if size(objectives, 1) < 2
%         % 如果只有一个解，它自身就是拐点
%         knee_point_solution = archive(1);
%         return;
%     end
% 
%     % 步骤 2: 将目标函数值归一化到 [0, 1] 区间
%     min_vals = min(objectives, [], 1);
%     max_vals = max(objectives, [], 1);
% 
%     % 处理分母为零的情况
%     range_vals = max_vals - min_vals;
%     range_vals(range_vals == 0) = 1;
% 
%     normalized_objectives = (objectives - min_vals) ./ range_vals;
% 
%     % 步骤 3: 计算每个点到理想点 [0, 0] 的欧几里得距离
%     % 对于最小化问题，理想点是 [0, 0]
%     distances = sqrt(sum(normalized_objectives.^2, 2));
% 
%     % 步骤 4: 找到距离最小的点的索引
%     [~, knee_point_index] = min(distances);
% 
%     % 步骤 5: 返回对应的完整解
%     knee_point_solution = archive(knee_point_index);
% end

% =========================================================================
% 最高 Tmax 分析脚本 (Tmax_Max_Analysis.m)
% 描述: 
%   本脚本分析每个算法 Pareto 最优解集中，Tmax 值最高的那个解。
%   这用于回答“哪个算法找到的最优解集，其最差情况下的延迟表现如何？”
%   这个问题，反映了解集的稳定性和风险。
% =========================================================================

clc;
clear;
close all;

%% 1. 配置区域 (与之前脚本相同)

% --- 数据文件的根目录 ---
base_path = 'E:\guthub-matlab\第二篇论文数据'; 

% --- 场景配置 ---
scenario_configs = {
    'S1_20_2_100x100',   '30-5-20-2-100-100.mat';
    % 'S2_20_2_200x200',   '30-5-20-200-200.mat';
    % 'S3_45_3_100x100',   '30-5-45-100-100.mat';
    % 'S4_45_3_200x200',   '30-5-45-200-200.mat';
    % 'S5_80_4_100x100',   '30-5-80-100-100.mat';
    % 'S6_80_4_200x200',   '30-5-80-200-200.mat'
};

% --- 算法配置 ---
alg_configs = {
    'Ours',      'FinalArchive_MyNSG...';
    'NSGA-II',   'FinalArchive_NSGA_Sl...';
    'A-MOPSO',   'FinalArchive_MOPSO_...';
    'ARG',       'FinalArchive_Baseline...'
};

%% 2. 初始化结果存储
num_scenarios = size(scenario_configs, 1);
num_algs = size(alg_configs, 1);
scenario_names = scenario_configs(:, 1);
alg_legend_names = alg_configs(:, 1)';
results_table_max_tmax = array2table(NaN(num_scenarios, num_algs), ...
    'RowNames', scenario_names, 'VariableNames', alg_legend_names);

%% 3. 主循环：遍历所有场景和算法

fprintf('开始处理数据 (寻找最高 Tmax)...\n');

for i = 1:num_scenarios
    scenario_name = scenario_configs{i, 1};
    file_name = scenario_configs{i, 2};
    file_path = fullfile(base_path, file_name);
    
    fprintf('\n--- 正在处理场景: %s ---\n', scenario_name);
    
    if ~exist(file_path, 'file')
        warning('文件未找到，跳过场景: %s', file_path);
        continue;
    end
    
    data = load(file_path);
    
    for j = 1:num_algs
        alg_legend_name = alg_configs{j, 1};
        alg_var_pattern = alg_configs{j, 2};
        
        all_vars = fieldnames(data);
        target_var_name = '';
        for k = 1:length(all_vars)
            if startsWith(all_vars{k}, strrep(alg_var_pattern, '...', ''))
                target_var_name = all_vars{k};
                break;
            end
        end
        
        if isempty(target_var_name)
            fprintf('  在文件 %s 中未找到算法 [%s] 的数据，跳过。\n', file_name, alg_legend_name);
            continue;
        end
        
        final_archive = data.(target_var_name);
        
        if isempty(final_archive)
            fprintf('  算法 [%s] 的解集为空，跳过。\n', alg_legend_name);
            continue;
        end
        
        % --- 核心步骤 (已修改) ---
        % 不再找最小值，而是直接从所有解中找到 Tmax 最大的那个
        try
            all_tmax_values = [final_archive.Tmax];
            max_tmax_value = max(all_tmax_values); % <--- 修改为 max()
            
            % 将结果存入表格
            results_table_max_tmax{scenario_name, alg_legend_name} = max_tmax_value;
            fprintf('  算法 [%-10s]: 找到最高 Tmax = %.4f\n', alg_legend_name, max_tmax_value);
            
        catch ME
            warning('在处理算法 [%s] 时出错: %s', alg_legend_name, ME.message);
        end
    end
end

%% 4. 显示最终结果

fprintf('\n\n================================================================\n');
fprintf('           最高 Tmax 性能比较 (遍历整个前沿)\n');
fprintf('================================================================\n');
disp(results_table_max_tmax);
