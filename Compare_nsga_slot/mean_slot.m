% % =========================================================================
% % 论文数据分析与绘图脚本 (整合版)
% % 新增：全局样式配置，确保所有图表风格统一
% % =========================================================================
% 
% clc;        % 清理命令行窗口
% clear;      % 清理工作区变量
% close all;  % 关闭所有已打开的图形窗口
% 
% %% 步骤 0: 全局样式与文本配置
% fprintf('步骤 0: 配置全局样式...\n');
% 
% % --- 在这里选择你想要的风格 ---
% % 'StyleA' -> 橙/蓝/绿/紫，不同标记，适合区分度高的场景
% % 'StyleB' -> 青/品红/橙/黄，统一圆点，适合简洁风格
% plot_style_choice = 'StyleB'; % <--- 修改这里来切换风格
% 
% % --- 全局文本配置 ---
% % 无论何种风格，都使用统一的文本标签
% global_title = {'S4'};
% global_xlabel ={ 'G_1 (s)'};
% global_ylabel = {'G_2 (J)'};
% 
% % --- 根据选择配置详细样式 ---
% if strcmpi(plot_style_choice, 'StyleA')
%     fprintf('使用风格: StyleA (多样标记)\n\n');
%     colors = [
%         0.8500, 0.3250, 0.0980; % 橙色 (Ours)
%         0, 0.4470, 0.7410;      % 蓝色 (NSGA-II)
%         0.4660, 0.6740, 0.1880; % 绿色 (A-MOPSO)
%         0.4940, 0.1840, 0.5560; % 紫色 (ARG)
%     ];
%     markers = {'^', 's', 'd', 'p'}; % 三角, 方块, 菱形, 五角星
%     knee_point_marker = 'o'; % 拐点用大圆圈
% 
% elseif strcmpi(plot_style_choice, 'StyleB')
%     fprintf('使用风格: StyleB (统一圆点)\n\n');
%     colors = [
%         [0, 0.8, 0.8];   % 青色 (Ours)
%         [1, 0, 1];       % 品红 (NSGA-II)
%         [1, 0.5, 0];     % 橙色 (A-MOPSO)
%         [1, 0.8, 0];     % 黄色 (ARG)
%     ];
%     markers = {'o', 'o', 'o', 'o'}; % 所有算法都用实心圆
%     knee_point_marker = 'p'; % 拐点用五角星
% else
%     error('未知的风格选项。请选择 "StyleA" 或 "StyleB"。');
% end
% 
% % 通用样式参数
% font_name = 'Arial';
% font_size = 12;
% line_width = 1.5;
% 
% %% 步骤 1: 配置与加载 S6 场景数据
% 
% fprintf('步骤 1: 配置并加载 S6 场景数据...\n');
% base_path = 'E:\guthub-matlab\第二篇论文数据'; 
% s6_file_name = '30-5-45-200-200.mat';
% full_s6_path = fullfile(base_path, s6_file_name);
% algonames_legend = {'Ours', 'NSGA-II', 'A-MOPSO', 'ARG'};
% algonames_data = {'MyNSGA_II', 'NSGA_II', 'MOPSO', 'Baseline'};
% 
% % 加载数据
% if ~isfile(full_s6_path)
%     error('错误: 在指定路径下找不到 S6 的数据文件: %s', full_s6_path);
% end
% fprintf('正在加载文件: %s\n', full_s6_path);
% s6_data = load(full_s6_path);
% fprintf('S6 数据加载成功！\n\n');
% 
% if ~isfield(s6_data, 'all_scenario_results') || ~isfield(s6_data, 'alg_names_for_results')
%     error('错误: 加载的.mat文件中未找到 "all_scenario_results" 或 "alg_names_for_results" 变量。');
% end
% 
% %% 步骤 2: 寻找每个算法的代表性运行 (基于最终时隙的IGD)
% fprintf('步骤 2: 寻找每个算法的代表性运行 (基于最终时隙的IGD)...\n');
% final_slot_idx = 5; 
% fprintf('分析的时隙: 第 %d 个时隙\n', final_slot_idx);
% all_scenario_results = s6_data.all_scenario_results;
% alg_names = s6_data.alg_names_for_results;
% num_algs = length(alg_names);
% representative_runs = struct();
% s_idx = 1;
% fprintf('------------------------------------------------------------\n');
% for i = 1:num_algs
%     alg_data_field = algonames_data{i}; 
%     igd_values_final_slot = all_scenario_results.(alg_data_field).IGD{s_idx}(:, final_slot_idx);
%     median_igd = median(igd_values_final_slot, 'omitnan');
%     distances_to_median = abs(igd_values_final_slot - median_igd);
%     [~, run_index] = min(distances_to_median);
%     representative_runs.(alg_data_field) = run_index;
%     fprintf('算法: %-12s | IGD中位数: %.4f | 代表性运行: 第 %2d 次 (IGD: %.4f)\n', ...
%             algonames_legend{i}, median_igd, run_index, igd_values_final_slot(run_index));
% end
% fprintf('------------------------------------------------------------\n');
% fprintf('\n代表性运行查找完毕！\n');
% fprintf('结果已保存在变量 "representative_runs" 中。\n\n');
% disp(representative_runs);
% 
% %% 步骤 3: 绘制代表性帕累托前沿对比图 (使用全局样式)
% fprintf('步骤 3: 正在绘制代表性帕累托前沿对比图 (使用全局样式)...\n');
% figure('Name', 'Representative Pareto Fronts Comparison');
% hold on; 
% grid on; 
% box on;
% if isfield(s6_data, 'igd_reference_front_obj')
%     pf_star = s6_data.igd_reference_front_obj;
%     plot(pf_star(:,1), pf_star(:,2), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 12, 'DisplayName', 'PF*');
% end
% all_fronts = s6_data.current_scenario_all_run_final_fronts;
% knee_points = zeros(num_algs, 2);
% for i = 1:num_algs
%     alg_data_field = algonames_data{i};
%     run_idx = representative_runs.(alg_data_field);
%     front = all_fronts{i, final_slot_idx, run_idx};
%     if isempty(front), continue; end
%     plot(front(:,1), front(:,2), ...
%         'Marker', markers{i}, ...
%         'MarkerEdgeColor', 'k', ...
%         'MarkerFaceColor', colors(i,:), ...
%         'LineStyle', 'none', ...
%         'MarkerSize', 8, ...
%         'DisplayName', algonames_legend{i});
%     knee_point = find_knee_point(front);
%     if ~isempty(knee_point)
%         plot(knee_point(1), knee_point(2), ...
%              knee_point_marker, ...
%              'MarkerEdgeColor', 'k', ...
%              'MarkerFaceColor', colors(i,:), ...
%              'MarkerSize', 12, 'LineWidth', line_width, ...
%              'HandleVisibility', 'off');
%         knee_points(i, :) = knee_point;
%     end
% end
% title(global_title, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', font_name);
% xlabel(global_xlabel, 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', font_name);
% ylabel(global_ylabel, 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', font_name);
% legend('show', 'Location', 'northeast', 'FontSize', font_size, 'FontName', font_name);
% set(gca, 'FontSize', font_size, 'LineWidth', 0.6, 'FontName', font_name);
% ax = gca;
% x_lim = get(ax, 'XLim');
% y_lim = get(ax, 'YLim');
% set(ax, 'XLim', [x_lim(1) - 0.05*diff(x_lim), x_lim(2) + 0.05*diff(x_lim)]);
% set(ax, 'YLim', [y_lim(1) - 0.05*diff(y_lim), y_lim(2) + 0.05*diff(y_lim)]);
% hold off;
% fprintf('前沿对比图绘制完成。\n\n');
% fprintf('--- 各算法代表性前沿的拐点方案对比 ---\n');
% fprintf('%-12s | %-12s | %-12s\n', 'Algorithm', 'G1 (s)', 'G2 (J)');
% fprintf('-----------------------------------------------\n');
% for i = 1:num_algs
%     fprintf('%-12s | %.4f       | %.4f\n', ...
%             algonames_legend{i}, knee_points(i,1), knee_points(i,2));
% end
% fprintf('-----------------------------------------------\n\n');
% 
% %% 步骤 4: 绘制性能随时间演进的折线图
% fprintf('步骤 4: 正在绘制性能随时间演进的折线图...\n');
% 
% % 提取时隙数量
% nSlots = size(all_scenario_results.MyNSGA_II.IGD{s_idx}, 2);
% 
% % 创建一个新的图形窗口
% figure('Name', 'Performance Evolution Over Time Slots');
% 
% % --- 绘制IGD演进图 ---
% 
% hold on; grid on; box on;
% 
% for i = 1:num_algs
%     alg_data_field = algonames_data{i};
%     % 计算每个时隙的平均IGD (对30次运行求均值)
%     mean_igd_over_slots = mean(all_scenario_results.(alg_data_field).IGD{s_idx}, 1, 'omitnan');
%     plot(1:nSlots, mean_igd_over_slots, ...
%          'Color', colors(i,:), ...
%          'Marker', markers{i}, ...
%          'MarkerFaceColor', colors(i,:), ... % 使标记实心
%          'LineStyle', '-', ...
%          'LineWidth', line_width, ...
%          'DisplayName', algonames_legend{i});
% end
% title('IGD Evolution', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', font_name);
% xlabel('Time Slot', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', font_name);
% ylabel('Average IGD', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', font_name);
% legend('show', 'Location', 'northeast', 'FontSize', font_size);
% set(gca, 'FontSize', font_size, 'LineWidth', 0.6, 'FontName', font_name);
% xticks(1:nSlots); % 确保x轴刻度为整数
% 
% % --- 绘制HV演进图 ---
% % 创建一个新的图形窗口
% figure('Name', 'Performance Evolution Over Time Slots');
% hold on; grid on; box on;
% 
% for i = 1:num_algs
%     alg_data_field = algonames_data{i};
%     % 计算每个时隙的平均HV
%     mean_hv_over_slots = mean(all_scenario_results.(alg_data_field).HV{s_idx}, 1, 'omitnan');
%     plot(1:nSlots, mean_hv_over_slots, ...
%          'Color', colors(i,:), ...
%          'Marker', markers{i}, ...
%          'MarkerFaceColor', colors(i,:), ...
%          'LineStyle', '-', ...
%          'LineWidth', line_width, ...
%          'DisplayName', algonames_legend{i});
% end
% title('HV Evolution', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', font_name);
% xlabel('Time Slot', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', font_name);
% ylabel('Average HV', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', font_name);
% legend('show', 'Location', 'southeast', 'FontSize', font_size);
% set(gca, 'FontSize', font_size, 'LineWidth', 0.6, 'FontName', font_name);
% xticks(1:nSlots);
% 
% fprintf('性能演进图绘制完成。\n\n');
% fprintf('================== 所有分析与绘图已完成 ==================\n');
% 
% %% 辅助函数：寻找拐点 (Knee Point)
% function knee_point = find_knee_point(front)
%     if isempty(front) || size(front, 1) < 3
%         knee_point = [];
%         return;
%     end
%     min_vals = min(front, [], 1);
%     max_vals = max(front, [], 1);
%     range_vals = max_vals - min_vals;
%     range_vals(range_vals == 0) = 1;
%     normalized_front = (front - min_vals) ./ range_vals;
%     [~, idx1] = min(normalized_front(:,1));
%     [~, idx2] = min(normalized_front(:,2));
%     extreme_point1 = normalized_front(idx1, :);
%     extreme_point2 = normalized_front(idx2, :);
%     v1 = extreme_point2 - extreme_point1;
%     distances = zeros(size(normalized_front, 1), 1);
%     for i = 1:size(normalized_front, 1)
%         v2 = normalized_front(i, :) - extreme_point1;
%         distances(i) = abs(v1(1)*v2(2) - v1(2)*v2(1)) / norm(v1);
%     end
%     [~, knee_idx] = max(distances);
%     knee_point = front(knee_idx, :);
% end

% =========================================================================
% 论文数据分析与绘图脚本 (LaTeX 渲染整合版)
% 新增：全局 LaTeX 样式配置，确保所有标签符合学术出版标准
% =========================================================================
clc; clear; close all;

%% 步骤 0: 全局样式与文本配置
fprintf('步骤 0: 配置全局样式...\n');

% --- 风格选择 ---
plot_style_choice = 'StyleB'; 

% --- 全局文本配置 (LaTeX 语法) ---
% 在 LaTeX 模式下，数学变量建议使用 $ $ 包裹，加粗使用 \textbf{}
global_title = 'Representative Pareto Fronts Comparison (S4)';
global_xlabel = '$G_1$ (s)';
global_ylabel = '$G_2$ (J)';

% --- 根据选择配置样式 ---
if strcmpi(plot_style_choice, 'StyleA')
    colors = [0.8500, 0.3250, 0.0980; 0, 0.4470, 0.7410; 0.4660, 0.6740, 0.1880; 0.4940, 0.1840, 0.5560];
    markers = {'^', 's', 'd', 'p'};
    knee_point_marker = 'o';
elseif strcmpi(plot_style_choice, 'StyleB')
    colors = [[0, 0.8, 0.8]; [1, 0, 1]; [1, 0.5, 0]; [1, 0.8, 0]];
    markers = {'o', 'o', 'o', 'o'};
    knee_point_marker = 'p';
end

% 通用样式参数
font_name = 'Times New Roman'; % 配合 LaTeX 渲染的常用衬线字体
font_size = 12;
line_width = 1.5;

%% 步骤 1: 配置与加载数据
fprintf('步骤 1: 加载数据文件...\n');
base_path = 'E:\guthub-matlab\第二篇论文数据'; 
s6_file_name = '30-5-45-200-200.mat';
full_s6_path = fullfile(base_path, s6_file_name);
algonames_legend = {'Ours', 'NSGA-II', 'A-MOPSO', 'ARG'};
algonames_data = {'MyNSGA_II', 'NSGA_II', 'MOPSO', 'Baseline'};

if ~isfile(full_s6_path)
    error('错误: 找不到数据文件: %s', full_s6_path);
end
s6_data = load(full_s6_path);

%% 步骤 2: 寻找代表性运行 (基于 IGD)
fprintf('步骤 2: 分析 IGD 中位数运行...\n');
final_slot_idx = 5; 
all_scenario_results = s6_data.all_scenario_results;
num_algs = length(algonames_data);
representative_runs = struct();
s_idx = 1;

for i = 1:num_algs
    alg_field = algonames_data{i};
    igd_vals = all_scenario_results.(alg_field).IGD{s_idx}(:, final_slot_idx);
    [~, run_idx] = min(abs(igd_vals - median(igd_vals, 'omitnan')));
    representative_runs.(alg_field) = run_idx;
end

%% 步骤 3: 绘制代表性帕累托前沿对比图 (LaTeX 渲染)
fprintf('步骤 3: 绘制 Pareto Front 对比图...\n');
figure('Color', 'w', 'Name', 'Pareto Fronts Comparison');
hold on; grid on; box on;

% 绘制 PF* (参考前沿)
if isfield(s6_data, 'igd_reference_front_obj')
    pf_star = s6_data.igd_reference_front_obj;
    plot(pf_star(:,1), pf_star(:,2), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 10, 'DisplayName', '$PF^*$');
end

all_fronts = s6_data.current_scenario_all_run_final_fronts;
knee_points = zeros(num_algs, 2);

for i = 1:num_algs
    alg_field = algonames_data{i};
    run_idx = representative_runs.(alg_field);
    front = all_fronts{i, final_slot_idx, run_idx};
    if isempty(front), continue; end
    
    % 绘制前沿点
    plot(front(:,1), front(:,2), 'Marker', markers{i}, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', colors(i,:), 'LineStyle', 'none', 'MarkerSize', 8, ...
        'DisplayName', algonames_legend{i});
    
    % 绘制拐点 (Knee Point)
    kp = find_knee_point(front);
    if ~isempty(kp)
        plot(kp(1), kp(2), knee_point_marker, 'MarkerEdgeColor', 'k', ...
             'MarkerFaceColor', colors(i,:), 'MarkerSize', 12, 'LineWidth', line_width, ...
             'HandleVisibility', 'off');
        knee_points(i, :) = kp;
    end
end

% 应用 LaTeX 解释器
title(global_title, 'FontSize', 14, 'Interpreter', 'latex', 'FontName', font_name);
xlabel(global_xlabel, 'FontSize', font_size, 'Interpreter', 'latex', 'FontName', font_name);
ylabel(global_ylabel, 'FontSize', font_size, 'Interpreter', 'latex', 'FontName', font_name);
legend('show', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', font_size, 'FontName', font_name);
hold off;

%% 步骤 4: 绘制性能随时间演进的折线图 (LaTeX 渲染)
nSlots = size(all_scenario_results.MyNSGA_II.IGD{s_idx}, 2);

% --- IGD 演进图 ---
figure('Color', 'w', 'Name', 'IGD Evolution');
hold on; grid on; box on;
for i = 1:num_algs
    mean_igd = mean(all_scenario_results.(algonames_data{i}).IGD{s_idx}, 1, 'omitnan');
    plot(1:nSlots, mean_igd, 'Color', colors(i,:), 'Marker', markers{i}, ...
         'MarkerFaceColor', colors(i,:), 'LineWidth', line_width, 'DisplayName', algonames_legend{i});
end
title('IGD Evolution over Time Slots', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time Slot $t$', 'Interpreter', 'latex');
ylabel('Average IGD', 'Interpreter', 'latex');
legend('show', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'XTick', 1:nSlots);

% --- HV 演进图 ---
figure('Color', 'w', 'Name', 'HV Evolution');
hold on; grid on; box on;
for i = 1:num_algs
    mean_hv = mean(all_scenario_results.(algonames_data{i}).HV{s_idx}, 1, 'omitnan');
    plot(1:nSlots, mean_hv, 'Color', colors(i,:), 'Marker', markers{i}, ...
         'MarkerFaceColor', colors(i,:), 'LineWidth', line_width, 'DisplayName', algonames_legend{i});
end
title('HV Evolution over Time Slots', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time Slot $t$', 'Interpreter', 'latex');
ylabel('Average HV', 'Interpreter', 'latex');
legend('show', 'Location', 'southeast', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'XTick', 1:nSlots);

fprintf('================== 所有分析与绘图已完成 (LaTeX 渲染) ==================\n');

%% 辅助函数：寻找拐点 (Knee Point)
function knee_point = find_knee_point(front)
    if isempty(front) || size(front, 1) < 3
        knee_point = []; return;
    end
    min_vals = min(front, [], 1);
    max_vals = max(front, [], 1);
    norm_front = (front - min_vals) ./ (max_vals - min_vals);
    [~, idx1] = min(norm_front(:,1));
    [~, idx2] = min(norm_front(:,2));
    p1 = norm_front(idx1, :); p2 = norm_front(idx2, :);
    v = p2 - p1;
    dist = abs(v(1)*(norm_front(:,2)-p1(2)) - v(2)*(norm_front(:,1)-p1(1))) / norm(v);
    [~, k_idx] = max(dist);
    knee_point = front(k_idx, :);
end