% clc; clear; close all;
% 
% %% --- 步骤 0: 配置信息 ---
% base_folder = "E:\guthub-matlab\第二篇论文数据\"; 
% file_name = "S5_V5.mat"; 
% alg_legend = {'Ours', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III'};
% alg_fields = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III'};
% metrics = {'IGD', 'HV', 'Spacing'}; 
% s_idx = 1; 
% 
% % 严格匹配示例图的配色
% colors = [
%     1.0, 0.5, 0.0;  % 橙色 (Ours)
%     0.0, 0.44, 0.74; % 深蓝色 (NSGA-II)
%     0.15, 0.65, 0.15; % 绿色 (MOEA/D)
%     0.0, 0.9, 0.9;   % 青色 (INSGA-II)
%     0.8, 0.0, 0.8;   % 紫色 (NSGA-III)
% ];
% % 匹配示例图的 Marker 形状
% markers = {'o', 's', '^', 'd', 'v'};
% 
% %% --- 步骤 1: 数据提取 ---
% full_path = fullfile(base_folder, file_name);
% if ~isfile(full_path)
%     error('未找到文件：%s', full_path);
% end
% data_struct = load(full_path);
% nSlots = size(data_struct.all_scenario_results.MyNSGA_II.IGD{s_idx}, 2);
% time_slots = 1:nSlots;
% 
% %% --- 步骤 2: 指标循环绘图 ---
% for m = 1:length(metrics)
%     curr_metric = metrics{m};
% 
%     % 调整画布比例，匹配示例图的方正感
%     fig = figure('Color', 'w', 'Position', [200, 200, 650, 500]);
%     hold on;
% 
%     % 强制开启网格并设置样式
%     grid on;
%     set(gca, 'GridLineStyle', '-', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.5);
% 
%     h = zeros(length(alg_fields), 1);
% 
%     for a = 1:length(alg_fields)
%         alg_name = alg_fields{a};
%         metric_cell = data_struct.all_scenario_results.(alg_name).(curr_metric);
% 
%         if length(metric_cell) >= s_idx && ~isempty(metric_cell{s_idx})
%             raw_data = metric_cell{s_idx};
%         else
%             raw_data = metric_cell{1};
%         end
% 
%         mean_trajectory = mean(raw_data, 1, 'omitnan');
% 
%         % 样式调整：加粗线条，Marker 大小调整为 8，MarkerEdgeColor 设为深色
%         h(a) = plot(time_slots, mean_trajectory, ...
%             'Color', colors(a, :), ...
%             'Marker', markers{a}, ...
%             'LineWidth', 2.5, ... 
%             'MarkerSize', 8, ...
%             'MarkerFaceColor', colors(a, :), ...
%             'MarkerEdgeColor', colors(a, :)*0.7, ... % 边缘稍深，增加立体感
%             'DisplayName', alg_legend{a});
%     end
% 
%     % 坐标轴标签美化
%     ylabel(sprintf('$\\mathrm{%s}$', curr_metric), 'Interpreter', 'latex', 'FontSize', 16);
%     xlabel('Time Slot Index $t$', 'Interpreter', 'latex', 'FontSize', 14);
% 
%     % 只有 HV 指标在左上角显示图注 (根据您的指示)
%     if strcmp(curr_metric, 'HV')
%         legend(h, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 11, 'Box', 'on');
%     else
%         legend('off');
%     end
% 
%     % 坐标轴边框与刻度设置 (模仿示例图的加粗效果)
%     set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 13, ...
%              'LineWidth', 1.5, 'Box', 'on', 'TickDir', 'in');
% 
%     % 动态调整 Y 轴范围，留出图注空间
%     y_vals = ylim;
%     ylim([y_vals(1)*0.9, y_vals(2)*1.2]); 
%     xlim([1, nSlots]);
%     xticks(1:nSlots);
% 
%     % 导出
%     save_name = sprintf('Trend_%s_S%d.png', curr_metric, s_idx);
%     exportgraphics(fig, save_name, 'Resolution', 600);
%     fprintf('成功导出趋势图: %s\n', save_name);
% end

clc; clear; close all;

%% --- 步骤 0: 配置信息 ---
% 数据读取路径
base_folder = "E:\guthub-matlab\ExperimentResults_Output\"; 
file_name = "S6.mat"; 
full_path = fullfile(base_folder, file_name);

% 图片保存路径
save_folder = "E:\guthub-matlab\第二篇论文数据\";
if ~exist(save_folder, 'dir'), mkdir(save_folder); end

% 算法配置 (加入 MODDPG)
alg_legend = {'Ours', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III', 'MODDPG'};
alg_fields = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III', 'DRL_Baseline'};
metrics = {'IGD', 'HV', 'Spacing'}; 

% 🌟 为三张图分配 (a)~(c) 标题前缀
prefixes = {'(a)', '(b)', '(c)'}; 
s_idx = 6; % 锁定 S6 场景

% 严格匹配示例图高颜值配色 + 增加 MODDPG 的沉稳棕红
colors = [
    1.00, 0.50, 0.00;  % 橙色 (Ours)
    0.00, 0.44, 0.74;  % 深蓝色 (NSGA-II)
    0.15, 0.65, 0.15;  % 绿色 (MOEA/D)
    0.00, 0.90, 0.90;  % 青色 (INSGA-II)
    0.80, 0.00, 0.80;  % 紫色 (NSGA-III)
    0.60, 0.20, 0.20;  % 棕红色 (MODDPG)
];

% 匹配示例图的 Marker 形状，并为 MODDPG 增加五角星 'p'
markers = {'o', 's', '^', 'd', 'v', 'p'};

%% --- 步骤 1: 数据提取 ---
if ~isfile(full_path)
    error('⚠️ 未找到文件：%s', full_path);
end
data_struct = load(full_path);

nSlots = size(data_struct.all_scenario_results.MyNSGA_II.IGD{end}, 2);
time_slots = 1:nSlots;

%% --- 步骤 2: 指标循环绘图 ---
for m = 1:length(metrics)
    curr_metric = metrics{m};

    % 🌟 统一尺寸缩放：宽度 450, 高度 360 (5:4 黄金比例，完美对齐你前沿图的单栏观感)
    fig = figure('Color', 'w', 'Position', [200, 200, 450, 360], 'Visible', 'off'); 
    hold on;

    % 强制开启网格并设置样式
    grid on;
    set(gca, 'GridLineStyle', '-', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.5);

    h = zeros(length(alg_fields), 1);

    for a = 1:length(alg_fields)
        alg_name = alg_fields{a};

        if ~isfield(data_struct.all_scenario_results, alg_name)
            continue;
        end

        metric_cell = data_struct.all_scenario_results.(alg_name).(curr_metric);

        if length(metric_cell) >= s_idx && ~isempty(metric_cell{s_idx})
            raw_data = metric_cell{s_idx};
        else
            raw_data = metric_cell{1};
        end

        mean_trajectory = mean(raw_data, 1, 'omitnan');

        % 样式调整
        h(a) = plot(time_slots, mean_trajectory, ...
            'Color', colors(a, :), ...
            'Marker', markers{a}, ...
            'LineWidth', 2, ... 
            'MarkerSize', 8, ...
            'MarkerFaceColor', colors(a, :), ...
            'MarkerEdgeColor', colors(a, :)*0.7, ... 
            'DisplayName', alg_legend{a});
    end

    h = h(h ~= 0); % 清理未画出的句柄

    % 🌟 新增：加入带 (a)~(c) 的子图 Title，并居中
    title_str = sprintf('%s %s Trend (S%d)', prefixes{m}, curr_metric, s_idx);
    title(title_str, 'Interpreter', 'latex', 'FontSize', 12);

    % 坐标轴标签美化
    ylabel(sprintf('$\\mathrm{%s}$', curr_metric), 'Interpreter', 'latex', 'FontSize', 11);
    xlabel('Time Slot Index $t$', 'Interpreter', 'latex', 'FontSize', 11);

    % 🌟 核心修改：仅在第一张图 (IGD) 保留图注！
    if strcmp(curr_metric, 'Spacing')
        % 放在右上角 (northeast)，双列排版防止遮挡曲线
        legend(h, 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 7, 'Box', 'on', 'NumColumns', 2);
    else
        legend('off');
    end

    % 坐标轴边框与刻度设置 (加粗效果)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11, ...
             'LineWidth', 1.1, 'Box', 'on', 'TickDir', 'in');

    % 动态调整 Y 轴范围
    y_vals = ylim;
    if strcmp(curr_metric, 'HV')
        ylim([y_vals(1)*0.9, y_vals(2)*1.1]); 
    else
        ylim([0, y_vals(2)*1.2]); % IGD 和 Spacing 从 0 起步
    end
    xlim([1, nSlots]);
    xticks(1:nSlots);

    % 导出 PNG 与 PDF
    save_name_png = sprintf('Trend_%s_S%d.png', curr_metric, s_idx);
    save_name_pdf = sprintf('Trend_%s_S%d.pdf', curr_metric, s_idx);

    exportgraphics(fig, fullfile(save_folder, save_name_png), 'Resolution', 600);
    exportgraphics(fig, fullfile(save_folder, save_name_pdf), 'ContentType', 'vector');

    fprintf('✅ 成功导出趋势图: %s\n', save_name_png);
    close(fig);
end
fprintf('🎉 所有 S6 时间序列指标图导出完毕！\n');

