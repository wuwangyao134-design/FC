clc; clear; close all;

%% --- 步骤 0: 全局环境与 LaTeX 引擎设置 ---
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

% 基础路径配置
base_folder = "E:\guthub-matlab\ExperimentResults_Output\"; 
file_name = "S6.mat"; 
full_path = fullfile(base_folder, file_name);

save_folder = "E:\guthub-matlab\第二篇论文数据\";
if ~exist(save_folder, 'dir'), mkdir(save_folder); end

% 算法配置
alg_legend = {'Ours', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III', 'MODDPG'};
alg_fields = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III', 'DRL_Baseline'};
metrics = {'IGD', 'HV', 'Spacing'}; 
sub_labels = {'(a)', '(b)', '(c)'};
s_idx = 6; % 锁定 S6 场景
target_slot = 10; % 锁定最后一个时间槽

% =========================================================================
% 🎨 学术配色方案
% =========================================================================
colors = [
    0.70, 0.10, 0.55; % 紫红色 (Ours)
    0.00, 0.40, 0.80; % 蓝色 (NSGA-II)
    0.00, 0.80, 0.10; % 绿色 (MOEA/D)
    0.00, 0.55, 0.55; % 孔雀蓝/深青色 (INSGA-II)
    1.00, 0.00, 1.00; % 紫色 (NSGA-III)
    0.60, 0.20, 0.20; % 棕红色 (MODDPG)
];
fontSize = 11; 

% 🌟 核心对齐参数：强制锁定所有图的内部主坐标轴位置与比例！
unified_main_axes_pos = [0.12, 0.15, 0.83, 0.75]; 

%% --- 步骤 1: 数据提取 ---
if ~isfile(full_path)
    error('⚠️ 未找到文件：%s', full_path);
end
data_struct = load(full_path);

%% --- 步骤 2: 指标循环绘图 ---
for m = 1:length(metrics)
    curr_metric = metrics{m};
    
    % 🌟 画布尺寸：固定单栏观感
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [2 2 4.5 3.5], 'Visible', 'off'); 
    
    % =====================================================================
    % 🌟 1. 绘制主图 (强制统一框线尺寸)
    % =====================================================================
    ax_main = axes('Position', unified_main_axes_pos); 
    hold(ax_main, 'on'); grid(ax_main, 'on'); box(ax_main, 'on');
    set(ax_main, 'GridLineStyle', '-', 'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.5);
    
    boxplot_data = [];
    group_idx = [];
    
    for a = 1:length(alg_fields)
        alg_name = alg_fields{a};
        if isfield(data_struct.all_scenario_results, alg_name)
            metric_cell = data_struct.all_scenario_results.(alg_name).(curr_metric);
            if length(metric_cell) >= s_idx && ~isempty(metric_cell{s_idx})
                run_data = metric_cell{s_idx}(:, target_slot);
            else
                run_data = metric_cell{1}(:, target_slot);
            end
            run_data(isnan(run_data)) = []; 
            
            % 🚨 兜底防爆设置 (配合画中画的 -10，这里设为 1e-11)
            if strcmp(curr_metric, 'IGD') || strcmp(curr_metric, 'Spacing')
                run_data(run_data <= 0) = 1e-11; 
            end
            
            % 🚨【重要修正】：Spacing 专属底层 log10 转换恢复！
            if strcmp(curr_metric, 'Spacing')
                run_data = log10(run_data);
            end
            
            boxplot_data = [boxplot_data; run_data];
            group_idx = [group_idx; repmat(a, length(run_data), 1)];
        end
    end
    
    if isempty(boxplot_data), close(fig); continue; end
    
    % 主图 Boxplot
    boxplot(ax_main, boxplot_data, group_idx, 'Notch', 'on', 'Symbol', 'r+', 'Widths', 0.5, 'Labels', alg_legend);
    
    % 主图颜色填充
    h_main = findobj(ax_main, 'Tag', 'Box');
    num_boxes = length(h_main);
    for j = 1:num_boxes
        alg_idx = num_boxes - j + 1;
        patch(get(h_main(j), 'XData'), get(h_main(j), 'YData'), colors(alg_idx, :), ...
            'FaceAlpha', 0.45, 'EdgeColor', colors(alg_idx, :), 'LineWidth', 1.2, 'Parent', ax_main);
    end
    
    % 主图标题美化
    title_str = sprintf('%s %s Boxplot (S%d)', sub_labels{m}, curr_metric, s_idx);
    title(ax_main, title_str, 'Interpreter', 'latex', 'FontSize', fontSize + 1);
    set(ax_main, 'TickLabelInterpreter', 'latex', 'FontSize', fontSize, 'LineWidth', 1.1, 'TickDir', 'in');
    xtickangle(ax_main, 25); 
    
   % 🚨 动态更新 Y 轴标签与 Spacing 的极限锁定
    if strcmp(curr_metric, 'Spacing')
        ylabel(ax_main, sprintf('$\\log_{10}(\\mathrm{%s})$ at $t=%d$', curr_metric, target_slot), 'Interpreter', 'latex', 'FontSize', fontSize);
        
        % 🚨【强制卡点】：Spacing 锁定在 [-6, 1.5]
        ylim(ax_main, [-6, 1.5]);
        
        % 🚨【完美步长】：以 1.5 为步长的优雅刻度 [-6, -4.5, -3, -1.5, 0, 1.5]
        yticks(ax_main, -6:1.5:1.5);
    else
        ylabel(ax_main, sprintf('$\\mathrm{%s}$ at $t=%d$', curr_metric, target_slot), 'Interpreter', 'latex', 'FontSize', fontSize);
    end

    % =====================================================================
    % 🌟 2. IGD 专属画中画【Ours 和 MOEA/D 单挑，范围 -10 到 1】
    % =====================================================================
    if strcmp(curr_metric, 'IGD')
        
        % 📐 定位左上方
        inset_pos = [0.18, 0.49, 0.35, 0.32]; 
        ax_inset = axes('Position', inset_pos);
        hold(ax_inset, 'on'); grid(ax_inset, 'on'); box(ax_inset, 'on');
        
        % 提取 Ours(1) 和 MOEA/D(3)
        selected_algs = [1, 3]; 
        inset_mask = (group_idx == 1) | (group_idx == 3);
        
        % 底层数据取 log10 数学转换
        inset_boxplot_data = log10(boxplot_data(inset_mask));
        inset_group_idx = group_idx(inset_mask);
        
        % 绘制画中画的 Boxplot
        boxplot(ax_inset, inset_boxplot_data, inset_group_idx, 'Notch', 'on', 'Symbol', 'r+', ...
            'Widths', 0.5, 'Labels', alg_legend(selected_algs));
        
        % 画中画颜色填充
        h_inset = findobj(ax_inset, 'Tag', 'Box');
        num_inset_boxes = length(h_inset);
        for j = 1:num_inset_boxes
            alg_idx = selected_algs(num_inset_boxes - j + 1); 
            patch(get(h_inset(j), 'XData'), get(h_inset(j), 'YData'), colors(alg_idx, :), ...
                'FaceAlpha', 0.6, 'EdgeColor', colors(alg_idx, :), 'LineWidth', 1.0, 'Parent', ax_inset);
        end
        
        % 🚨【极限卡点】：Y 轴死锁在 -10 到 1
        ylim(ax_inset, [-12, 1]); 
        
        % 🚨【清爽刻度】：-10, -8, -6, -4, -2, 0, 1
        yticks(ax_inset, [-12,-10, -8, -6, -4, -2, 0, 1]);
        
        % 画中画网格线与格式美化
        set(ax_inset, 'GridLineStyle', '-', 'GridColor', [0.7 0.7 0.7]);
        set(ax_inset, 'TickLabelInterpreter', 'latex', 'FontSize', 8, 'LineWidth', 0.8);
        title(ax_inset, '$\log_{10}(\mathrm{IGD}) \in [-12, 1]$', 'Interpreter', 'latex', 'FontSize', 9);
    end
    % =====================================================================
    
    % --- 步骤 4: 高保真纯矢量 PDF 导出 ---
    save_name_pdf = sprintf('Boxplot_%s_S%d.pdf', curr_metric, s_idx);
    full_save_path = fullfile(save_folder, save_name_pdf);
    
    exportgraphics(fig, full_save_path, 'ContentType', 'vector');
    fprintf('✅ 成功导出 PDF: %s\n', save_name_pdf);
    close(fig);
end
fprintf('🎉 终极完全体优化完毕！(Spacing 恢复全局 log10，所有尺寸与画中画完美适配)！\n');