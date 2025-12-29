% =========================================================================
% 算法性能分布分析脚本 (单场景独立版 - 扩展 Spacing 输出)
% 功能：为 6 个场景分别生成 IGD, HV 和 Spacing 的箱线图（共 18 张图）
% =========================================================================
clc; clear; close all;

%% --- 步骤 0: 配置信息 ---
base_folder = "E:\guthub-matlab\第二篇论文数据\";
alg_legend = {'Ours', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III'};
alg_fields = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III'};

% --- 核心修改：增加 Spacing 指标 ---
metrics = {'IGD', 'HV', 'Spacing'}; 
target_slot = 10; 

% 学术配色 (保持与趋势图一致)
colors = [
    1.0, 0.5, 0.0;  % 橙色 (Ours)
    0.0, 0.4, 0.8;  % 蓝色 (NSGA-II)
    0.0, 0.8, 0.1;  % 绿色 (MOEA/D)
    0.0, 1.0, 1.0;  % 青色 (INSGA-II)
    1.0, 0.0, 1.0;  % 紫色 (NSGA-III)
];

%% --- 步骤 1: 场景大循环 (S1 - S6) ---
for s_idx = 1:6
    file_path = fullfile(base_folder, sprintf('S%d_V5.mat', s_idx));
    if ~isfile(file_path)
        fprintf('跳过：未找到文件 S%d\n', s_idx);
        continue; 
    end
    data_struct = load(file_path);
    
    % --- 步骤 2: 指标循环 (IGD, HV, Spacing) ---
    for m = 1:length(metrics)
        curr_metric = metrics{m};
        
        % 画布大小适配单场景
        fig = figure('Color', 'w', 'Position', [200, 200, 550, 480], 'Visible', 'off'); 
        hold on; grid on; box on;
        
        boxplot_data = [];
        group_idx = [];
        
        % 提取数据
        for a = 1:length(alg_fields)
            alg_name = alg_fields{a};
            try
                if isfield(data_struct.all_scenario_results, alg_name)
                    metric_cell = data_struct.all_scenario_results.(alg_name).(curr_metric);
                    
                    % 动态兼容索引逻辑 (适配单个场景存为一个 cell 或多个的情况)
                    if length(metric_cell) >= s_idx && ~isempty(metric_cell{s_idx})
                        run_data = metric_cell{s_idx}(:, target_slot);
                    else
                        run_data = metric_cell{1}(:, target_slot);
                    end
                    
                    boxplot_data = [boxplot_data; run_data];
                    group_idx = [group_idx; repmat(a, length(run_data), 1)];
                end
            catch
                fprintf('数据缺失: S%d-%s-%s\n', s_idx, alg_name, curr_metric);
            end
        end
        
        if isempty(boxplot_data), close(fig); continue; end
        
        % --- 步骤 3: 绘制带缺口的学术箱线图 ---
        % 'Notch', 'on' 展示中位数的置信区间
        boxplot(boxplot_data, group_idx, 'Notch', 'on', 'Symbol', 'r+', 'Widths', 0.6, 'Labels', alg_legend);
        
        % 美化：颜色填充
        h = findobj(gca, 'Tag', 'Box');
        num_boxes = length(h);
        for j = 1:num_boxes
            % 这里的索引需要反转，因为 findobj 找到的顺序是倒序
            alg_idx = num_boxes - j + 1;
            patch(get(h(j), 'XData'), get(h(j), 'YData'), colors(alg_idx, :), ...
                'FaceAlpha', 0.45, 'EdgeColor', colors(alg_idx, :), 'LineWidth', 1.2);
        end
        
        % 坐标轴美化
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12, 'LineWidth', 1.1);
        xtickangle(35); 
        
        % 动态 Y 轴标签
        ylabel(sprintf('$\\mathrm{%s}$ at $t=10$', curr_metric), 'Interpreter', 'latex', 'FontSize', 14);
        title(sprintf('Statistical Analysis: %s (S%d)', curr_metric, s_idx), 'Interpreter', 'latex', 'FontSize', 14);
        
        % --- 步骤 4: 高清导出 ---
        save_name = sprintf('Boxplot_%s_S%d.png', curr_metric, s_idx);
        exportgraphics(fig, save_name, 'Resolution', 300);
        
        close(fig); % 及时关闭窗口释放内存
    end
    fprintf('场景 S%d 处理完成。\n', s_idx);
end
fprintf('>>> 全部 18 张独立箱线图已生成并保存至当前文件夹。\n');