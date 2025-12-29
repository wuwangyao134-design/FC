% =========================================================================
% HMD-NSGA-II 消融实验绘图 - 包含 Spacing 指标专业版
% =========================================================================

%% 1. 全局环境与 LaTeX 设置
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

% --- 修改点：增加 Spacing 指标配置 ---
metrics = {'IGD', 'HV', 'Spacing', 'NumFeasibleSolutions'};
y_labels = {'Average $\mathrm{IGD}$ Value', ...
            'Average $\mathrm{HV}$ Value', ...
            'Average $\mathrm{Spacing}$ Value', ...
            'Feasibility Rate ($\%$)'};
file_names = {'Ablation_Trend_IGD_Updated.png', ...
              'Ablation_Trend_HV_Updated.png', ...
              'Ablation_Trend_Spacing_Updated.png', ...
              'Ablation_Trend_Feas_Updated.png'};

% 算法配置与顶刊配色
alg_list = {'HMD_no_Hybrid', 'HMD_no_ISMM', 'HMD_no_AMS', 'HMD_Full', 'Standard_NSGAII'};
labels = {'Proposed', 'Variant A (w/o ISMM)', 'Variant B (w/o AMS)', ...
          'Variant C (w/o HO)', 'Baseline'};
colors = [
    0.90, 0.12, 0.06; % Proposed: 活力红
    0.12, 0.47, 0.71; % Variant A: 宝石蓝
    0.17, 0.63, 0.17; % Variant B: 森林绿
    1.00, 0.50, 0.00; % Variant C: 琥珀橙
    0.42, 0.24, 0.60; % Baseline: 灰调紫
];
markers = {'o', 's', '^', 'd', 'x'};

if ~exist('nSlots', 'var'), nSlots = 10; end
slots = 1:nSlots;

%% 2. 批处理绘图
for m = 1:length(metrics)
    curr_metric = metrics{m};
    
    % 创建画布
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [2, 2, 6, 4.5]);
    hold on; grid on;
    plot_handles = [];
    
    for i = 1:length(alg_list)
        name = alg_list{i};
        
        % 直接从修复后的 all_results 中提取数据
        data = all_results.(name).(curr_metric);
        
        % 统计计算
        avg_val = mean(data, 1);
        std_val = std(data, 0, 1);
        
        % 绘制标准差阴影
        fill([slots, fliplr(slots)], [avg_val-std_val, fliplr(avg_val+std_val)], ...
            colors(i,:), 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % 绘制主曲线
        h = plot(slots, avg_val, ['-', markers{i}], 'Color', colors(i,:), ...
            'LineWidth', 1.8, 'MarkerSize', 5, 'MarkerFaceColor', colors(i,:));
        plot_handles = [plot_handles, h];
    end
    
    % --- 坐标轴样式优化 ---
    set(gca, 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'Box', 'on', 'LineWidth', 1.2);
    xlim([0.9, nSlots + 0.5]);
    set(gca, 'XTick', [1, 2, 5, 10]);
    xticklabels({'$10^0$', '2', '5', '$10^1$'}); 
    
    ylabel(y_labels{m}, 'Interpreter', 'latex', 'FontSize', 12, 'Rotation', 90, 'VerticalAlignment', 'bottom');
    xlabel('Time Slot Index $t$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['Ablation Analysis: $\mathrm{', strrep(curr_metric,'_','\_'), '}$ Performance'], 'Interpreter', 'latex', 'FontSize', 12);
    
    % 图例仅在 HV 或 Spacing 图中展示（你可以根据需要修改）
    if strcmp(curr_metric, 'HV')
        legend(plot_handles, labels, 'Interpreter', 'latex', 'FontSize', 8, ...
               'Location', 'best', 'NumColumns', 1);
    end
    
    % 导出
    exportgraphics(fig, file_names{m}, 'Resolution', 300);
    fprintf('指标 %s 趋势图已保存至: %s\n', curr_metric, file_names{m});
end

%% 3. 输出汇总统计表 (包含 Spacing)
% --- 修改点：增加 Spacing 列，加宽表格 ---
fprintf('\n%s\n', repmat('=', 1, 115));
fprintf('%-20s | %-18s | %-18s | %-18s | %-12s\n', ...
    'Algorithm Variant', 'IGD (Mean±Std)', 'HV (Mean±Std)', 'Spacing (Mean±Std)', 'Feasible (%)');
fprintf('%s\n', repmat('-', 1, 115));

for i = 1:length(alg_list)
    n = alg_list{i};
    fIGD = all_results.(n).IGD(:, end); 
    fHV = all_results.(n).HV(:, end); 
    fSpacing = all_results.(n).Spacing(:, end); % 提取 Spacing 数据
    fFeas = all_results.(n).NumFeasibleSolutions(:, end);
    
    % 打印包含 Spacing 的完整数据
    fprintf('%-20s | %.4f ± %.4f | %.4f ± %.4f | %.4f ± %.4f | %.2f%%\n', ...
        labels{i}, mean(fIGD), std(fIGD), mean(fHV), std(fHV), mean(fSpacing), std(fSpacing), mean(fFeas));
end
fprintf('%s\n', repmat('=', 1, 115));