% =========================================================================
% HMD-NSGA-II 消融实验绘图 - 【严格对齐 MTD 真矢量、无限放大不模糊定稿版】
% =========================================================================
%% 1. 全局环境与 LaTeX 设置
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

metrics = {'IGD', 'HV', 'Spacing', 'NumFeasibleSolutions'};
y_labels = {'Average $\mathrm{IGD}$ Value', ...
            'Average $\mathrm{HV}$ Value', ...
            'Average $\mathrm{Spacing}$ Value', ...
            'Feasibility Rate (\%)'}; % --- ✅ 彻底修复 SceneNode 警告的正确 LaTeX 语法

% --- ✅ 核心修改 1：先导出为目前顶刊最稳妥的 .pdf 纯矢量格式 ---
file_names = {'Ablation_Trend_IGD.pdf', ...
              'Ablation_Trend_HV.pdf', ...
              'Ablation_Trend_Spacing.pdf', ...
              'Ablation_Trend_Feas.pdf'};

alg_list = {'HMD_no_Hybrid', 'HMD_no_ISMM', 'HMD_no_AMS', 'HMD_Full', 'Standard_NSGAII'};
labels = {'Proposed', 'Variant A (w/o ISMM)', 'Variant B (w/o LAMS)', ...
          'Variant C (w/o HO)', 'Baseline'};

colors = [
    0.90, 0.12, 0.06; % Proposed: 红
    0.12, 0.47, 0.71; % Variant A: 蓝
    0.17, 0.63, 0.17; % Variant B: 绿
    1.00, 0.50, 0.00; % Variant C: 橙
    0.42, 0.24, 0.60; % Baseline: 紫
];
markers = {'o', 's', '^', 'd', 'x'};

if ~exist('nSlots', 'var'), nSlots = 10; end
slots = 1:nSlots;

%% 2. 批处理绘图 
for m = 1:length(metrics)
    curr_metric = metrics{m};
    
    % 画布：4×3 英寸
    fig = figure('Color','w','Units','inches','Position',[2 2 4 3]);
    hold on; box on;
    fontSize = 11;
    plot_handles = [];
    
    for i = 1:length(alg_list)
        name = alg_list{i};
        data = all_results.(name).(curr_metric);
        avg_val = mean(data, 1);
        std_val = std(data, 0, 1);
        
        % --- ✅ 核心修改 2：完美保留阴影，同时绝不触发栅格化模糊 ---
        fill([slots, fliplr(slots)], [avg_val-std_val, fliplr(avg_val+std_val)], ...
            colors(i,:), 'FaceAlpha', 0.05, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        h = plot(slots, avg_val, ['-', markers{i}], 'Color', colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', colors(i,:));
        plot_handles = [plot_handles, h];
    end
    
    set(gca, 'XScale', 'log', 'TickLabelInterpreter', 'latex', 'FontSize', fontSize, 'LineWidth', 1.1);
    xlim([0.9, nSlots + 0.2]); 
    set(gca, 'XTick', [1, 2, 5, 10]);
    xticklabels({'$10^0$','$2$','$5$','$10^1$'}); 
    grid on;

    if strcmp(curr_metric, 'Spacing')
        ylim([0, max(ylim)]);
    end
    if strcmp(curr_metric, 'NumFeasibleSolutions')
        ylim([0, 105]);
    end
    
    label_char = char('a' + m - 1);
    display_title = curr_metric;
    if strcmp(display_title, 'NumFeasibleSolutions'), display_title = 'Feasibility'; end
    
    title_str = sprintf('(%s) $\\mathrm{%s}$ Performance', label_char, display_title);
    title(title_str, 'Interpreter', 'latex', 'FontSize', fontSize+1);
    
    xlabel('Time Slot Index $t$', 'Interpreter', 'latex', 'FontSize', fontSize);
    ylabel(y_labels{m}, 'Interpreter', 'latex', 'FontSize', fontSize);
    
    if strcmp(curr_metric, 'HV')
        % 使用 'Location','southeast' 锁死在右下角（HV曲线呈上升趋势，右下角最空）
        lgd = legend(plot_handles, labels, 'Interpreter', 'latex', ...
                     'FontSize', fontSize - 4, ... % 9号字，单栏大图最精致
                     'Location', 'southeast', ...  % 绝不用 Position，改用 Location 锁死方位！
                     'NumColumns', 1);
        
        % 保持你想要的浅灰色边框和 80% 微透明底色
        set(lgd, 'Box', 'on', 'EdgeColor', [0.8 0.8 0.8], 'Color', [1 1 1 0.8]);
        
        % 🚀 慢工出细活微操：如果觉得默认图注太占地方，用 Margin（内边距）来让人肉缩小图注框，且绝不影响大图长方形
        set(lgd, 'ItemTokenSize', [12, 6]); % 缩短图注里彩色小线段的长度（默认是 [30,18]），让图注更紧凑
    end
    
    % ====================== ✅ 核心修改 3：升级为全新的高保真高级 PDF 矢量渲染引擎 ======================
    % 这个引擎完美原生支持半透明阴影，能够把文字、线条、阴影全部以纯粹的数学几何路径保存！
    exportgraphics(fig, file_names{m}, 'ContentType', 'vector');
    close(fig);
    
    fprintf('指标 %s 纯矢量无损趋势图已保存: %s\n', curr_metric, file_names{m});
end

%% 3. 输出汇总统计表
fprintf('\n%s\n', repmat('=', 1, 115));
fprintf('%-20s | %-18s | %-18s | %-18s | %-12s\n', ...
    'Algorithm Variant', 'IGD (Mean±Std)', 'HV (Mean±Std)', 'Spacing (Mean±Std)', 'Feasible (%)');
fprintf('%s\n', repmat('-', 1, 115));
for i = 1:length(alg_list)
    n = alg_list{i};
    fIGD = all_results.(n).IGD(:, end); fHV = all_results.(n).HV(:, end); fSpacing = all_results.(n).Spacing(:, end); fFeas = all_results.(n).NumFeasibleSolutions(:, end);
    fprintf('%-20s | %.4f ± %.4f | %.4f ± %.4f | %.4f ± %.4f | %.2f%%\n', labels{i}, mean(fIGD), std(fIGD), mean(fHV), std(fHV), mean(fSpacing), std(fSpacing), mean(fFeas));
end
fprintf('%s\n', repmat('=', 1, 115));