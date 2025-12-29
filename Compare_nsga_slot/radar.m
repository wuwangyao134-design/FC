% =========================================================================
% 顶级现代几何雷达图 (Contrast Enhanced Academic Version)
% 优化：强化刻度/标签对比度，全图 LaTeX 深度渲染，精修 Ours 特效
% =========================================================================
clc; clear; close all;

%% --- 步骤 0: 配置与数据提取 ---
base_folder = "E:\guthub-matlab\第二篇论文数据\";
alg_legend = {'Ours', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III'};
alg_fields = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III'};
target_slot = 10;

% 算法配色：提高饱和度以增强区分
colors = [
    0.95, 0.40, 0.10;  % 炽橙 (Ours)
    0.05, 0.40, 0.85;  % 宝石蓝
    0.10, 0.60, 0.10;  % 深草绿
    0.00, 0.75, 0.75;  % 暗青
    0.70, 0.00, 0.70;  % 紫罗兰
];

hv_matrix = zeros(length(alg_fields), 6); 
for s_idx = 1:6
    file_path = fullfile(base_folder, sprintf('S%d_V5.mat', s_idx));
    if ~isfile(file_path), continue; end
    data_struct = load(file_path);
    for a = 1:length(alg_fields)
        try
            metric_cell = data_struct.all_scenario_results.(alg_fields{a}).Spacing;
            idx = min(s_idx, length(metric_cell));
            run_data = metric_cell{idx}(:, target_slot);
            hv_matrix(a, s_idx) = median(run_data, 'omitnan');
        catch, hv_matrix(a, s_idx) = 0; end
    end
end

% 归一化缩放：提高中心起始值 (0.35) 确保中心区域文字清晰
norm_hv = (hv_matrix - min(hv_matrix)) ./ (max(hv_matrix) - min(hv_matrix) + 1e-6);
norm_hv = norm_hv * 0.65 + 0.35; 

%% --- 步骤 1: 绘图初始化 ---
fig = figure('Color', 'w', 'Position', [100, 100, 800, 850]);
% 调整绘图区位置，留足外部标签空间
ax = axes('Position', [0.12, 0.15, 0.76, 0.75]); 
hold on; axis equal; axis off;

num_scenes = 6;
num_algs = length(alg_fields);
scene_angles = linspace(0, 2*pi, num_scenes + 1) + pi/2;
scene_angles(end) = [];

% 1. 绘制背景圆环 (强化颜色：[0.7 0.7 0.7])
r_ticks = 0.4:0.2:1.0; % 减少冗余刻度，强化关键刻度
theta_fine = linspace(0, 2*pi, 500);
for r = r_ticks
    plot(r*cos(theta_fine), r*sin(theta_fine), 'Color', [0.3 0.3 0.3], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    % 刻度数字：颜色加深至 [0.3 0.3 0.3]，加粗
    text(0.02, r, sprintf('$%.1f$', r), 'Interpreter', 'latex', 'FontSize', 11, ...
        'Color', [0.1 0.1 0.1], 'FontWeight', 'bold');
end

% 2. 绘制放射轴线 (强化颜色：[0.8 0.8 0.8])
for k = 1:num_scenes
    line([0.1*cos(scene_angles(k)), 1.05*cos(scene_angles(k))], ...
         [0.1*sin(scene_angles(k)), 1.05*sin(scene_angles(k))], ...
         'Color', [0.1 0.1 0.1], 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

%% --- 步骤 2: 绘制算法包络 (强化对比版) ---
h_plots = gobjects(num_algs, 1);

for a = 1:num_algs
    r_vals = norm_hv(a, :);
    x_plot = [r_vals .* cos(scene_angles), r_vals(1) * cos(scene_angles(1))];
    y_plot = [r_vals .* sin(scene_angles), r_vals(1) * sin(scene_angles(1))];
    
    if a == 1 % --- Ours: 核心强调，增加“实体感” ---
        % 填充背景，增加不透明度 (0.12)
        patch(x_plot, y_plot, colors(a,:), 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        % 外发光描边 (多层渲染)
        for w = 1:2
            plot(x_plot, y_plot, 'Color', [colors(a,:) 0.15/w], 'LineWidth', 4 + w*4, 'HandleVisibility', 'off');
        end
        % 加粗主线 (LineWidth: 3.5)
        h_plots(a) = plot(x_plot, y_plot, 'Color', colors(a,:), 'LineWidth', 3.5, ...
            'Marker', 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(a,:));
    else % --- 其他算法: 极简虚线，腾出视觉空间 ---
        patch(x_plot, y_plot, colors(a,:), 'FaceAlpha', 0.04, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        % 使用细实线，降低透明度，区分 Ours
        h_plots(a) = plot(x_plot, y_plot, 'Color', [colors(a,:) 0.6], 'LineWidth', 1.5, 'Marker', 'none');
    end
end

%% --- 步骤 3: 标签与标题 (极黑字体，增强清晰度) ---
for k = 1:num_scenes
    % 使用 mathit 渲染，颜色改为纯黑 [0 0 0]，字号加大到 14
    lbl = sprintf('$\\mathit{Scenario\\ S%d}$', k); 
    tx = 1.25 * cos(scene_angles(k));
    ty = 1.25 * sin(scene_angles(k));
    text(tx, ty, lbl, 'Interpreter', 'latex', 'FontSize', 14, ...
        'Color', [0 0 0], 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

% 标题：加深加粗
%title('$\mathbf{Global\ Robustness\ Analysis\ (Normalized\ HV)}$', ...
  %  'Interpreter', 'latex', 'FontSize', 20, 'Units', 'normalized', 'Position', [0.5, 1.08], 'Color', [0 0 0]);

% 图例：移除边框，增加间距
lgd = legend(h_plots, alg_legend, 'Location', 'southoutside', ...
    'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 13, 'Interpreter', 'latex');
lgd.Units = 'normalized';
lgd.Position(2) = 0.03; 

% 导出 (600 DPI 保证印刷级清晰)
exportgraphics(fig, 'Ultra_Contrast_Radar_Spacing.png', 'Resolution', 600);