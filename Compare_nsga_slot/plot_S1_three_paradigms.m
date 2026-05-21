% =========================================================================
% 论文发表级：Scenario S1 三大范式格式【绝对统一】看板 (Transactions 定稿版)
% 描述: 载入 S1.mat 和 精确解.mat，提取 Run 27 数据。
% 应对审稿人: 
%   1. 解决 "small fonts": 引入相对物理尺寸 (inches) 与 fontSize 分级机制。
%   2. 解决 "dense visualization": 采用纯矢量 PDF 导出，无限放大不模糊。
%   3. 解决 "missing legends": 全局 LaTeX 图注强制置顶。
% =========================================================================
clc; clear; close all;

%% --- 步骤 0: 全局环境与 LaTeX 设置 (严格对齐消融实验框架) ---
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

% 基础路径配置
base_folder = "E:\guthub-matlab";
s1_path     = fullfile(base_folder, "ExperimentResults_Output\S1.mat");
exact_path  = fullfile(base_folder, "Exact_Pruned_Baselines_20260517_2109.mat");

algonames_data   = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III', 'DRL_Baseline'};
algonames_legend = {'Ours (HMD-NSGA-II)', '', '', '', '', 'MODDPG'}; % 严格对齐正文 MODDPG
final_slot_idx   = 10; 
target_run_idx   = 27; % 严格锁定中位数运行 Run 27

% =========================================================================
% 🎨 颜色方案 (紫红、纯黑、棕红)
% =========================================================================
Color_Ours  = [0.70, 0.10, 0.55]; % 紫红色 (Ours)
Color_Exact = [0.00, 0.00, 0.00]; % 纯黑 (HAES)
Color_DRL   = [0.60, 0.20, 0.20]; % 棕红 (MODDPG)

colors = [Color_Ours; [0,0,0]; [0,0,0]; [0,0,0]; [0,0,0]; Color_DRL];

%% --- 步骤 1: 数据加载 ---
if ~isfile(s1_path), error('❌ 未找到 S1 结果文件: %s', s1_path); end
if ~isfile(exact_path), error('❌ 未找到精确解结果文件: %s', exact_path); end

data_struct = load(s1_path);
exact_struct = load(exact_path);
exact_pf = exact_struct.all_scenarios_pruned_exact_pfs{1}{target_run_idx};

%% --- 步骤 2: 核心看板绘制 (引入 Inches 物理尺寸) ---
% 🌟 采用 4.5 x 3.5 英寸物理画幅，确保插入 IEEE 双栏排版时字体依然硕大清晰
fig = figure('Color', 'w', 'Units', 'inches', 'Position', [2 2 4.5 3.5], 'Visible', 'on'); 
hold on; grid on; box on;
all_visible_points = []; 
plot_handles = []; % 用于精确控制图注顺序

% -------------------------------------------------------------------------
% 1. 🌟 传统剪枝精确求解器前沿 (Exact HAES) -> 黑色、空心
% -------------------------------------------------------------------------
if ~isempty(exact_pf)
    [~, sort_idx_exact] = sort(exact_pf(:,1));
    exact_pf = exact_pf(sort_idx_exact, :);
    all_visible_points = [all_visible_points; exact_pf];
    
    % 纯黑虚线 + 纯白底空心圆圈
    h_exact = plot(exact_pf(:,1), exact_pf(:,2), '--o', 'Color', Color_Exact, ...
        'LineWidth', 1.3, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', Color_Exact, ...
        'DisplayName', 'HAES');
    plot_handles = [plot_handles, h_exact];
    
    % 空心五角星拐点
    kp_exact = find_knee_point(exact_pf);
    if ~isempty(kp_exact)
        plot(kp_exact(1), kp_exact(2), 'p', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', Color_Exact, ...
            'LineWidth', 1.2, 'HandleVisibility', 'off'); 
    end
end

% -------------------------------------------------------------------------
% 2. 🌟 定向绘制 Ours 和 MODDPG -> 严格采用空心统一逻辑
% -------------------------------------------------------------------------
all_fronts_data = data_struct.current_scenario_all_run_final_fronts;
target_indices = [1, 6]; 

for i = target_indices
    front = all_fronts_data{i, final_slot_idx, target_run_idx};
    if isempty(front), continue; end
    front(any(isnan(front) | isinf(front), 2), :) = [];
    
    [~, sort_idx] = sort(front(:,1));
    front = front(sort_idx, :);
    all_visible_points = [all_visible_points; front];
    
    % 虚线(--) + 圆圈(o) + 空心(MarkerFaceColor='w')
    h_alg = plot(front(:,1), front(:,2), '--o', 'Color', colors(i,:), ...
        'LineWidth', 1.3, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(i,:), ...
        'DisplayName', algonames_legend{i});
    plot_handles = [plot_handles, h_alg];
    
    % 空心五角星拐点
    kp = find_knee_point(front);
    if ~isempty(kp)
        plot(kp(1), kp(2), 'p', 'MarkerSize', 10, ...
            'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(i,:), ...
            'LineWidth', 1.2, 'HandleVisibility', 'off'); 
    end
end

%% --- 步骤 3: 坐标轴自动适配与学术美化 (严格字号管理) ---
if ~isempty(all_visible_points)
    x_min = min(all_visible_points(:,1)); x_max = max(all_visible_points(:,1));
    y_min = min(all_visible_points(:,2)); y_max = max(all_visible_points(:,2));
    dx = (x_max - x_min) * 0.08; dy = (y_max - y_min) * 0.08;
    xlim([x_min - dx, x_max + dx]); ylim([y_min - dy, y_max + dy]);
end

% 🌟 动态字号基准 (应对 reviewer "small fonts" 质疑)
fontSize = 11; 

set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on', 'LineWidth', 1.1, ...
    'FontSize', fontSize);
xlabel('$\mathrm{G}_1$ (Average Latency / s)', 'Interpreter', 'latex', 'FontSize', fontSize);
ylabel('$\mathrm{G}_2$ (Total Energy Consumption / J)', 'Interpreter', 'latex', 'FontSize', fontSize);
title('Pareto Frontiers Juxtaposition', 'Interpreter', 'latex', 'FontSize', fontSize + 1); % 简化标题避免拥挤

% 🌟 图注：设定半透明白底避免遮挡，右上角排列
lgd = legend(plot_handles, {'HAES', 'Ours (HMD-NSGA-II)', 'MODDPG'}, ...
             'Interpreter', 'latex', 'FontSize', fontSize - 4, ...
             'Location', 'northeast', 'NumColumns', 1);
set(lgd, 'Box', 'on', 'EdgeColor', [0.8 0.8 0.8], 'Color', [1 1 1 0.8], 'ItemTokenSize', [15, 6]);

%% --- 步骤 4: 高保真纯矢量 PDF 导出 (无限放大不模糊) ---
save_folder = 'E:\guthub-matlab\第二篇论文数据';
if ~exist(save_folder, 'dir'), mkdir(save_folder); end

% 🌟 强制保存为 PDF 格式
full_save_path = fullfile(save_folder, 'PF_Result_S1.pdf');

% 🌟 调用全新纯矢量引擎
exportgraphics(fig, full_save_path, 'ContentType', 'vector');
fprintf('\n📊 【矢量定稿】完美统一格式、字号达标的 PDF 纯矢量图已导出至：\n  %s\n', full_save_path);

%% --- 辅助函数：寻找拐点 ---
function knee_point = find_knee_point(front)
    if isempty(front) || size(front, 1) < 3, knee_point = []; return; end
    f_min = min(front); f_max = max(front);
    range = f_max - f_min; range(range == 0) = 1e-6;
    norm_f = (front - f_min) ./ range;
    [~, i1] = min(norm_f(:,1)); [~, i2] = min(norm_f(:,2));
    p1 = norm_f(i1,:); p2 = norm_f(i2,:);
    vec = p2 - p1;
    dist = abs(vec(1)*(norm_f(:,2)-p1(2)) - vec(2)*(norm_f(:,1)-p1(1))) / max(norm(vec), 1e-6);
    [~, k_idx] = max(dist);
    knee_point = front(k_idx, :);
end