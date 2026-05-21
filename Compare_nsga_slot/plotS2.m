% =========================================================================
% 论文发表级：Scenario S1, S3, S5, S6, S7, S8 前沿批量导出 (带 a-f 自动编号)
% =========================================================================
clc; clear; close all;

%% --- 步骤 0: 全局环境与 LaTeX 引擎设置 ---
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');

% 基础路径配置
base_folder = "E:\guthub-matlab\ExperimentResults_Output";
save_folder = "E:\guthub-matlab\第二篇论文数据";
if ~exist(save_folder, 'dir'), mkdir(save_folder); end

% 算法配置
algonames_legend = {'Ours (HMD-NSGA-II)', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III', 'MODDPG'};
algonames_data   = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III', 'DRL_Baseline'};
final_slot_idx   = 10; 

% =========================================================================
% 🎨 学术配色方案
% =========================================================================
colors = [
    0.70, 0.10, 0.55; % 紫红色 (Ours)
    0.00, 0.40, 0.80; % 蓝色 (NSGA-II)
    0.00, 0.80, 0.10; % 绿色 (MOEA/D)
    0.00, 0.55, 0.55; % 孔雀蓝/深青色 (INSGA-II)
    1.00, 0.00, 1.00; % 紫色 (NSGA-III)
    0.60, 0.20, 0.20; % 棕红色 (DRL-Baseline)
];
% 全局字号基准
fontSize = 11; 

%% --- 步骤 1: 锁定目标场景与自动编号 ---
target_scenarios = [1, 3, 5, 6, 7, 8];
sub_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};

for loop_idx = 1:length(target_scenarios)
    s_idx = target_scenarios(loop_idx);
    current_label = sub_labels{loop_idx};
    
    file_name = sprintf("S%d.mat", s_idx);
    full_path = fullfile(base_folder, file_name);
    
    if ~isfile(full_path)
        fprintf('⚠️ 跳过：文件 %s 不存在。\n', full_path);
        continue;
    end
    
    fprintf('▶️ 正在处理场景 %s S%d (文件: %s)...\n', current_label, s_idx, file_name);
    data_struct = load(full_path);
    
    % 动态寻找数据真实的存储索引
    if isfield(data_struct, 'all_reference_fronts') && length(data_struct.all_reference_fronts) >= s_idx && ~isempty(data_struct.all_reference_fronts{s_idx})
        data_s_idx = s_idx; 
    else
        data_s_idx = 1;     
    end
    
    % --- 1.1 寻找各个算法的代表性运行 (依据 HV 中位数) ---
    num_algs = length(algonames_data);
    representative_runs = zeros(num_algs, 1);
    
    for i = 1:num_algs
        alg_field = algonames_data{i};
        try
            hv_list = data_struct.all_scenario_results.(alg_field).HV{data_s_idx}(:, final_slot_idx);
            med_hv = median(hv_list, 'omitnan');
            [~, run_idx] = min(abs(hv_list - med_hv));
            representative_runs(i) = run_idx;
        catch
            representative_runs(i) = 1; 
        end
    end
    
    target_run_idx = representative_runs(1); 
    fprintf('  - 锁定代表性运行 (环境): Run %d\n', target_run_idx);
    
    % --- 1.2 核心看板绘制 ---
    fig = figure('Color', 'w', 'Units', 'inches', 'Position', [2 2 4.5 3.5], 'Visible', 'off'); 
    hold on; grid on; box on;
    all_visible_points = []; 
    plot_handles = []; 
    lgd_labels = {};
    
    % 🌟 绘制参考前沿 (Reference PF*)
    if isfield(data_struct, 'all_reference_fronts') && length(data_struct.all_reference_fronts) >= data_s_idx
        pf_star = data_struct.all_reference_fronts{data_s_idx}{target_run_idx, final_slot_idx};
        if ~isempty(pf_star)
            academic_red = [0.8, 0.1, 0.1]; 
            h_ref = plot(pf_star(:,1), pf_star(:,2), '-', 'Color', academic_red, ...
                'LineWidth', 1.5, 'DisplayName', 'Reference $\mathrm{PF}^*$');
            
            % 半透明空心套靶圆圈
            scatter(pf_star(:,1), pf_star(:,2), 60, 'Marker', 'o', ...
                'MarkerEdgeColor', academic_red, 'MarkerEdgeAlpha', 0.5, ...
                'MarkerFaceColor', 'none', ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
                
            all_visible_points = [all_visible_points; pf_star];
            plot_handles = [plot_handles, h_ref];
            lgd_labels{end+1} = 'Reference $\mathrm{PF}^*$';
        end
    end
    
    % 🌟 遍历并绘制各算法前沿
    all_fronts_data = data_struct.current_scenario_all_run_final_fronts;
    for i = 1:num_algs
        front = all_fronts_data{i, final_slot_idx, target_run_idx};
        if isempty(front), continue; end
        front(any(isnan(front) | isinf(front), 2), :) = [];
        
        [~, sort_idx] = sort(front(:,1));
        front = front(sort_idx, :);
        all_visible_points = [all_visible_points; front];
        
        h_alg = plot(front(:,1), front(:,2), '--o', 'Color', colors(i,:), ...
            'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(i,:), ...
            'DisplayName', algonames_legend{i});
            
        plot_handles = [plot_handles, h_alg];
        lgd_labels{end+1} = algonames_legend{i};
        
        kp = find_knee_point(front);
        if ~isempty(kp)
            plot(kp(1), kp(2), 'p', 'MarkerSize', 10, ...
                'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(i,:), ...
                'LineWidth', 1.2, 'HandleVisibility', 'off'); 
        end
    end
    
    % --- 1.3 坐标轴自动适配与学术美化 ---
    if ~isempty(all_visible_points)
        x_min = min(all_visible_points(:,1)); x_max = max(all_visible_points(:,1));
        y_min = min(all_visible_points(:,2)); y_max = max(all_visible_points(:,2));
        dx = (x_max - x_min) * 0.08; dy = (y_max - y_min) * 0.08;
        xlim([x_min - dx, x_max + dx]); ylim([y_min - dy, y_max + dy]);
    end
    
    set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on', 'LineWidth', 1.1, ...
        'FontSize', fontSize);
    xlabel('$\mathrm{G}_1$ (Latency / s)', 'Interpreter', 'latex', 'FontSize', fontSize);
    ylabel('$\mathrm{G}_2$ (Energy / J)', 'Interpreter', 'latex', 'FontSize', fontSize);
    
    % 🌟 【核心修改】：自带序号的 Title，例如 "(a) Pareto Frontiers Comparison (S1)"
    title(sprintf('%s Pareto Frontiers Comparison (S%d)', current_label, s_idx), ...
          'Interpreter', 'latex', 'FontSize', fontSize + 1);
    
    % 🌟 统一图注管理 (只在 S1 显示)
    if s_idx == 1
        lgd = legend(plot_handles, lgd_labels, 'Interpreter', 'latex', ...
            'FontSize', fontSize - 4, 'Location', 'northeast'); 
        set(lgd, 'Box', 'on', 'EdgeColor', [0.8 0.8 0.8], 'Color', [1 1 1 0.8], 'ItemTokenSize', [15, 6]);
    end
    
    % --- 1.4 高保真纯矢量 PDF 导出 ---
    save_name = sprintf('PF_Result_S%d.pdf', s_idx);
    full_save_path = fullfile(save_folder, save_name);
    
    exportgraphics(fig, full_save_path, 'ContentType', 'vector');
    fprintf('✅ 成功导出带编号纯矢量 PDF： %s\n\n', full_save_path);
    
    close(fig); 
end
fprintf('🎉 --- 6个精选场景图片批量导出完毕 (自带 LaTeX 编号) ---\n');

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