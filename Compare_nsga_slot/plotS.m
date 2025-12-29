% =========================================================================
% 场景代表性帕累托前沿分析脚本 (S1-S6 自动批量生成版)
% 功能：自动循环 6 个 .mat 文件，匹配 IGD 中位数运行，高清导出
% =========================================================================
clc; clear; close all;

%% --- 步骤 0: 配置信息 ---
% 基础路径（根据你的实际路径调整前缀）
base_folder = "E:\guthub-matlab\第二篇论文数据\";
algonames_legend = {'Ours', 'NSGA-II', 'MOEA/D', 'INSGA-II', 'NSGA-III'};
algonames_data = {'MyNSGA_II', 'NSGA_II', 'MOEA_D', 'INSGA_II', 'NSGA_III'};

final_slot_idx = 10; % 固定选择第 10 时隙

% 学术配色方案
colors = [
    1, 0.5, 0;      % 橙色 (Ours)
    0, 0.4, 0.8;    % 蓝色 (NSGA-II)
    0, 0.8, 0.1;    % 绿色 (MOEA/D)
    0, 1, 1;        % 青色 (INSGA-II)
    1, 0, 1;        % 紫色 (NSGA-III)
];

%% --- 步骤 1: 循环处理 S1 到 S6 ---
for s_idx = 1:6
    % 构建当前场景的文件完整路径
    file_name = sprintf("S%d_V5.mat", s_idx);
    full_path = fullfile(base_folder, file_name);
    
    if ~isfile(full_path)
        fprintf('跳过：文件 %s 不存在。\n', full_path);
        continue;
    end
    
    fprintf('正在处理场景 S%d (文件: %s)...\n', s_idx, file_name);
    data_struct = load(full_path);
    
    % --- 1.1 寻找该文件中的代表性运行 (HV 中位数) ---
    num_algs = length(algonames_data);
    representative_runs = zeros(num_algs, 1);
    for i = 1:num_algs
        alg_field = algonames_data{i};
        % 注意：此处假设每个文件中 s_idx 始终对应当前的场景索引或为 1
        % 如果每个文件内只存了该场景的数据，可能需要将 s_idx 设为 1
        target_s = s_idx; 
        if ~isfield(data_struct.all_scenario_results.(alg_field).HV, 's' + string(target_s)) && ...
           iscell(data_struct.all_scenario_results.(alg_field).HV)
           % 自动适配 cell 格式或其它存储格式
        end
        
        try
            igd_list = data_struct.all_scenario_results.(alg_field).IGD{s_idx}(:, final_slot_idx);
            med_igd = median(igd_list, 'omitnan');
            [~, run_idx] = min(abs(igd_list - med_igd));
            representative_runs(i) = run_idx;
        catch
            representative_runs(i) = 1; 
        end
    end

    % --- 1.2 绘图与美化 ---
    fig = figure('Color', 'w', 'Position', [100 100 800 600], 'Visible', 'off'); % 后台运行，不弹窗
    hold on; grid on; box on;
    all_visible_points = []; 

    % 1. 绘制参考前沿 (Reference PF)
    if isfield(data_struct, 'igd_reference_front_obj')
        pf_star = data_struct.igd_reference_front_obj;
        
        % 定义一种学术绯红色
        academic_red = [1, 0.1, 0.1]; 
        
        % 使用单层 scatter，设置较低的 Alpha 值达到“透”的效果
        scatter(pf_star(:,1), pf_star(:,2), 35, ... 
            'Marker', 'o', ...
            'MarkerFaceColor', academic_red, ...
            'MarkerEdgeColor', [0.2, 0.2, 0.2], ...
            'MarkerFaceAlpha', 0.5, ...           % 透明度设为 0.35，呈现透亮的质感
            'DisplayName', 'Reference PF');
            
        all_visible_points = [all_visible_points; pf_star];
    end
    % 2. 遍历并绘制各算法前沿
    all_fronts_data = data_struct.current_scenario_all_run_final_fronts;
    for i = 1:num_algs
        run_idx = representative_runs(i);
        % 假设数据结构为 {alg_idx, slot_idx, run_idx}
        front = all_fronts_data{i, final_slot_idx, run_idx};
        
        if isempty(front), continue; end
        front(any(isnan(front) | isinf(front), 2), :) = [];
        
        scatter(front(:,1), front(:,2), 65, ...
            'Marker', 'o', 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', [0.15 0.15 0.15], 'MarkerFaceAlpha', 0.7, ...
            'LineWidth', 0.5, 'DisplayName', algonames_legend{i});
        
        all_visible_points = [all_visible_points; front];
        
        % 绘制拐点 (Knee Point)
        kp = find_knee_point(front);
        if ~isempty(kp)
            plot(kp(1), kp(2), 'p', 'MarkerSize', 15, ...
                'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.0, 'HandleVisibility', 'off'); 
        end

        % 绘制并打印拐点 (Knee Point)
        if ~isempty(kp)
            % 【新增】：在此处打印拐点坐标信息
            fprintf('  - %-10s 的拐点坐标: [G1: %.4f, G2: %.4f]\n', algonames_legend{i}, kp(1), kp(2));
            
            plot(kp(1), kp(2), 'p', 'MarkerSize', 15, ...
                'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', ...
                'LineWidth', 1.0, 'HandleVisibility', 'off'); 
        end
    end

    % --- 1.3 坐标轴自动适配 ---
    if ~isempty(all_visible_points)
        x_min = min(all_visible_points(:,1)); x_max = max(all_visible_points(:,1));
        y_min = min(all_visible_points(:,2)); y_max = max(all_visible_points(:,2));
        dx = (x_max - x_min) * 0.08; dy = (y_max - y_min) * 0.08;
        xlim([x_min - dx, x_max + dx]); ylim([y_min - dy, y_max + dy]);
    end

    set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'on', 'LineWidth', 1.2, ...
        'FontName', 'Times New Roman', 'FontSize', 12);
    xlabel('$\mathrm{G}_1$ (Latency/s)', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel('$\mathrm{G}_2$ (Energy/J)', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('Pareto Front Comparison - S%d', s_idx), 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);

    % 【核心修改】：仅 S1 显示图注
    if s_idx == 1
        legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 10, 'EdgeColor', 'none');
    else
        legend('off');
    end

    % --- 1.4 导出图像 ---
    save_name = sprintf('PF_Result_S%d.png', s_idx);
    exportgraphics(fig, save_name, 'Resolution', 300);
    fprintf('成功保存： %s\n', save_name);
    
    close(fig); % 及时关闭窗口释放内存
end

fprintf('--- 所有场景处理完毕 ---\n');

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