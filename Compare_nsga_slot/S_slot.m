clc; clear; close all;

% --- 配置路径 ---
base_path = 'E:\guthub-matlab\第二篇论文数据'; 
output_dir = fullfile(base_path, 'output_figures');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

file_names_list = {
    '30-5-20-2-100-100.mat',
    '30-5-20-200-200.mat',
    '30-5-45-100-100.mat',
    '30-5-45-200-200.mat',
    '30-5-80-100-100.mat',
    '30-5-80-200-200.mat'
};

algo_name_to_plot = 'Ours';
algo_field_name = 'MyNSGA_II';

metrics_to_plot = {'IGD', 'HV'};
metric_colors = [
    0, 1, 1;   % Cyan - IGD
    1, 0, 1    % Magenta - HV
];

slots_to_show = 5;
num_runs = 30;
offset = 0.12;

% --- 遍历六个场景 ---
for file_idx = 1:length(file_names_list)
    file_to_analyze = fullfile(base_path, file_names_list{file_idx});
    if ~isfile(file_to_analyze)
        warning('文件缺失: %s，跳过', file_to_analyze);
        continue;
    end

    data = load(file_to_analyze);
    result = data.all_scenario_results;

    % 检查字段
    for i = 1:2
        m = metrics_to_plot{i};
        if ~isfield(result.(algo_field_name), m)
            warning('缺少指标 "%s"，跳过', m);
            continue;
        end
    end

    igd_data = result.(algo_field_name).IGD{1}(:, 1:slots_to_show);
    hv_data  = result.(algo_field_name).HV{1}(:, 1:slots_to_show);

    % 构建绘图数据
    positions = [];
    plot_data = [];
    group = [];
    slot_labels = {};
    for s = 1:slots_to_show
        pos_igd = s - offset;
        pos_hv  = s + offset;

        plot_data = [plot_data; igd_data(:, s); hv_data(:, s)];
        group = [group; repmat(pos_igd, num_runs, 1); repmat(pos_hv, num_runs, 1)];
        positions = [positions, pos_igd, pos_hv];

        slot_labels{end+1} = sprintf('Solt %d', s);
    end

    % 绘图
    figure('Position', [100, 100, 900, 450]);
    hold on;
    boxplot(plot_data, group, 'Positions', positions, 'Symbol', '+', 'Widths', 0.2);

    % 上色
    h = findobj(gca,'Tag','Box');
    for j = 1:length(h)
        if mod(j, 2) == 1
            cidx = 2;
        else
            cidx = 1;
        end
        patch(get(h(j),'XData'), get(h(j),'YData'), metric_colors(cidx,:), 'FaceAlpha', 0.6);
    end

    set(gca, 'XTick', 1:slots_to_show, 'XTickLabel', slot_labels, 'FontSize', 12);
    ylabel('Value', 'FontSize', 14);

    % ✅ 自动设置 Y 轴范围
    y_max = max(plot_data, [], 'omitnan');
    y_min = min(plot_data, [], 'omitnan');
    y_margin = 0.05 * (y_max - y_min);
    ylim([max(0, y_min - y_margin), y_max + y_margin]);

    title(sprintf('S%d', file_idx));


    % 图例
    h_legend(1) = patch(NaN, NaN, metric_colors(1,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    h_legend(2) = patch(NaN, NaN, metric_colors(2,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k');
    legend(h_legend, metrics_to_plot, 'Location', 'northwest', 'FontSize', 12);

    grid on; box on; hold off;

    % 可选：保存图像（如有需要）
    % saveas(gcf, fullfile(output_dir, sprintf('SlotPlot_S%d.png', file_idx)));

end

