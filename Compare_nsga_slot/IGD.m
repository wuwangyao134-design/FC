clc; clear; close all;

% --- 配置区域 (更易于维护) ---
% 使用相对路径或让用户选择
base_path = 'E:\guthub-matlab\第二篇论文数据'; 

file_names = {
    '30-5-20-2-100-100.mat'
    '30-5-20-200-200.mat'
    '30-5-45-100-100.mat'
    '30-5-45-200-200.mat'
    '30-5-80-100-100.mat'
    '30-5-80-200-200.mat'
};

% 下面这行现在可以正常工作了
file_list = fullfile(base_path, file_names); 

% 确保图例名称和数据字段名一一对应
algonames_legend = {'Ours', 'NSGA-II', 'A-MOPSO', 'ARG'}; % 用于图例的名称
algonames_data = {'MyNSGA_II', 'NSGA_II', 'MOPSO', 'Baseline'}; % 用于访问数据的字段名

% 检查配置是否一致
if length(algonames_legend) ~= length(algonames_data)
    error('图例名称和数据字段名的数量不匹配!');
end

num_scenarios = length(file_list);
num_algorithms = length(algonames_legend);
num_runs = 30; % 假设每次运行都是30次
slot_to_plot = 5; % 将要绘制的列明确定义为变量

% --- 初始化 ---
all_data = [];
group = [];
positions = [];

colors = [
    0, 1, 1;       % cyan - Ours
    1, 0, 1;       % magenta - NSGA-II
    1, 0.5, 0;     % orange - A-MOPSO
    1, 1, 0        % yellow - ARG
];
% --- 数据处理循环 ---
for s = 1:num_scenarios
    % 健壮性检查：文件是否存在
    if ~isfile(file_list{s})
        warning('文件不存在，已跳过: %s', file_list{s});
        continue;
    end
    data = load(file_list{s});
    
    % 健壮性检查：变量是否存在
    if ~isfield(data, 'all_scenario_results')
        warning('在文件 %s 中未找到变量 "all_scenario_results"，已跳过。', file_list{s});
        continue;
    end
    result = data.all_scenario_results;

    % 提取当前场景的所有算法数据
    scenario_igd_data = zeros(num_runs, num_algorithms);
    for a = 1:num_algorithms
        field_name = algonames_data{a};
        % 健壮性检查：字段和数据结构是否正确
        if isfield(result, field_name) && isfield(result.(field_name), 'NumFeasibleSolutions') && ...
           iscell(result.(field_name).NumFeasibleSolutions) && ~isempty(result.(field_name).NumFeasibleSolutions) && ...
           size(result.(field_name).NumFeasibleSolutions{1}, 2) >= slot_to_plot
            
            current_data = result.(field_name).NumFeasibleSolutions{1}(:, slot_to_plot);
            % 检查运行次数是否匹配
            if size(current_data, 1) ~= num_runs
                warning('文件 %s 中算法 %s 的运行次数 (%d) 与预设值 (%d) 不符', ...
                        file_list{s}, field_name, size(current_data, 1), num_runs);
            end
            scenario_igd_data(:, a) = current_data;
        else
            warning('在文件 %s 中算法 %s 的数据结构不完整，将使用 NaN 填充。', file_list{s}, field_name);
            scenario_igd_data(:, a) = NaN(num_runs, 1);
        end
    end
    
    % 构建 boxplot 所需的向量
    for a = 1:num_algorithms
        % 定义位置：每个场景占一组，组内算法有小的偏移
        current_pos = s + (a - (num_algorithms + 1) / 2) * 0.18; % 居中对齐
        positions = [positions; current_pos]; % 正确构建 positions 向量
        all_data = [all_data; scenario_igd_data(:, a)];
        group = [group; repmat(current_pos, num_runs, 1)];
    end
end

% --- 绘图 ---
figure('Position', [100, 100, 1200, 500]);
hold on;

% 绘制 boxplot
boxplot(all_data, group, 'positions', positions, ...
        'colors', 'k', 'symbol', '+', 'widths', 0.15);

% 使用更可靠的方法为箱线图上色
h = findobj(gca,'Tag','Box');
% findobj 返回的句柄顺序通常是逆向的，所以我们倒序上色
color_indices = repmat(1:num_algorithms, 1, num_scenarios);
if length(h) ~= length(color_indices)
    warning('箱线图数量与预期不符，颜色可能不正确。');
else
    for j = 1:length(h)
        patch(get(h(j),'XData'), get(h(j),'YData'), colors(color_indices(length(h)-j+1),:), 'FaceAlpha', 0.6);
    end
end

% --- 格式化图表 ---
set(gca, 'XTick', 1:num_scenarios, 'XTickLabel', {'S1','S2','S3','S4','S5','S6'}, 'FontSize', 12);
ylabel('Spacing', 'FontSize', 14);
title(['Spacing Boxplot (Slot ', num2str(slot_to_plot), ')'], 'FontSize', 16, 'FontWeight', 'bold');
grid on;
box on;

% 创建自定义图例
h_legend = gobjects(num_algorithms, 1);
for i = 1:num_algorithms
    h_legend(i) = patch(NaN, NaN, colors(i,:), 'FaceAlpha', 0.6, 'EdgeColor', 'k');
end
legend(h_legend, algonames_legend, 'Location', 'northwest', 'FontSize', 12);

hold off;