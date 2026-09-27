%% Summarize_S1_S8_Metrics
% 统计 S1--S8 第10时隙的 IGD、HV、FFR、Tmax、Nnds 和 Runtime。
% 输出：
%   1) S1_S8_Metric_Summary.csv      长表，可用于复核和二次制表；
%   2) S1_S8_Metric_Summary.mat      MATLAB汇总结构；
%   3) S1_S8_LaTeX_Table.tex         可直接复制的LaTeX表格。
%
% [统计修改-01 | 2026-09-24]
% 原因：汇总最终S1--S8的30次独立运行结果，以FFR替换Spacing，并正确
% 处理极端场景中空前沿导致的NaN。场景由I/M/区域尺寸识别，不依赖旧文件
% 中可能不准确的S编号。
%
% [统计修改-02 | 2026-09-24]
% 原因：正式数据已整理为S1.mat--S8.mat。优先按这八个明确文件读取，并
% 逐文件核验内部scenario_id及场景配置；若简洁文件不存在，再兼容原始长文件名。
%
% [统计修改-03 | 2026-09-24]
% 原因：按照论文原有tabularx模板输出表格，保持Ours/N-II/IN-II/N-III/
% M/D/MODDPG的列顺序，并以FFR替换Spacing。最优值依据表内显示精度自动加粗。

clc;

%% 1. 路径和基本设置
if ~exist('base_folder', 'var') || isempty(base_folder)
    script_dir = fileparts(mfilename('fullpath'));
    package_dir = fileparts(script_dir);
    base_folder = fullfile(package_dir, 'data', 'formal_results');
end
if ~exist('output_folder', 'var') || isempty(output_folder)
    output_folder = fullfile(package_dir, 'results', 'generated', ...
        'S1_S8_Metric_Summary');
end
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

evaluated_slot = 10;

% [Scenario ID, I, M, width, height]
scenario_catalog = [ ...
    1,  20,  2, 100, 100; ...
    2,  20,  2, 200, 200; ...
    3,  45,  3, 100, 100; ...
    4,  45,  3, 200, 200; ...
    5,  80,  4, 100, 100; ...
    6,  80,  4, 200, 200; ...
    7, 150,  8, 300, 300; ...
    8, 300, 15, 500, 500];

algorithm_fields = { ...
    'MyNSGA_II', 'NSGA_II', 'MOEA_D', ...
    'INSGA_II', 'NSGA_III', 'DRL_Baseline'};
algorithm_labels = { ...
    'HMD-NSGA-II', 'NSGA-II', 'MOEA/D', ...
    'INSGA-II', 'NSGA-III', 'MODDPG'};

%% 2. 读取S1.mat--S8.mat；必要时兼容原始长文件名
selected_paths = strings(8,1);
selected_dates = -inf(8,1);

for scenario_id = 1:8
    compact_path = fullfile(base_folder, ...
        sprintf('S%d_Formal.mat', scenario_id));
    if isfile(compact_path)
        selected_paths(scenario_id) = string(compact_path);
    end
end

% 兼容未重命名的原始结果文件；仅为缺失场景补位。
if any(strlength(selected_paths) == 0)
    matched_files = dir(fullfile(base_folder, '**', ...
        'Experiment_Statistical_Results_*.mat'));
    names = {matched_files.name};
    keep = ~contains(names, '_NormalizedMetrics', 'IgnoreCase', true) & ...
           ~contains(names, 'RAW_CHECKPOINT', 'IgnoreCase', true);
    matched_files = matched_files(keep);

    for file_idx = 1:numel(matched_files)
        file_name = matched_files(file_idx).name;
        tokens = regexp(file_name, ...
            '_I(\d+)_M(\d+)_R2(\d+)x(\d+)_Slot', 'tokens', 'once');
        if isempty(tokens)
            continue;
        end
        config = str2double(tokens);
        catalog_row = find( ...
            scenario_catalog(:,2) == config(1) & ...
            scenario_catalog(:,3) == config(2) & ...
            scenario_catalog(:,4) == config(3) & ...
            scenario_catalog(:,5) == config(4), 1);
        if isempty(catalog_row)
            continue;
        end

        scenario_id = scenario_catalog(catalog_row,1);
        if strlength(selected_paths(scenario_id)) == 0 && ...
                matched_files(file_idx).datenum > selected_dates(scenario_id)
            selected_dates(scenario_id) = matched_files(file_idx).datenum;
            selected_paths(scenario_id) = fullfile( ...
                matched_files(file_idx).folder, matched_files(file_idx).name);
        end
    end
end

missing_scenarios = find(strlength(selected_paths) == 0);
if ~isempty(missing_scenarios)
    error('Summary:MissingScenarios', ...
        '未找到以下场景的结果文件：%s', mat2str(missing_scenarios'));
end

fprintf('将使用以下结果文件：\n');
for scenario_id = 1:8
    fprintf('  S%d: %s\n', scenario_id, selected_paths(scenario_id));
end

%% 3. 逐场景、逐算法计算统计量
summary_rows = repmat(struct( ...
    'Scenario', '', 'Algorithm', '', ...
    'IGD_Mean', NaN, 'IGD_Std', NaN, 'IGD_ValidRuns', 0, ...
    'HV_Mean', NaN, 'HV_Std', NaN, 'HV_ValidRuns', 0, ...
    'FFR_Percent', NaN, 'SuccessfulRuns', 0, 'TotalRuns', 0, ...
    'Tmax_Mean', NaN, 'Tmax_Std', NaN, ...
    'Nnds_Mean', NaN, 'Nnds_Std', NaN, ...
    'Runtime_Mean', NaN, 'Runtime_Std', NaN, ...
    'SourceFile', ''), 8*numel(algorithm_fields), 1);

latex_values = cell(8, numel(algorithm_fields), 6);
metric_means = nan(8, numel(algorithm_fields), 6);
metric_stds = nan(8, numel(algorithm_fields), 6);
row_idx = 0;

fprintf('\n================ S1--S8 指标汇总（时隙 %d） ================\n', ...
    evaluated_slot);

for scenario_id = 1:8
    result_data = load(selected_paths(scenario_id), ...
        'all_scenario_results', 'nSlots', 'num_stat_runs', ...
        'experimental_scenarios', 'scenario_id');

    validate_result_identity(result_data, scenario_id, ...
        scenario_catalog(scenario_id,2:5), selected_paths(scenario_id));

    if evaluated_slot > result_data.nSlots
        error('Summary:InvalidSlot', ...
            'S%d仅包含%d个时隙，无法读取时隙%d。', ...
            scenario_id, result_data.nSlots, evaluated_slot);
    end

    % 当前正式文件通常只含一个场景，因此数据索引为1；兼容旧多场景文件。
    data_s_idx = resolve_scenario_index(result_data, scenario_catalog(scenario_id,2:5));

    fprintf('\nS%d (I=%d, M=%d, %dx%d m^2)\n', ...
        scenario_id, scenario_catalog(scenario_id,2), ...
        scenario_catalog(scenario_id,3), scenario_catalog(scenario_id,4), ...
        scenario_catalog(scenario_id,5));

    for algorithm_idx = 1:numel(algorithm_fields)
        row_idx = row_idx + 1;
        field_name = algorithm_fields{algorithm_idx};
        display_name = algorithm_labels{algorithm_idx};

        if ~isfield(result_data.all_scenario_results, field_name)
            error('Summary:MissingAlgorithm', ...
                'S%d结果中缺少算法字段%s。', scenario_id, field_name);
        end
        block = result_data.all_scenario_results.(field_name);

        igd = get_slot_vector(block.IGD{data_s_idx}, evaluated_slot);
        hv = get_slot_vector(block.HV{data_s_idx}, evaluated_slot);
        tmax = get_slot_vector(block.Tmax{data_s_idx}, evaluated_slot);
        nnds = get_slot_vector(block.NumFeasibleSolutions{data_s_idx}, evaluated_slot);
        runtime = get_slot_vector(block.Runtime{data_s_idx}, evaluated_slot);

        total_runs = numel(nnds);
        successful_mask = isfinite(nnds) & nnds > 0;
        successful_runs = sum(successful_mask);
        ffr = 100 * successful_runs / total_runs;

        [igd_mean, igd_std, igd_n] = finite_stats(igd);
        [hv_mean, hv_std, hv_n] = finite_stats(hv);
        [tmax_mean, tmax_std] = finite_stats(tmax);
        [nnds_mean, nnds_std] = finite_stats(nnds);
        [runtime_mean, runtime_std] = finite_stats(runtime);

        summary_rows(row_idx).Scenario = sprintf('S%d', scenario_id);
        summary_rows(row_idx).Algorithm = display_name;
        summary_rows(row_idx).IGD_Mean = igd_mean;
        summary_rows(row_idx).IGD_Std = igd_std;
        summary_rows(row_idx).IGD_ValidRuns = igd_n;
        summary_rows(row_idx).HV_Mean = hv_mean;
        summary_rows(row_idx).HV_Std = hv_std;
        summary_rows(row_idx).HV_ValidRuns = hv_n;
        summary_rows(row_idx).FFR_Percent = ffr;
        summary_rows(row_idx).SuccessfulRuns = successful_runs;
        summary_rows(row_idx).TotalRuns = total_runs;
        summary_rows(row_idx).Tmax_Mean = tmax_mean;
        summary_rows(row_idx).Tmax_Std = tmax_std;
        summary_rows(row_idx).Nnds_Mean = nnds_mean;
        summary_rows(row_idx).Nnds_Std = nnds_std;
        summary_rows(row_idx).Runtime_Mean = runtime_mean;
        summary_rows(row_idx).Runtime_Std = runtime_std;
        summary_rows(row_idx).SourceFile = char(selected_paths(scenario_id));

        latex_values{scenario_id,algorithm_idx,1} = format_pm(igd_mean,igd_std,4);
        latex_values{scenario_id,algorithm_idx,2} = format_pm(hv_mean,hv_std,4);
        latex_values{scenario_id,algorithm_idx,3} = sprintf('%.1f',ffr);
        latex_values{scenario_id,algorithm_idx,4} = format_pm(tmax_mean,tmax_std,4);
        latex_values{scenario_id,algorithm_idx,5} = format_pm(nnds_mean,nnds_std,1);
        latex_values{scenario_id,algorithm_idx,6} = format_pm(runtime_mean,runtime_std,3);
        metric_means(scenario_id,algorithm_idx,:) = [ ...
            igd_mean, hv_mean, ffr, tmax_mean, nnds_mean, runtime_mean];
        metric_stds(scenario_id,algorithm_idx,:) = [ ...
            igd_std, hv_std, NaN, tmax_std, nnds_std, runtime_std];

        fprintf(['  %-12s IGD=%-20s HV=%-20s FFR=%5.1f%% (%d/%d) ' ...
                 'Tmax=%-20s Nnds=%-16s Runtime=%s s\n'], ...
            display_name, latex_values{scenario_id,algorithm_idx,1}, ...
            latex_values{scenario_id,algorithm_idx,2}, ffr, ...
            successful_runs,total_runs, ...
            latex_values{scenario_id,algorithm_idx,4}, ...
            latex_values{scenario_id,algorithm_idx,5}, ...
            latex_values{scenario_id,algorithm_idx,6});
    end
end

%% 4. 保存CSV和MAT汇总
summary_table = struct2table(summary_rows);
csv_path = fullfile(output_folder, 'S1_S8_Metric_Summary.csv');
mat_path = fullfile(output_folder, 'S1_S8_Metric_Summary.mat');
writetable(summary_table, csv_path);
save(mat_path, 'summary_table', 'selected_paths', 'scenario_catalog', ...
    'algorithm_fields', 'algorithm_labels', 'evaluated_slot');

%% 5. 生成可直接复制的LaTeX表格
tex_path = fullfile(output_folder, 'S1_S8_LaTeX_Table.tex');
fid = fopen(tex_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('Summary:CannotWriteTex', '无法创建LaTeX文件：%s', tex_path);
end
fprintf(fid, '%% Automatically generated by Summarize_S1_S8_Metrics.m\n');
fprintf(fid, '\\begin{table*}[!t]\n');
fprintf(fid, '\t\\centering\n');
fprintf(fid, ['\t\\caption{Comprehensive Benchmark Results Across All ' ...
    'Experimental Scenarios (S1 to S8). Except for FFR, results are ' ...
    'presented as $\\text{Mean} \\pm \\text{Std}$ with condensed ' ...
    'precision for visual clarity; FFR is reported as a percentage. ' ...
    'The best results within each scenario and metric are bolded. ' ...
    '(N-II: NSGA-II, IN-II: INSGA-II, N-III: NSGA-III, ' ...
    'M/D: MOEA/D).}\n']);
fprintf(fid, '\t\\label{tab:comprehensive_global_benchmark}\n');
fprintf(fid, '\t\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid, '\t\\setlength{\\tabcolsep}{5pt}\n');
fprintf(fid, ['\t\\begin{tabularx}{\\textwidth}' ...
    '{@{} l l *{6}{>{\\centering\\arraybackslash}X} @{}}\n']);
fprintf(fid, '\t\t\\toprule\n');
fprintf(fid, ['\t\t\\textbf{Scenario} & \\textbf{Metric} & ' ...
    '\\textbf{Ours} & \\textbf{N-II} & \\textbf{IN-II} & ' ...
    '\\textbf{N-III} & \\textbf{M/D} & \\textbf{MODDPG} \\\\\n']);
fprintf(fid, '\t\t\\midrule\n');

% 以下内容通过fprintf的%s参数写入，因此每个LaTeX命令仅保留一个反斜杠。
metric_tex = {'$\mathrm{IGD} \downarrow$', ...
    '$\mathrm{HV} \uparrow$', ...
    '$\mathrm{FFR}\,(\%) \uparrow$', ...
    '$T_{\max}\,\mathrm{(s)} \downarrow$', ...
    '$N_{\mathrm{nds}} \uparrow$', ...
    '$\mathrm{Runtime}\,\mathrm{(s)} \downarrow$'};
% 与原论文模板一致：均值保留较高精度，标准差采用压缩精度以保留列间留白。
mean_digits = [4, 4, 1, 4, 1, 3];
std_digits = [2, 2, 0, 2, 1, 2];
larger_is_better = [false, true, true, false, true, false];
% 内部保存顺序为HMD, NSGA-II, MOEA/D, INSGA-II, NSGA-III, MODDPG；
% 论文模板顺序为Ours, N-II, IN-II, N-III, M/D, MODDPG。
latex_order = [1, 2, 4, 5, 3, 6];

for scenario_id = 1:8
    fprintf(fid, '\t\t%% --- SCENARIO S%d ---\n', scenario_id);
    for metric_idx = 1:numel(metric_tex)
        displayed_means = round( ...
            squeeze(metric_means(scenario_id,:,metric_idx)), ...
            mean_digits(metric_idx));
        finite_mask = isfinite(displayed_means);
        if larger_is_better(metric_idx)
            best_value = max(displayed_means(finite_mask));
        else
            best_value = min(displayed_means(finite_mask));
        end

        if metric_idx == 1
            fprintf(fid, '\t\t\\multirow{6}{*}{\\textbf{S%d}}', scenario_id);
        else
            fprintf(fid, '\t\t');
        end
        fprintf(fid, ' & %s', metric_tex{metric_idx});
        for column_idx = 1:numel(latex_order)
            algorithm_idx = latex_order(column_idx);
            mu = metric_means(scenario_id,algorithm_idx,metric_idx);
            sigma = metric_stds(scenario_id,algorithm_idx,metric_idx);
            is_best = isfinite(mu) && ...
                round(mu,mean_digits(metric_idx)) == best_value;
            cell_text = format_latex_cell(mu, sigma, ...
                mean_digits(metric_idx), std_digits(metric_idx), ...
                is_best, metric_idx == 3);
            fprintf(fid, ' & %s', cell_text);
        end
        fprintf(fid, ' \\\\\n');
    end
    if scenario_id < 8
        fprintf(fid, '\t\t\\cmidrule(lr){1-8}\n\n');
    end
end

fprintf(fid, '\t\t\\bottomrule\n');
fprintf(fid, '\t\\end{tabularx}\n');
fprintf(fid, '\\end{table*}\n\n');
fprintf(fid, ['%% Note: A dash (--) indicates that the corresponding metric ' ...
    'could not be computed because no valid front was obtained.\n']);
fclose(fid);

fprintf('\n汇总完成：\n');
fprintf('  CSV:   %s\n', csv_path);
fprintf('  MAT:   %s\n', mat_path);
fprintf('  LaTeX: %s\n', tex_path);

%% Local functions
function data_s_idx = resolve_scenario_index(result_data, target_config)
    data_s_idx = 1;
    if ~isfield(result_data, 'experimental_scenarios') || ...
            numel(result_data.experimental_scenarios) <= 1
        return;
    end
    for idx = 1:numel(result_data.experimental_scenarios)
        cfg = result_data.experimental_scenarios{idx};
        if iscell(cfg)
            cfg = cell2mat(cfg);
        end
        if numel(cfg) >= 4 && isequal(double(cfg(1:4)), double(target_config))
            data_s_idx = idx;
            return;
        end
    end
end

function validate_result_identity(result_data, expected_id, expected_config, file_path)
    if isfield(result_data, 'scenario_id') && ...
            ~isequal(double(result_data.scenario_id), double(expected_id))
        error('Summary:ScenarioIdMismatch', ...
            '%s内部scenario_id=%g，但文件被指定为S%d。', ...
            file_path, double(result_data.scenario_id), expected_id);
    end

    if isfield(result_data, 'current_scenario_config')
        actual_config = flatten_numeric_config( ...
            result_data.current_scenario_config);
        if numel(actual_config) < 4 || ...
                ~isequal(actual_config(1:4), double(expected_config))
            error('Summary:ScenarioConfigMismatch', ...
                '%s内部场景配置%s与S%d预期配置%s不一致。', ...
                file_path, mat2str(actual_config), expected_id, ...
                mat2str(expected_config));
        end
    end
end

function values = flatten_numeric_config(config)
    while iscell(config) && isscalar(config)
        config = config{1};
    end
    if iscell(config)
        values = cellfun(@double, config);
    else
        values = double(config);
    end
    values = values(:)';
end

function x = get_slot_vector(matrix, slot_idx)
    if size(matrix,2) < slot_idx
        error('Summary:MetricSlotMissing', ...
            '指标矩阵仅有%d个时隙，无法读取时隙%d。',size(matrix,2),slot_idx);
    end
    x = matrix(:,slot_idx);
end

function [mu,sigma,n] = finite_stats(x)
    x = x(isfinite(x));
    n = numel(x);
    if isempty(x)
        mu = NaN;
        sigma = NaN;
    else
        mu = mean(x);
        sigma = std(x);
    end
end

function text = format_pm(mu,sigma,digits)
    if ~isfinite(mu) || ~isfinite(sigma)
        text = '--';
        return;
    end
    % 避免sprintf将LaTeX命令\pm误解析为MATLAB转义序列。
    number_format = sprintf('%%.%df',digits);
    text = [sprintf(number_format,mu), ' $\pm$ ', ...
        sprintf(number_format,sigma)];
end

function text = format_latex_cell(mu,sigma,mean_digits,std_digits,is_best,is_ffr)
    if ~isfinite(mu)
        text = '$\mathrm{--}$';
        return;
    end
    mean_format = sprintf('%%.%df',mean_digits);
    mean_text = sprintf(mean_format,mu);
    if is_ffr
        if is_best
            mean_text = ['\mathbf{', mean_text, '}'];
        end
        text = ['$', mean_text, '$'];
    else
        std_format = sprintf('%%.%df',std_digits);
        std_text = sprintf(std_format,sigma);
        if is_best
            % 整个最优单元格加粗，使均值、正负号和标准差均清晰突出。
            text = ['$\mathbf{', mean_text, '}\,\boldsymbol{\pm}\,', ...
                '\mathbf{', std_text, '}$'];
        else
            text = ['$', mean_text, '\,\pm\,', std_text, '$'];
        end
    end
end
