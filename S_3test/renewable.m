% =========================================================================
% 指标修复与补全脚本 (v2)：增加 Spacing 计算
% =========================================================================
fprintf('正在基于原始解集重新计算所有指标 (含 Spacing)...\n');

for t = 1:nSlots
    % 1. 动态构建该时隙的 PF* 和 hv_ref
    objs_collector = [];
    for a_idx = 1:length(alg_names)
        for r = 1:num_stat_runs
            objs_collector = [objs_collector; current_run_archives{a_idx, t, r}];
        end
    end
    
    valid_mask = ~any(objs_collector >= 1e9 | isnan(objs_collector), 2);
    objs_collector = objs_collector(valid_mask, :);
    
    if ~isempty(objs_collector)
        pf_idx = FindNonDominatedSolutions(objs_collector);
        pf_star = objs_collector(pf_idx, :);
        % 参考点设为最大值的 1.1 倍，确保 HV 计算基准统一且为正
        hv_ref = max(objs_collector, [], 1) * 1.1; 
    else
        pf_star = [];
        hv_ref = [1.1, 1.1];
    end

    % 2. 遍历算法进行重算
    for a_idx = 1:length(alg_names)
        name = alg_names{a_idx};
        for r = 1:num_stat_runs
            objs = current_run_archives{a_idx, t, r};
            
            if isempty(objs)
                all_results.(name).IGD(r, t) = 1;
                all_results.(name).HV(r, t) = 0;
                all_results.(name).Spacing(r, t) = 1;
                continue;
            end
            
            % --- A. 计算 HV (调用专业函数) ---
            all_results.(name).HV(r, t) = calculateHV(objs, hv_ref);
            
            % --- B. 计算 IGD (归一化以收窄阴影) ---
            norm_objs = objs ./ hv_ref;
            norm_pf = pf_star ./ hv_ref;
            all_results.(name).IGD(r, t) = calculateIGD(norm_objs, norm_pf);
            
            % --- C. 计算 Spacing (衡量分布均匀性) ---
            % 算法逻辑：计算每个解到其他解的最小距离 di，再求 di 的标准差
            if size(objs, 1) > 1
                % 计算目标空间的点对点距离矩阵
                d = pdist2(objs, objs); 
                d(d==0) = inf; % 排除自身距离
                di = min(d, [], 2); % 找到每个点到其他点的最短距离
                % 计算 Spacing 值：均方根偏差
                all_results.(name).Spacing(r, t) = sqrt(sum((mean(di)-di).^2)/(size(objs,1)-1));
            else
                all_results.(name).Spacing(r, t) = 1; % 单点无意义，设为 1
            end
        end
    end
    fprintf('Slot %d 计算完成。\n', t);
end
fprintf('所有指标已修复！Spacing 数据已准备就绪。\n');