function metrics = CalculateMetricsOnly(archive, igd_reference_front, hv_reference_point)
    metrics = struct('IGD', NaN, 'HV', NaN, 'Spacing', NaN, 'Spread', NaN, 'Runtime', NaN, 'NumFeasibleSolutions', 0);
    
    if isempty(archive)
        return;
    end
    
    % RunTime 字段应在主循环中附加到 archive struct 的第一个元素上
    if isfield(archive, 'RunTime') && ~isempty(archive) && ~isempty(archive(1).RunTime)
        metrics.Runtime = archive(1).RunTime;
    end
    
    objectives_final = vertcat(archive.Objectives);
    % 过滤惩罚值和无效解，确保只对可行解计算指标
    valid_obj_idx = ~any(objectives_final >= 1e9 | isnan(objectives_final) | isinf(objectives_final), 2);
    filtered_objects = objectives_final(valid_obj_idx, :); 
    
    metrics.NumFeasibleSolutions = size(filtered_objects, 1);

    if ~isempty(filtered_objects)
        % HV 计算
        if ~isempty(hv_reference_point) && ~any(isnan(hv_reference_point)) && ~any(isinf(hv_reference_point))
            metrics.HV = calculateHV(filtered_objects, hv_reference_point);
        end
        % IGD 计算
        if ~isempty(igd_reference_front) && size(igd_reference_front,1) >= 2
            metrics.IGD = calculateIGD(filtered_objects, igd_reference_front);
        end
        % Spread计算 - 使用新的calculateSpread函数
        if ~isempty(igd_reference_front) && size(filtered_objects,1) >= 2
            metrics.Spread = calculateSpread(filtered_objects, igd_reference_front);
        end
    end

    if size(filtered_objects, 1) >= 2
        % Spacing 计算
        N_PF = size(filtered_objects, 1);
        distances = zeros(N_PF, 1);
        for i_pf = 1:N_PF
            current_point = filtered_objects(i_pf,:);
            min_dist = inf;
            for j_pf = 1:N_PF
                if i_pf == j_pf, continue; end
                dist = norm(current_point - filtered_objects(j_pf,:));
                if dist < min_dist, min_dist = dist; end
            end
            distances(i_pf) = min_dist;
        end
        d_bar = mean(distances);
        metrics.Spacing = sqrt(sum((distances - d_bar).^2) / (N_PF - 1));
    end
end
