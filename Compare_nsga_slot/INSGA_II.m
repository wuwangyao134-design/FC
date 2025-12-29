function Pop = INSGA_II(problem, params, ~)
    %% 1. 初始化
    tic;
    % 统一粒子模板，确保字段名与 NSGA-II 逻辑严谨对齐
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'CrowdingDistance', [], 'Rank', []); 
    N = params.N; 
    MaxIt = params.T_max;
    
    % ======================== 核心修复代码段开始 ========================
    % 动态初始化部署边界 (确保维度为 2 x D)
    if ~isfield(problem.bounds, 'deployment') || size(problem.bounds.deployment, 2) ~= problem.nFogNodes * 2
        % 强制构造: 第一行下界, 第二行上界 [cite: 67]
        x_min = problem.area(1,1); x_max = problem.area(1,2);
        y_min = problem.area(2,1); y_max = problem.area(2,2);
        
        lower_bound_row = repmat([x_min, y_min], 1, problem.nFogNodes);
        upper_bound_row = repmat([x_max, y_max], 1, problem.nFogNodes);
        
        problem.bounds.deployment = [lower_bound_row; upper_bound_row];
    end
    
    Pop = repmat(particle_template, N, 1);
    perturb_strength = 0.05 * (problem.area(1,2) - problem.area(1,1)); 
    
    % 提取边界行向量，确保后面运算兼容 [cite: 180]
    LB_dep = problem.bounds.deployment(1, :);
    UB_dep = problem.bounds.deployment(2, :);
    
    for i = 1:N
        if i == 1
             dep_pos = problem.initial_fog_deployment_flat(:)'; % 强制行向量
        else
             % 生成扰动并确保为行向量
             rand_noise = (rand(1, problem.nFogNodes * 2) - 0.5);
             dep_pos = problem.initial_fog_deployment_flat(:)' + perturb_strength .* rand_noise;
        end
        
        % 修正后的边界检查：确保所有输入均为行向量 [cite: 67, 70]
        Pop(i).Position.deployment = max(min(dep_pos, UB_dep), LB_dep);
        
        % 其余初始化保持不变...
        Pop(i).Position.bandwidth = problem.bounds.bandwidth(1,:) + rand(1, problem.nTerminals) .* ...
                                    (problem.bounds.bandwidth(2,:) - problem.bounds.bandwidth(1,:));
        Pop(i).Position.offloading = randi([1, problem.nFogNodes], 1, problem.nTerminals);
        
        EvalRes = feval(problem.objFunc, Pop(i).Position, problem);
        Pop(i).Objectives = EvalRes.Objectives;
        Pop(i).Tmax = EvalRes.Tmax;
    end

    %% 2. 主循环
    for t = 1:MaxIt
        % 关键改进：自适应权重因子 C (从 1 降至 0.02) [cite: 146]
        C = 1 - 0.98 * (t / MaxIt); 
        
        Offspring = repmat(particle_template, N, 1);
        for i = 1:N

                pos = CombinePos(Pop(i).Position);
                % 计算分层步长掩码 (Mask)
                % 连续部分权重设为 1.0, 离散部分权重设为 params.pm_disc_coeff
                len_cont = problem.nFogNodes * 2 + problem.nTerminals;
                len_disc = problem.nTerminals;
                coeff_mask = [ones(1, len_cont), ones(1, len_disc) * (params.pm_disc_coeff/params.pm_cont_coeff)];
        
                if rand < C
                    L = Levy(numel(pos));
                    % 将系数掩码应用到步长上，实现协同探索
                    new_pos = pos + C * L .* pos .* coeff_mask; 
                else
                    G = randn(1, numel(pos)); 
                    new_pos = pos + (1 - C) * G .* pos .* coeff_mask;
                end
                Offspring(i).Position = SeparatePos(new_pos, problem);
        end
        
        % 批量评估子代
        OffResults = feval(problem.objFunc, Offspring, problem);
        for i = 1:N
            Offspring(i).Objectives = OffResults.Objectives(i,:);
            Offspring(i).Tmax = OffResults.Tmax(i);
        end
        
        % 环境选择：合并 -> 非支配排序 -> 拥挤度距离 [cite: 45, 90, 177]
        Combined = [Pop; Offspring];
        Fronts = FindAllFronts(Combined); 
        Pop = [];
        for f = 1:numel(Fronts)
            current_front = Fronts{f};
            if isempty(current_front), continue; end
            
            % 关键修正：确保在选择前重新计算该层拥挤度距离 [cite: 83, 178]
            current_front = CalculateCrowdingDistanceStruct(current_front);
            
            if numel(Pop) + numel(current_front) <= N
                for i_f = 1:numel(current_front), current_front(i_f).Rank = f; end
                Pop = [Pop; current_front];
            else
                num_needed = N - numel(Pop);
                % 按拥挤度距离降序排列，保留分布更广的解 [cite: 45, 85]
                [~, idx] = sort([current_front.CrowdingDistance], 'descend');
                temp_front = current_front(idx(1:num_needed));
                for i_f = 1:numel(temp_front), temp_front(i_f).Rank = f; end
                Pop = [Pop; temp_front];
                break;
            end
        end
        
        if mod(t, 100) == 0
            fprintf('INSGA-II 迭代: %d/%d, 第一层解数: %d\n', t, MaxIt, sum([Pop.Rank]==1));
        end
    end
    toc;
end

function front = CalculateCrowdingDistanceStruct(front)
    % 修正版：直接操作结构体数组，计算拥挤度距离 [cite: 83]
    nPop = numel(front);
    if nPop == 0, return; end
    for i=1:nPop, front(i).CrowdingDistance = 0; end
    
    objs = vertcat(front.Objectives);
    nObj = size(objs, 2);
    
    for j = 1:nObj
        [sorted_vals, sorted_idx] = sort(objs(:,j));
        front(sorted_idx(1)).CrowdingDistance = inf; % 边界处理 [cite: 85]
        front(sorted_idx(end)).CrowdingDistance = inf;
        
        range = sorted_vals(end) - sorted_vals(1);
        if range == 0, continue; end
        
        for i = 2:nPop-1
            front(sorted_idx(i)).CrowdingDistance = front(sorted_idx(i)).CrowdingDistance + ...
                (sorted_vals(i+1) - sorted_vals(i-1)) / range;
        end
    end
end

function pos_vector = CombinePos(PositionStruct)
    % 向量化封装：处理混合编码变量 [cite: 132]
    pos_vector = [PositionStruct.deployment(:)', ...
                  PositionStruct.bandwidth(:)', ...
                  PositionStruct.offloading(:)'];
end

function NewPositionStruct = SeparatePos(pos_vector, problem)
    % 1. 物理长度切片
    len_dep = problem.nFogNodes * 2;
    len_bw = problem.nTerminals;
    
    dep = pos_vector(1 : len_dep);
    bw = pos_vector(len_dep + 1 : len_dep + len_bw);
    off = pos_vector(len_dep + len_bw + 1 : end);
    
    % 2. 连续变量：维持高精度裁剪
    NewPositionStruct.deployment = max(min(dep, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
    NewPositionStruct.bandwidth = max(min(bw, problem.bounds.bandwidth(2,:)), problem.bounds.bandwidth(1,:));
    
    % 3. 离散变量：强化协同约束
    % 直接取整可能导致大量重复，加入一个微小的离散扰动防止搜索停滞
    off_rounded = round(off);
    % 确保索引不越界，且为整数
    off_rounded = max(1, min(problem.nFogNodes, off_rounded));
    NewPositionStruct.offloading = off_rounded;
end

function L = Levy(dim)
    % 莱维分布步长生成 (文献式 5-7) [cite: 117, 118]
    gamma_val = 1.5;
    sigma_u = (gamma(1+gamma_val) * sin(pi*gamma_val/2) / ...
              (gamma((1+gamma_val)/2) * gamma_val * 2^((gamma_val-1)/2)))^(1/gamma_val);
    u = randn(1, dim) * sigma_u;
    v = randn(1, dim);
    L = u ./ abs(v).^(1/gamma_val);
end

