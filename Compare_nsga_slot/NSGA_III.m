% =========================================================================
% NSGA-III 算法主函数 (ONSGA_III.m) - 统一封装版
% =========================================================================
function Pop = NSGA_III(problem, params, ~)
    %% 1. 初始化
    tic;
    % 统一粒子模板
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'Rank', [], 'Association', []); 
    
    % 初始化部署边界
    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2);
    for k = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k-1)*2 + 1) = problem.area(1,1);
        problem.bounds.deployment(2, (k-1)*2 + 1) = problem.area(1,2);
        problem.bounds.deployment(1, (k-1)*2 + 2) = problem.area(2,1);
        problem.bounds.deployment(2, (k-1)*2 + 2) = problem.area(2,2);
    end

    % --- NSGA-III 核心：生成参考点 ---
    % M 为目标数，p 为划分数
    M = 2; % 假设为 Delay 和 Energy 两个目标
    Z = GenerateReferencePoints(M, params.p); 
    params.N = size(Z, 1); % 种群大小自动对齐参考点数

    % 初始化种群位置 (复用 HMD 的逻辑)
    Pop = repmat(particle_template, params.N, 1);
    perturb = 0.05 * (problem.area(1,2) - problem.area(1,1));
    for i = 1:params.N
        if i == 1
            dep = problem.initial_fog_deployment_flat;
        else
            dep = problem.initial_fog_deployment_flat + perturb .* (rand(1, problem.nFogNodes*2) - 0.5);
        end
        Pop(i).Position.deployment = max(min(dep, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
        Pop(i).Position.bandwidth = problem.bounds.bandwidth(1,:) + rand(1, problem.nTerminals) .* ...
                                    (problem.bounds.bandwidth(2,:) - problem.bounds.bandwidth(1,:));
        Pop(i).Position.offloading = randi([1, problem.nFogNodes], 1, problem.nTerminals);
    end

    % 初始批量评估
    EvalRes = feval(problem.objFunc, Pop, problem);
    for i = 1:params.N
        Pop(i).Objectives = EvalRes.Objectives(i,:);
        Pop(i).Tmax = EvalRes.Tmax(i);
    end

    %% 2. 主循环
    for gen = 1:params.T_max
        % 1. 父代选择 (NSGA-III 通常采用简单的随机选择或基于 Rank 的选择)
        selected_indices = randi(params.N, 1, params.N);
        Parents = Pop(selected_indices);
        
        % 2. 生成子代 (复用 Crossover & Mutate)
        Offspring = repmat(particle_template, params.N, 1);
        for i = 1:2:params.N-1
            [c1_pos, c2_pos] = Crossover(Parents(i).Position, Parents(i+1).Position, problem, params);
            child1 = particle_template; child1.Position = c1_pos;
            child2 = particle_template; child2.Position = c2_pos;
            % 这里可以按需加入自适应变异逻辑
            child1 = Mutate(child1, problem, params);
            child2 = Mutate(child2, problem, params);
            Offspring(i) = child1; Offspring(i+1) = child2;
        end
        
        % 批量评估子代
        OffResults = feval(problem.objFunc, Offspring, problem);
        for i = 1:params.N
            Offspring(i).Objectives = OffResults.Objectives(i,:);
            Offspring(i).Tmax = OffResults.Tmax(i);
        end
        
        % 3. 环境选择 (合并 -> 非支配排序 -> 参考点小生境选择)
        CombinedPop = [Pop; Offspring];
        Pop = EnvironmentalSelection_NSGA3(CombinedPop, params.N, Z);
        
        if mod(gen, 100) == 0
            fprintf('NSGA-III 迭代: %d/%d, 第一层解数: %d\n', gen, params.T_max, sum([Pop.Rank]==1));
        end
    end
    toc;
end

%% ========================== NSGA-III 核心函数 ==========================
function Pop_next = EnvironmentalSelection_NSGA3(R, N, Z)
    Fronts = FindAllFronts(R);
    Pop_next = []; % 初始化为空结构体数组
    last_front_idx = 0;
    
    for f = 1:numel(Fronts)
        current_front = Fronts{f};
        if isempty(current_front), continue; end
        
        % 标记 Rank
        for i = 1:numel(current_front)
            current_front(i).Rank = f;
        end
        
        % 检查当前 Pop_next 的规模
        num_in_next = numel(Pop_next);
        if num_in_next + numel(current_front) <= N
            Pop_next = [Pop_next; current_front];
        else
            last_front_idx = f;
            break;
        end
    end
    
    % 如果已经刚好凑够 N 个解，直接返回
    if numel(Pop_next) == N
        return;
    end
    
    % --- 关键修复：如果 Pop_next 依然是空的，需要初始化为对应结构的空值 ---
    if isempty(Pop_next)
        % 创建一个带字段但长度为0的结构体，防止 AssociateAndSelect 报错
        Pop_next = R(1:0); 
    end
    
    % 处理临界层
    if last_front_idx > 0
        LastFront = Fronts{last_front_idx};
        num_needed = N - numel(Pop_next);
        % 调用关联选择
        SelectedFromLast = AssociateAndSelect(Pop_next, LastFront, Z, num_needed);
        Pop_next = [Pop_next; SelectedFromLast];
    end
end

%% ========================== NSGA-III 高级核心算子 ==========================
function SelectedFromLast = AssociateAndSelect(Pop_next, LastFront, Z, num_needed)
    % 1. 鲁棒性提取目标值
    PopNextObjs = vertcat(Pop_next.Objectives);
    LastFrontObjs = vertcat(LastFront.Objectives);
    
    if isempty(LastFrontObjs)
        SelectedFromLast = [];
        return;
    end
    
    AllObjs = [PopNextObjs; LastFrontObjs];
    nNext = size(PopNextObjs, 1);
    nLast = size(LastFrontObjs, 1);
    nObj = size(AllObjs, 2);

    % 2. 目标归一化 (Normalization)
    ideal_point = min(AllObjs, [], 1);
    f_prime = AllObjs - ideal_point;
    
    % 找到极值点
    extreme_points = zeros(nObj, nObj);
    for j = 1:nObj
        [~, idx] = min(max(f_prime ./ (ones(size(f_prime, 1), 1) * (1e-6 + (1:nObj == j))), [], 2));
        extreme_points(j, :) = f_prime(idx, :);
    end
    
    % --- 修正点 1: 计算拦截点 a 并增加回退机制 ---
    % 拦截点计算的鲁棒化处理
    hyp = pinv(extreme_points) * ones(nObj, 1);
    a = 1 ./ (hyp + 1e-10);
    
    % --- 协同保障：防止惩罚项破坏参考点关联 ---
    % 如果所有解都违规，或者拦截点计算异常，退退到 NADIR 点
    if any(a < 1e-6) || any(isnan(a)) || any(isinf(a))
        a = max(f_prime, [], 1)'; 
    end
    
    f_double_prime = f_prime ./ (a' + 1e-10);

    % --- 修正点 2: 参考点单位化 ---
    % 确保每个参考点向量长度为 1，使垂直距离计算在同一刻度下
    Z_unit = Z ./ (sqrt(sum(Z.^2, 2)) + 1e-10);
    nRef = size(Z_unit, 1);

    % 3. 关联参考点 (Association)
    nPop = size(f_double_prime, 1);
    all_dist = zeros(nPop, nRef);
    all_rho = zeros(nPop, 1);
    
    for i = 1:nPop
        w = Z_unit; % 使用单位化后的参考点
        w_norm = sum(w.^2, 2); % 此时 w_norm 理论上应接近 1
        inner_product = f_double_prime(i, :) * w';
        proj = (inner_product' ./ (w_norm + 1e-10)) .* w;
        dist_to_rays = sqrt(sum((f_double_prime(i, :) - proj).^2, 2));
        
        all_dist(i, :) = dist_to_rays';
        [~, min_idx] = min(dist_to_rays);
        all_rho(i) = min_idx;
    end
    
    rho_next = all_rho(1:nNext);
    rho_last = all_rho(nNext+1:end);
    dist_last = all_dist(nNext+1:end, :);
    
    % 4. 小生境保留策略 (Niche Preservation)
    niche_count = zeros(1, nRef);
    for j = 1:nRef
        niche_count(j) = sum(rho_next == j);
    end
    
    SelectedFromLast = [];
    is_available = true(nLast, 1);
    
    while numel(SelectedFromLast) < num_needed
        min_niche = min(niche_count);
        relevant_Z_indices = find(niche_count == min_niche);
        j_bar = relevant_Z_indices(randi(numel(relevant_Z_indices)));
        
        candidates = find((rho_last(:) == j_bar) & is_available);
        
        if ~isempty(candidates)
            if niche_count(j_bar) == 0
                [~, min_idx] = min(dist_last(candidates, j_bar));
                best_cand_idx = candidates(min_idx);
            else
                best_cand_idx = candidates(randi(numel(candidates)));
            end
            
            SelectedFromLast = [SelectedFromLast; LastFront(best_cand_idx)];
            is_available(best_cand_idx) = false;
            niche_count(j_bar) = niche_count(j_bar) + 1;
        else
            niche_count(j_bar) = inf;
        end
    end
end





function Z = GenerateReferencePoints(M, p)
    % 生成均匀分布的参考点 (Das & Dennis 方案简易版)
    Z = combinations(M, p) / p;
end




function C = combinations(M, p)
    if M == 1, C = p; return; end
    C = [];
    for i = 0:p
        C = [C; [i*ones(size(combinations(M-1, p-i),1),1), combinations(M-1, p-i)]];
    end
end

% 注：AssociateAndSelect 包含目标归一化、垂直距离计算和小生境计数逻辑
% 由于代码较长，建议将之前 DNSGA-II 的 Crossover/Mutate 放入同级目录下调用

% =========================================================================
% 辅助函数 - 交叉 (Crossover.m)
% 模拟二进制交叉 (SBX) for continuous, 均匀交叉 for discrete
% =========================================================================
function [child1_pos, child2_pos] = Crossover(parent1_pos, parent2_pos, problem, params)
    child1_pos = parent1_pos; % 默认复制父代
    child2_pos = parent2_pos;


    % 1. 连续部分交叉 (Simulated Binary Crossover - SBX)
    % 部署位置
    if rand() <= params.pc % 根据交叉概率决定是否交叉
        [child1_dep, child2_dep] = SBX(parent1_pos.deployment, parent2_pos.deployment, params.mu);
        child1_pos.deployment = child1_dep;
        child2_pos.deployment = child2_dep;
        [child1_bw, child2_bw] = SBX(parent1_pos.bandwidth, parent2_pos.bandwidth, params.mu);
        child1_pos.bandwidth = child1_bw;
        child2_pos.bandwidth = child2_bw;
    end


    % 2. 离散部分交叉 (Uniform Crossover for discrete part)
    if rand() <= params.pc
        for d = 1:problem.nTerminals
            if rand() < 0.5 % 50%概率交换
                temp = child1_pos.offloading(d);
                child1_pos.offloading(d) = child2_pos.offloading(d);
                child2_pos.offloading(d) = temp;
            end
        end
    end
    
    % 边界检查 (在变异后统一进行更合理，但为了安全在交叉后也进行一次)
    child1_pos.deployment = max(child1_pos.deployment, problem.bounds.deployment(1,:));
    child1_pos.deployment = min(child1_pos.deployment, problem.bounds.deployment(2,:));
    child1_pos.bandwidth = max(child1_pos.bandwidth, problem.bounds.bandwidth(1,:));
    child1_pos.bandwidth = min(child1_pos.bandwidth, problem.bounds.bandwidth(2,:));


    child2_pos.deployment = max(child2_pos.deployment, problem.bounds.deployment(1,:));
    child2_pos.deployment = min(child2_pos.deployment, problem.bounds.deployment(2,:));
    child2_pos.bandwidth = max(child2_pos.bandwidth, problem.bounds.bandwidth(1,:));
    child2_pos.bandwidth = min(child2_pos.bandwidth, problem.bounds.bandwidth(2,:));
end


% 辅助函数：Simulated Binary Crossover (SBX)
function [child1_val, child2_val] = SBX(parent1_val, parent2_val, mu)
    r = rand(size(parent1_val));
    beta = zeros(size(parent1_val));
    beta(r <= 0.5) = (2 * r(r <= 0.5)).^(1/(mu + 1));
    beta(r > 0.5) = (1 ./ (2 * (1 - r(r > 0.5)))).^(1/(mu + 1));
    
    child1_val = 0.5 * (((1 + beta) .* parent1_val) + ((1 - beta) .* parent2_val));
    child2_val = 0.5 * (((1 - beta) .* parent1_val) + ((1 + beta) .* parent2_val));
end


% =========================================================================
% 辅助函数 - 变异 (Mutate.m)
% 多项式变异 (Polynomial Mutation) for continuous, 随机变异 for discrete
% =========================================================================
function child = Mutate(parent, problem, params)
    child = parent;
    
    % --- 协同核心：分层参数 ---
    pm_cont = params.pm;
    pm_disc = params.pm;
    if isfield(params, 'pm_cont_coeff'), pm_cont = params.pm * params.pm_cont_coeff; end
    if isfield(params, 'pm_disc_coeff'), pm_disc = params.pm * params.pm_disc_coeff; end
    
    % 1. 连续变量协同变异 (部署位置 + 带宽)
    % 连续变量通常需要较低的变异频率以维持稳定性
    if rand() <= pm_cont
        % 部署位置变异 (多项式变异)
        for d = 1:numel(child.Position.deployment)
            child.Position.deployment(d) = PolynomialMutation(child.Position.deployment(d), ...
                problem.bounds.deployment(1,d), problem.bounds.deployment(2,d), params.mum);
        end
        % 带宽分配变异
        for d = 1:numel(child.Position.bandwidth)
            child.Position.bandwidth(d) = PolynomialMutation(child.Position.bandwidth(d), ...
                problem.bounds.bandwidth(1,d), problem.bounds.bandwidth(2,d), params.mum);
        end
    end
    
    % 2. 离散变量协同变异 (卸载决策)
    % 离散变量通常赋予更高的变异频率 (pm_disc_coeff > 1)，以增强探索能力
    for d = 1:problem.nTerminals
        if rand() <= pm_disc
            choices = 1:problem.nFogNodes;
            choices(choices == child.Position.offloading(d)) = [];
            if ~isempty(choices)
                child.Position.offloading(d) = choices(randi(numel(choices)));
            end
        end
    end
    
    % 边界检查 (确保协同后的变量依然在物理界限内)
    child.Position.deployment = max(min(child.Position.deployment, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
    child.Position.bandwidth = max(min(child.Position.bandwidth, problem.bounds.bandwidth(2,:)), problem.bounds.bandwidth(1,:));
end



% 辅助函数：Polynomial Mutation
function val = PolynomialMutation(val, lower_bound, upper_bound, mum)
    delta = rand();
    if delta < 0.5
        delta_q = (2 * delta)^(1/(mum + 1)) - 1;
    else
        delta_q = 1 - (2 * (1 - delta))^(1/(mum + 1));
    end
    val = val + delta_q * (upper_bound - lower_bound);
    
    % 边界检查
    val = max(val, lower_bound);
    val = min(val, upper_bound);
end
