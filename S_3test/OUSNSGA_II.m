% =========================================================================
% HMD-NSGA-II 算法主函数 (OUSNSGA_II.m) - 消融实验支持版
% 功能：支持记忆机制 (ISMM)、自适应策略 (AMS) 和混合编码 (Hybrid) 的独立开关
% =========================================================================
function Pop = OUSNSGA_II(problem, params, LastSlotArchive)
    %% 1. 初始化
    tic;
    % NSGA-II 粒子模板
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'CrowdingDistance', [], 'Rank', []); 
    
    % 定义部署区域边界
    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2);
    for k = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k-1)*2 + 1) = problem.area(1,1); % x_min
        problem.bounds.deployment(2, (k-1)*2 + 1) = problem.area(1,2); % x_max
        problem.bounds.deployment(1, (k-1)*2 + 2) = problem.area(2,1); % y_min
        problem.bounds.deployment(2, (k-1)*2 + 2) = problem.area(2,2); % y_max
    end

    % --- 新增：分层参数默认值（不影响原有消融实验） ---
    if ~isfield(params, 'pm_cont_coeff'), params.pm_cont_coeff = 0.8; end % 连续变异系数
    if ~isfield(params, 'pm_disc_coeff'), params.pm_disc_coeff = 1.2; end % 离散变异系数
    if ~isfield(params, 'mum_cont_coeff'), params.mum_cont_coeff = 1.1; end % 连续变异指数系数
    if ~isfield(params, 'mum_disc_coeff'), params.mum_disc_coeff = 0.9; end % 离散变异指数系数

    % --- 机制 1: ISMM (记忆机制) ---
    Pop = InitializePopulation_NSGAII(problem, params, particle_template, LastSlotArchive);
    
    % 初始排序：利用环境选择逻辑为初始种群计算 Rank 和 CrowdingDistance
    Pop = EnvironmentalSelection(Pop, params.N, particle_template);

    %% 2. 主循环
    for gen = 1:params.T_max
        % 1. 父代选择 (二元锦标赛)
        selected_indices = BinaryTournamentSelection(Pop, params.N);
        Parents = Pop(selected_indices);
        
        % 2. 生成子代 (Crossover + Mutation)
        Offspring = repmat(particle_template, params.N, 1);
        for i = 1:2:params.N-1
            p1 = Parents(i); p2 = Parents(i+1);
            
            % --- 机制 2: Hybrid Crossover (混合交叉算子) ---
            [c1_pos, c2_pos] = Crossover(p1.Position, p2.Position, problem, params);
            
            child1 = particle_template; child1.Position = c1_pos;
            child2 = particle_template; child2.Position = c2_pos;
            
            % --- 机制 3: AMS (自适应变异策略) ---
            child1 = Mutate(child1, problem, params, p1.CrowdingDistance, gen);
            child2 = Mutate(child2, problem, params, p2.CrowdingDistance, gen);
            
            Offspring(i) = child1;
            Offspring(i+1) = child2;
        end
        
        % 批量评估子代 (提速核心)
        OffResults = feval(problem.objFunc, Offspring, problem);
        for i = 1:params.N
            Offspring(i).Objectives = OffResults.Objectives(i, :);
            Offspring(i).Tmax = OffResults.Tmax(i);
        end

        % 3. 环境选择 (父代 + 子代合并)
        CombinedPop = [Pop; Offspring];
        % 过滤无效解 (Objectives 为空或包含 NaN/Inf 的解)
        valid_mask = ~cellfun(@(x) isempty(x) || any(isnan(x)) || any(isinf(x)), {CombinedPop.Objectives});
        CombinedPop = CombinedPop(valid_mask);
        
        if isempty(CombinedPop), error('所有解均无效，请检查评估函数或约束。'); end
        
        Pop = EnvironmentalSelection(CombinedPop, params.N, particle_template);
        
        if mod(gen, 100) == 0 || gen == params.T_max
            % 获取当前第一前沿的数量
            tmp_fronts = FindAllFronts(Pop);
            fprintf('迭代: %d/%d, 种群大小: %d, 第一前沿解数: %d\n', gen, params.T_max, numel(Pop), numel(tmp_fronts{1}));
        end
    end
    toc;
end

%% ========================== 核心机制辅助函数 ==========================

% --- 初始化函数 (ISMM) ---
function Pop = InitializePopulation_NSGAII(problem, params, template, LastSlotArchive)
    Pop = repmat(template, params.N, 1);
    nInherited = 0;
    
    if isfield(params, 'memory_ratio') && params.memory_ratio > 0 && ~isempty(LastSlotArchive)
        nInherited = min(floor(params.memory_ratio * params.N), numel(LastSlotArchive));
        idx = randperm(numel(LastSlotArchive), nInherited);
        for i = 1:nInherited
            Pop(i).Position = LastSlotArchive(idx(i)).Position;
        end
    end
    
    perturb = 0.05 * (problem.area(1,2) - problem.area(1,1));
    for i = (nInherited + 1):params.N
        if i == (nInherited + 1)
            dep = problem.initial_fog_deployment_flat;
        else
            dep = problem.initial_fog_deployment_flat + perturb .* (rand(1, problem.nFogNodes*2) - 0.5);
            dep = max(min(dep, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
        end
        Pop(i).Position.deployment = dep;
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
end

% --- 交叉函数 (Hybrid Switch) ---
function [c1_pos, c2_pos] = Crossover(p1_pos, p2_pos, problem, params)
    c1_pos = p1_pos; c2_pos = p2_pos;
    isHybrid = ~isfield(params, 'hybrid_enabled') || params.hybrid_enabled;
    
    if isHybrid
        % 混合模式：连续 SBX + 离散单点交叉
        if rand() <= params.pc
            [c1_pos.deployment, c2_pos.deployment] = SBX(p1_pos.deployment, p2_pos.deployment, params.mu);
            [c1_pos.bandwidth, c2_pos.bandwidth] = SBX(p1_pos.bandwidth, p2_pos.bandwidth, params.mu);
        end
        if rand() <= params.pc
            cp = randi(problem.nTerminals - 1);
            c1_pos.offloading = [p1_pos.offloading(1:cp), p2_pos.offloading(cp+1:end)];
            c2_pos.offloading = [p2_pos.offloading(1:cp), p1_pos.offloading(cp+1:end)];
        end
    else
        % 单一模式：全部视为连续变量处理 (消融对比)
        if rand() <= params.pc
            [c1_pos.deployment, c2_pos.deployment] = SBX(p1_pos.deployment, p2_pos.deployment, params.mu);
            [c1_pos.bandwidth, c2_pos.bandwidth] = SBX(p1_pos.bandwidth, p2_pos.bandwidth, params.mu);
            [c1_pos.offloading, c2_pos.offloading] = SBX(p1_pos.offloading, p2_pos.offloading, params.mu);
        end
    end
    c1_pos = BoundaryCheck(c1_pos, problem);
    c2_pos = BoundaryCheck(c2_pos, problem);
end

% --- 变异函数 (AMS Switch) ---
function child = Mutate(parent, problem, params, crowd_dist, gen)
    child = parent;
    % --- 修改：分层自适应参数计算 ---
    if isfield(params, 'adaptive_enabled') && params.adaptive_enabled
        xi_g = 1 - (gen / params.T_max)^2;
        phi_j = min(1, max(0.01, crowd_dist)); 
        if isinf(crowd_dist), phi_j = 0.5; end
        
        % 基础自适应参数（原有逻辑）
        base_pm = params.pm_min + (params.pm_max - params.pm_min) * xi_g * phi_j;
        base_mum = params.mum_max - (params.mum_max - params.mum_min) * xi_g * phi_j;
        
        % 分层参数：连续变量保守，离散变量活跃
        curr_pm_cont = base_pm * params.pm_cont_coeff; % 连续变异概率 = 基础 * 0.8
        curr_pm_disc = base_pm * params.pm_disc_coeff; % 离散变异概率 = 基础 * 1.2
        curr_mum_cont = base_mum * params.mum_cont_coeff; % 连续变异指数 = 基础 * 1.1（步长更小）
        curr_mum_disc = base_mum * params.mum_disc_coeff; % 离散变异指数（仅用于消融模式）
    else
        % 非自适应模式：同样分层
        curr_pm_cont = params.pm * params.pm_cont_coeff;
        curr_pm_disc = params.pm * params.pm_disc_coeff;
        curr_mum_cont = params.mum * params.mum_cont_coeff;
        curr_mum_disc = params.mum * params.mum_disc_coeff;
    end
    
    % --- 连续变量变异 (部署+带宽)：使用分层后的连续参数 ---
    if rand() <= curr_pm_cont
        for d = 1:numel(child.Position.deployment)
            child.Position.deployment(d) = PolynomialMutation(child.Position.deployment(d), ...
                problem.bounds.deployment(1,d), problem.bounds.deployment(2,d), curr_mum_cont);
        end
        for d = 1:numel(child.Position.bandwidth)
            child.Position.bandwidth(d) = PolynomialMutation(child.Position.bandwidth(d), ...
                problem.bounds.bandwidth(1,d), problem.bounds.bandwidth(2,d), curr_mum_cont);
        end
    end
    
    % --- 离散变量变异 (卸载)：使用分层后的离散参数 ---
    for d = 1:problem.nTerminals
        if rand() <= curr_pm_disc % 离散用更高的变异概率
            choices = 1:problem.nFogNodes;
            choices(choices == child.Position.offloading(d)) = [];
            if ~isempty(choices), child.Position.offloading(d) = choices(randi(numel(choices))); end
        end
    end
    
    % 边界检查
    child.Position.deployment = max(min(child.Position.deployment, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
    child.Position.bandwidth = max(min(child.Position.bandwidth, problem.bounds.bandwidth(2,:)), problem.bounds.bandwidth(1,:));
end

%% ========================== NSGA-II 核心逻辑函数 ==========================

function selected_idx = BinaryTournamentSelection(Pop, k)
    selected_idx = zeros(1, k);
    nPop = numel(Pop);
    for i = 1:k
        p1 = randi(nPop); p2 = randi(nPop);
        while p1 == p2, p2 = randi(nPop); end
        % 锦标赛规则：比 Rank，Rank 同则比 CrowdingDistance
        if Pop(p1).Rank < Pop(p2).Rank
            selected_idx(i) = p1;
        elseif Pop(p2).Rank < Pop(p1).Rank
            selected_idx(i) = p2;
        else
            if Pop(p1).CrowdingDistance > Pop(p2).CrowdingDistance
                selected_idx(i) = p1;
            else
                selected_idx(i) = p2;
            end
        end
    end
end

function Pop_next = EnvironmentalSelection(R, N, template)
    Fronts = FindAllFronts(R);
    Pop_next = repmat(template, 0, 1);
    for f = 1:numel(Fronts)
        curr_front = Fronts{f};
        if isempty(curr_front), continue; end
        
        % 标记 Rank
        for i = 1:numel(curr_front), curr_front(i).Rank = f; end
        
        % 计算拥挤度
        if numel(curr_front) > 1
            dists = CalculateCrowdingDistance(vertcat(curr_front.Objectives));
            for i = 1:numel(curr_front), curr_front(i).CrowdingDistance = dists(i); end
        else
            curr_front(1).CrowdingDistance = inf;
        end
        
        % 择优录取到下一代
        if (numel(Pop_next) + numel(curr_front)) <= N
            Pop_next = [Pop_next; curr_front];
        else
            [~, idx] = sort([curr_front.CrowdingDistance], 'descend');
            num_needed = N - numel(Pop_next);
            Pop_next = [Pop_next; curr_front(idx(1:num_needed))];
            break;
        end
    end
end

%% ========================== 基础算子 (SBX/Mutation/Check) ==========================

function [c1, c2] = SBX(p1, p2, mu)
    u = rand(size(p1));
    beta = zeros(size(p1));
    beta(u <= 0.5) = (2*u(u<=0.5)).^(1/(mu+1));
    beta(u > 0.5) = (1./(2*(1-u(u>0.5)))).^(1/(mu+1));
    c1 = 0.5 * ((1+beta).*p1 + (1-beta).*p2);
    c2 = 0.5 * ((1-beta).*p1 + (1+beta).*p2);
end

function val = PolynomialMutation(val, lb, ub, mum)
    u = rand();
    if u < 0.5
        delta = (2*u)^(1/(mum+1)) - 1;
    else
        delta = 1 - (2*(1-u))^(1/(mum+1));
    end
    val = val + delta * (ub - lb);
    val = max(min(val, ub), lb);
end

function pos = BoundaryCheck(pos, problem)
    pos.deployment = max(min(pos.deployment, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
    pos.bandwidth = max(min(pos.bandwidth, problem.bounds.bandwidth(2,:)), problem.bounds.bandwidth(1,:));
    pos.offloading = max(1, min(problem.nFogNodes, round(pos.offloading)));
end

% =========================================================================
% CalculateCrowdingDistance.m - 计算拥挤度距离
% =========================================================================
function crowd_dist = CalculateCrowdingDistance(objectives)
    nObj = size(objectives, 2); % 目标函数数量
    nPop = size(objectives, 1); % 解的数量
    
    if nPop <= 2 % 少于等于2个解，边界解拥挤距离为无穷大
        crowd_dist = inf(nPop, 1);
        return;
    end 


    crowd_dist = zeros(nPop, 1);
    
    % 对每个目标函数进行排序
    for j = 1:nObj
        [sorted_obj, sorted_indices] = sort(objectives(:,j));
        
        % 边界解的拥挤距离设为无穷大
        crowd_dist(sorted_indices(1)) = inf;
        crowd_dist(sorted_indices(end)) = inf;
        
        % 计算中间解的拥挤距离 (归一化范围)
        obj_range = sorted_obj(end) - sorted_obj(1);
        if obj_range == 0, obj_range = 1; end % 避免除以0，如果是平坦的，则不贡献距离
        
        for i = 2:nPop-1
            crowd_dist(sorted_indices(i)) = crowd_dist(sorted_indices(i)) + ...
                (sorted_obj(i+1) - sorted_obj(i-1)) / obj_range;
        end
    end
end
