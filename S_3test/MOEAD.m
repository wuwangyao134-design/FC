% =========================================================================
% MOEA/D 算法主函数 (MOEAD.m)
% 描述: 按照 NSGA_II.m 的接口和结构进行封装，以便于在 main_compare.m 中直接调用。
% 作者: Gemini (根据用户代码定制) — 修正版
% 变更要点：
% 1) 修正邻域距离矩阵的欧氏距离实现
% 2) 健壮化 params.T、N=1 权重、初始化 z 的 NaN/Inf 处理
% 3) bandwidth 边界与尺寸断言
% 4) 可选 delta 概率全局父代抽样（默认 0.9）
% =========================================================================
function Pop = MOEAD(problem, params, ~)
    %% 1. 初始化 (批量化修改)
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', []); 
    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2);
    for k = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k-1)*2 + 1) = problem.area(1,1);
        problem.bounds.deployment(2, (k-1)*2 + 1) = problem.area(1,2);
        problem.bounds.deployment(1, (k-1)*2 + 2) = problem.area(2,1);
        problem.bounds.deployment(2, (k-1)*2 + 2) = problem.area(2,2);
    end
    
    Pop = repmat(particle_template, params.N, 1);
    initial_deployment_perturb_strength = 0.05 * (problem.area(2,2) - problem.area(2,1)); 
    
    % --- 批量初始化位置 ---
    for i = 1:params.N
        if i == 1
            deployment_pos = problem.initial_fog_deployment_flat;
        else
            deployment_pos = problem.initial_fog_deployment_flat + ...
                initial_deployment_perturb_strength .* (rand(1, problem.nFogNodes * 2) - 0.5);
            deployment_pos = max(deployment_pos, problem.bounds.deployment(1,:));
            deployment_pos = min(deployment_pos, problem.bounds.deployment(2,:));
        end
        bandwidth = problem.bounds.bandwidth(1,:) + rand(1, problem.nTerminals) .* (problem.bounds.bandwidth(2,:) - problem.bounds.bandwidth(1,:));
        offloading_plan = randi([1, problem.nFogNodes], 1, problem.nTerminals);
        
        Pop(i).Position.deployment = deployment_pos;
        Pop(i).Position.bandwidth  = bandwidth;
        Pop(i).Position.offloading = offloading_plan;
    end

    % === 【提速点 1：初始种群批量评估】 ===
    EvalResults = feval(problem.objFunc, Pop, problem);
    for i = 1:params.N
        Pop(i).Objectives = EvalResults.Objectives(i, :);
        Pop(i).Tmax = EvalResults.Tmax(i);
    end

    % --- 初始化权重、邻域与理想点 ---
    weights = init_weights(params.N, 2);
    D = permute(weights,[1 3 2]) - permute(weights,[3 1 2]);
    dist_matrix = sqrt(sum(D.^2, 3));
    [~, B] = sort(dist_matrix, 2);
    B = B(:, 1:params.T);
    z = min(vertcat(Pop.Objectives), [], 1);

    %% 2. 主循环
    for gen = 1:params.T_max
        % --- 步骤 A: 批量生成子代位置 ---
        Offspring = repmat(particle_template, params.N, 1);
        for i = 1:params.N
            if rand < params.delta, pool = B(i,:); else, pool = 1:params.N; end
            p_indices = pool(randperm(numel(pool), 2));
            parent1 = Pop(p_indices(1));
            parent2 = Pop(p_indices(2));
            
            if rand <= params.pc
                [child_pos, ~] = Crossover(parent1.Position, parent2.Position, problem, params);
            else
                child_pos = parent1.Position;
            end
            temp_child = particle_template;
            temp_child.Position = child_pos;
            Offspring(i) = Mutate(temp_child, problem, params);
        end

        % === 【提速点 2：子代批量评估】 ===
        OffResults = feval(problem.objFunc, Offspring, problem);
        for i = 1:params.N
            Offspring(i).Objectives = OffResults.Objectives(i, :);
            Offspring(i).Tmax = OffResults.Tmax(i);
        end

        % --- 步骤 B: 统一更新理想点与邻域 ---
        for i = 1:params.N
            child = Offspring(i);
            if ~all(isfinite(child.Objectives)), continue; end

            
            if all(child.Objectives < 1e8) % 仅当是可行解时更新理想点
                z = min(z, child.Objectives); % 更新理想点
            end
            
            % 更新邻域内个体的解
            for j = 1:params.T
                neighbor_idx = B(i, j);
                % 计算 Tchebycheff 聚合函数值
                g_old = max(abs(Pop(neighbor_idx).Objectives - z) .* weights(neighbor_idx, :));
                g_new = max(abs(child.Objectives - z) .* weights(neighbor_idx, :));
                if g_new < g_old
                    Pop(neighbor_idx) = child;
                end
            end
        end

        if mod(gen, 100) == 0
            fprintf('MOEA/D 迭代: %d/%d, 当前 z: [%.2f, %.2f]\n', gen, params.T_max, z(1), z(2));
        end
    end
    Pop = add_dummy_fields(Pop);
end

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
    
    % --- 协同核心：分层变异概率 ---
    pm_cont = params.pm;
    pm_disc = params.pm;
    if isfield(params, 'pm_cont_coeff'), pm_cont = params.pm * params.pm_cont_coeff; end
    if isfield(params, 'pm_disc_coeff'), pm_disc = params.pm * params.pm_disc_coeff; end
    
    % 1. 连续变量变异 (部署位置 + 带宽)
    if rand() <= pm_cont
        for d = 1:numel(child.Position.deployment)
            child.Position.deployment(d) = PolynomialMutation(child.Position.deployment(d), ...
                problem.bounds.deployment(1,d), problem.bounds.deployment(2,d), params.mum);
        end
        for d = 1:numel(child.Position.bandwidth)
            child.Position.bandwidth(d) = PolynomialMutation(child.Position.bandwidth(d), ...
                problem.bounds.bandwidth(1,d), problem.bounds.bandwidth(2,d), params.mum);
        end
    end
    
    % 2. 离散变量变异 (卸载决策)
    for d = 1:problem.nTerminals
        if rand() <= pm_disc
            choices = 1:problem.nFogNodes;
            choices(choices == child.Position.offloading(d)) = [];
            if ~isempty(choices)
                child.Position.offloading(d) = choices(randi(numel(choices)));
            end
        end
    end
    % 边界检查 (确保在变异后仍在范围内)
    child.Position.deployment = max(child.Position.deployment, problem.bounds.deployment(1,:));
    child.Position.deployment = min(child.Position.deployment, problem.bounds.deployment(2,:));
    child.Position.bandwidth = max(child.Position.bandwidth, problem.bounds.bandwidth(1,:));
    child.Position.bandwidth = min(child.Position.bandwidth, problem.bounds.bandwidth(2,:));
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

% =========================================================================
% 辅助函数
% =========================================================================

function W = init_weights(N, M)
    % 生成 M=2 的权重向量；对 N=1 做防护
    if M ~= 2
        error('Weight generation for M > 2 is not implemented in this MOEAD version.');
    end
    if N == 1
        W = [0.5 0.5];
        return;
    end
    W = zeros(N, 2);
    for i = 1:N
        W(i,1) = (i - 1) / (N - 1);
        W(i,2) = 1 - W(i,1);
    end
end

function Pop_out = add_dummy_fields(Pop_in)
    % 为输出的种群添加 NSGA-II 特有的字段，以确保兼容性
    particle_template_final = struct('Position', [], 'Objectives', [], 'Tmax', [], 'CrowdingDistance', [], 'Rank', []);
    Pop_out = repmat(particle_template_final, numel(Pop_in), 1);
    for i = 1:numel(Pop_in)
        Pop_out(i).Position = Pop_in(i).Position;
        Pop_out(i).Objectives = Pop_in(i).Objectives;
        Pop_out(i).Tmax = Pop_in(i).Tmax;
        Pop_out(i).Rank = NaN;               % MOEA/D 不使用 Rank
        Pop_out(i).CrowdingDistance = NaN;   % MOEA/D 不使用 CrowdingDistance
    end
end
