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
    %% 0. 参数与断言
    if ~isfield(params,'N') || params.N<=0
        error('params.N 必须为正整数');
    end
    if ~isfield(params,'T_max') || params.T_max<=0
        error('params.T_max 必须为正整数');
    end
    if ~isfield(params,'T'), params.T = min(20, params.N); end          % 默认邻域大小
    params.T = min(max(1, params.T), params.N);                          % 夹取到 [1,N]
    if ~isfield(params,'pc'), params.pc = 0.9; end                       % 默认交叉概率
    if ~isfield(params,'delta'), params.delta = 0.9; end                 % 邻域采样概率（其余为全局）
    if ~isfield(problem,'nObj') || problem.nObj~=2
        error('当前版本仅支持 nObj=2 的 MOEA/D。');
    end
    if ~isfield(problem,'nFogNodes') || ~isfield(problem,'nTerminals') || ~isfield(problem,'area')
        error('problem 缺少 nFogNodes / nTerminals / area 定义。');
    end
    % bandwidth 边界断言（需外部构造）
    assert(isfield(problem,'bounds') && isfield(problem.bounds,'bandwidth') ...
        && all(size(problem.bounds.bandwidth)==[2 problem.nTerminals]), ...
        'problem.bounds.bandwidth 缺失或尺寸不匹配，应为 2 x nTerminals。');

    %% 1. 初始化
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', []); 

    % --- 部署变量的边界（复用 NSGA-II 初始化逻辑） ---
    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2);
    for k = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k-1)*2 + 1) = problem.area(1,1); % x_min
        problem.bounds.deployment(2, (k-1)*2 + 1) = problem.area(1,2); % x_max
        problem.bounds.deployment(1, (k-1)*2 + 2) = problem.area(2,1); % y_min
        problem.bounds.deployment(2, (k-1)*2 + 2) = problem.area(2,2); % y_max
    end

    Pop = repmat(particle_template, params.N, 1);
    initial_deployment_perturb_strength = 0.05 * (problem.area(2,2) - problem.area(2,1)); 

    for i = 1:params.N
        if i == 1
            deployment_pos = problem.initial_fog_deployment_flat;
        else
            deployment_pos = problem.initial_fog_deployment_flat + ...
                initial_deployment_perturb_strength .* (rand(1, problem.nFogNodes * 2) - 0.5);
            deployment_pos = max(deployment_pos, problem.bounds.deployment(1,:));
            deployment_pos = min(deployment_pos, problem.bounds.deployment(2,:));
        end

        bw_lb = problem.bounds.bandwidth(1,:);
        bw_ub = problem.bounds.bandwidth(2,:);
        bandwidth = bw_lb + rand(1, problem.nTerminals) .* (bw_ub - bw_lb);

        offloading_plan = randi([1, problem.nFogNodes], 1, problem.nTerminals);

        Pop(i).Position.deployment = deployment_pos;
        Pop(i).Position.bandwidth  = bandwidth;
        Pop(i).Position.offloading = offloading_plan;

        EvalResults = feval(problem.objFunc, Pop(i).Position, problem);
        Pop(i).Objectives = EvalResults.Objectives;
        Pop(i).Tmax = EvalResults.Tmax;
    end

    % --- MOEA/D 特定部分的初始化 ---
    % a. 初始化权重向量（M=2）
    weights = init_weights(params.N, problem.nObj);

    % b. 初始化邻域 B（无工具箱欧氏距离，修正版）
    % dist_matrix: N x N，dist(i,j)=||w_i - w_j||_2
    D = permute(weights,[1 3 2]) - permute(weights,[3 1 2]); % N x N x M
    dist_matrix = sqrt(sum(D.^2, 3));
    [~, B] = sort(dist_matrix, 2);
    B = B(:, 1:params.T); % 邻域索引（含自身）

    % c. 初始化理想点 z（过滤 NaN/Inf）
    A = vertcat(Pop.Objectives);
    valid = all(isfinite(A), 2);
    if ~any(valid)
        error('MOEA/D 初始化失败：所有初始解的目标函数值均无效（Inf/NaN）。');
    end
    z = min(A(valid, :), [], 1);

    %% 2. 主循环
    for gen = 1:params.T_max
        for i = 1:params.N
            % 1) 父代选择：以 delta 概率在邻域采样，否则全局采样
            if rand < params.delta
                pool = B(i,:);
            else
                pool = 1:params.N;
            end
            if numel(pool) >= 2
                p_indices = pool(randperm(numel(pool), 2));
            else
                % 退化情形保护（N=1）
                p_indices = [i i];
            end
            parent1 = Pop(p_indices(1));
            parent2 = Pop(p_indices(2));

            % 2) 交叉 + 变异
            if rand <= params.pc
                [child_pos, ~] = Crossover(parent1.Position, parent2.Position, problem, params);
            else
                child_pos = parent1.Position;
            end

            child = particle_template;
            child.Position = child_pos;
            child = Mutate(child, problem, params);

            % 3) 评估
            EvalResults_child = feval(problem.objFunc, child.Position, problem);
            child.Objectives = EvalResults_child.Objectives;
            child.Tmax = EvalResults_child.Tmax;

            % 若出现 NaN/Inf，跳过更新（也可在此加入罚函数/可行优先）
            if ~all(isfinite(child.Objectives))
                continue;
            end

            % 4) 更新理想点 z
            z = min(z, child.Objectives);

            % 5) 邻域内更新（Tchebycheff）
            for j = 1:params.T
                neighbor_idx = B(i, j);
                g_old = max(abs(Pop(neighbor_idx).Objectives - z) .* weights(neighbor_idx, :));
                g_new = max(abs(child.Objectives - z) .* weights(neighbor_idx, :));
                if g_new < g_old
                    Pop(neighbor_idx) = child;
                end
            end
        end

        % 日志（每 100 代/起始/结束打印）
        if gen == 1 || gen == params.T_max || mod(gen, 100) == 0
            all_objs = vertcat(Pop.Objectives);
            valid_gen = all(isfinite(all_objs),2);
            if any(valid_gen)
                best_g1 = min(all_objs(valid_gen,1));
                fprintf('迭代: %d/%d, 种群: %d, 当前最优 G1: %.6f\n', ...
                    gen, params.T_max, numel(Pop), best_g1);
            else
                fprintf('迭代: %d/%d, 警告：本代所有个体目标无效。\n', gen, params.T_max);
            end
        end
    end

    %% 3. 输出格式化（与 NSGA-II 兼容）
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
    child = parent; % 默认复制父代


    % 1. 连续部分变异 (Polynomial Mutation)
    % 部署位置
    dim = numel(child.Position.deployment);
    for d = 1:dim
        if rand() <= params.pm
            child.Position.deployment(d) = PolynomialMutation(child.Position.deployment(d), ...
                problem.bounds.deployment(1,d), problem.bounds.deployment(2,d), params.mum);
        end
    end
    % 带宽
    dim_bw = numel(child.Position.bandwidth);
    for d = 1:dim_bw
        if rand() <= params.pm
            child.Position.bandwidth(d) = PolynomialMutation(child.Position.bandwidth(d), ...
                problem.bounds.bandwidth(1,d), problem.bounds.bandwidth(2,d), params.mum);
        end
    end


    % 2. 离散部分变异 (Random Mutation)
    for d = 1:problem.nTerminals
        if rand() <= params.pm % 使用统一变异概率，也可以设置独立概率
            current_choice = child.Position.offloading(d);
            possible_new_nodes = 1:problem.nFogNodes;
            possible_new_nodes(possible_new_nodes == current_choice) = []; % 移除当前选择
            
            if ~isempty(possible_new_nodes) % 确保有其他选项
                child.Position.offloading(d) = possible_new_nodes(randi(numel(possible_new_nodes)));
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
