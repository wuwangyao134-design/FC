% =========================================================================
% NSGA-II 算法主函数 (NSGA_II.m)
% =========================================================================
function Pop = NSGA_II(problem, params,~)
    %% 1. 初始化
    % NSGA-II的particle_template需要Rank和CrowdingDistance字段
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'CrowdingDistance', [], 'Rank', []); 

    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2); % <--- 修正：行数固定为2，且使用problem.nFogNodes
    
    for k = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k-1)*2 + 1) = problem.area(1,1); % x_min = 0
        problem.bounds.deployment(2, (k-1)*2 + 1) = problem.area(1,2); % x_max = 100
        problem.bounds.deployment(1, (k-1)*2 + 2) = problem.area(2,1); % y_min = 0
        problem.bounds.deployment(2, (k-1)*2 + 2) = problem.area(2,2); % y_max = 100
    end

    Pop = repmat(particle_template, params.N, 1); % 初始化 Pop
    initial_deployment_perturb_strength = 0.05 * (problem.area(2,2) - problem.area(2,1)); 


    for i = 1:params.N
        % === 修正：基于统一的初始雾节点位置进行初始化 ===
        if i == 1 % 第一个粒子使用统一的初始雾节点位置
             deployment_pos = problem.initial_fog_deployment_flat;
        else % 剩余粒子在统一位置附近进行小范围随机扰动
             deployment_pos = problem.initial_fog_deployment_flat + ...
                         initial_deployment_perturb_strength .* (rand(1, problem.nFogNodes * 2) - 0.5); % rand-0.5 使得扰动在正负之间
        
             % 确保扰动后的部署位置仍在有效边界内 (使用 problem.bounds.deployment 定义的全局边界)
             deployment_pos = max(deployment_pos, problem.bounds.deployment(1,:));
             deployment_pos = min(deployment_pos, problem.bounds.deployment(2,:));
        end
             % === 核心修正：删除下面这行，因为它会覆盖上面 if/else 块的计算结果 ===
             % deployment_pos = problem.bounds.deployment(1,:) + rand(1, problem.nFogNodes * 2) .* (problem.bounds.deployment(2,:) - problem.bounds.deployment(1,:));
    
             % 初始化连续部分 (带宽)
             bandwidth = problem.bounds.bandwidth(1,:) + rand(1, problem.nTerminals) .* (problem.bounds.bandwidth(2,:) - problem.bounds.bandwidth(1,:));
    
             % 初始化离散部分 (卸载)
             offloading_plan = randi([1, problem.nFogNodes], 1, problem.nTerminals);
    
             % 将计算好的部署位置、带宽和卸载方案赋值给当前粒子
             Pop(i).Position.deployment = deployment_pos;
             Pop(i).Position.bandwidth = bandwidth;
             Pop(i).Position.offloading = offloading_plan;
    
             % 评估初始粒子
             EvalResults = feval(problem.objFunc, Pop(i).Position, problem);
             Pop(i).Objectives = EvalResults.Objectives;
             Pop(i).Tmax = EvalResults.Tmax;
    end


    % 对初始种群进行第一次非支配排序和拥挤度距离计算
    Fronts = FindAllFronts(Pop); 
    
    for f_idx = 1:numel(Fronts)
        current_front = Fronts{f_idx};
        
        if isempty(current_front), continue; end
        
        % 计算拥挤度距离
        if isscalar(current_front) % 修正：使用numel而不是isscalar
            crowd_dist_vals = inf;
        else
            crowd_dist_vals = CalculateCrowdingDistance(vertcat(current_front.Objectives));
        end
        
        % 赋值Rank和CrowdingDistance给原Pop中的对应粒子
        for i_mem = 1:numel(current_front)
            original_idx = find(arrayfun(@(p) isequal(p.Position, current_front(i_mem).Position) && ...
                                           isequal(p.Objectives, current_front(i_mem).Objectives), Pop), 1, 'first');
            if ~isempty(original_idx)
                Pop(original_idx).Rank = f_idx;
                Pop(original_idx).CrowdingDistance = crowd_dist_vals(i_mem);
            end
        end
    end


    %% 2. 主循环
    for gen = 1:params.T_max
        % 1. 选择父代 (Binary Tournament Selection)
        Parents = repmat(Pop(1), params.N, 1); % 初始化 Parents 结构
        selected_indices = BinaryTournamentSelection(Pop, params.N); % 批量选择所有父代的索引
        Parents = Pop(selected_indices); % 根据索引获取父代粒子


        % 2. 交叉 (Crossover) 和 变异 (Mutation) 生成子代 (Offspring)
        Offspring = repmat(Pop(1), params.N, 1); % 初始化 Offspring 结构
        for i = 1:2:params.N-1 % 每次生成两个子代
            parent1 = Parents(i);
            parent2 = Parents(i+1);
            
            % 交叉操作 (Simulated Binary Crossover for continuous parts)
            if rand() <= params.pc
                [child1_pos, child2_pos] = Crossover(parent1.Position, parent2.Position, problem, params);
                child1 = particle_template; child1.Position = child1_pos; % 创建新粒子，只赋Position
                child2 = particle_template; child2.Position = child2_pos;
            else
                % 不交叉则直接复制父代的位置
                child1 = particle_template; child1.Position = parent1.Position;
                child2 = particle_template; child2.Position = parent2.Position;
            end
            
            % 变异操作 (Polynomial Mutation for continuous, custom for discrete)
            child1 = Mutate(child1, problem, params);
            child2 = Mutate(child2, problem, params);
            
            % 评估子代 (在交叉和变异后)
            EvalResults_child1 = feval(problem.objFunc, child1.Position, problem);
            child1.Objectives = EvalResults_child1.Objectives;
            child1.Tmax = EvalResults_child1.Tmax;
            
            EvalResults_child2 = feval(problem.objFunc, child2.Position, problem);
            child2.Objectives = EvalResults_child2.Objectives;
            child2.Tmax = EvalResults_child2.Tmax;
            
            Offspring(i) = child1;
            Offspring(i+1) = child2;
        end
        % 针对params.N为奇数的情况，最后一个父代可能没有配对，这里处理最后一个父代
        if mod(params.N, 2) == 1
            last_parent = Parents(params.N);
            mutated_child = particle_template;
            mutated_child.Position = last_parent.Position; % 继承最后一个父代的位置
            mutated_child = Mutate(mutated_child, problem, params); % 变异
            EvalResults_last = feval(problem.objFunc, mutated_child.Position, problem);
            mutated_child.Objectives = EvalResults_last.Objectives;
            mutated_child.Tmax = EvalResults_last.Tmax;
            Offspring(params.N) = mutated_child;
        end


        % 3. 合并父代和子代，进行非支配排序和拥挤度距离计算
        R = [Pop; Offspring];
        
        valid_indices_R = true(numel(R), 1);
        for i = 1:numel(R)
            if isempty(R(i).Objectives) || any(isinf(R(i).Objectives)) || any(isnan(R(i).Objectives))
                valid_indices_R(i) = false;
            end
        end
        R = R(valid_indices_R);
        if isempty(R), Pop = []; break; end % 如果R为空，算法中止


        % 找到所有非支配前沿并计算拥挤度距离
        Fronts = FindAllFronts(R);
        
        % 4. 选择下一代种群
        Pop_next = repmat(particle_template, 0, 1); % 初始化空数组
        current_pop_size = 0;
        
        for f_idx = 1:numel(Fronts)
            current_front = Fronts{f_idx};
            
            if isempty(current_front), continue; end
            
            % 赋值Rank
            for i_mem = 1:numel(current_front)
                current_front(i_mem).Rank = f_idx;
            end


            % 计算拥挤度距离
            if isscalar(current_front)
                crowd_dist_vals = inf;
            else
                crowd_dist_vals = CalculateCrowdingDistance(vertcat(current_front.Objectives));
            end
            
            for i_mem = 1:numel(current_front)
                current_front(i_mem).CrowdingDistance = crowd_dist_vals(i_mem);
            end
            
            % 根据拥挤度距离降序排序
            [~, sorted_crowd_idx] = sort([current_front.CrowdingDistance], 'descend');
            current_front_sorted = current_front(sorted_crowd_idx);
            
            num_to_add = min(numel(current_front_sorted), params.N - current_pop_size);
            Pop_next = [Pop_next; current_front_sorted(1:num_to_add)];
            current_pop_size = numel(Pop_next);
            
            if current_pop_size >= params.N
                break; 
            end
        end
        Pop = Pop_next; % 更新为下一代种群

        if mod(gen, 100) == 0 || gen == 1 || gen == params.T_max
            fprintf('迭代: %d/%d, 种群大小: %d, 第一前沿大小: %d\n', gen, params.T_max, numel(Pop), numel(Fronts{1}));
        end
    end
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
% 辅助函数 - 二元锦标赛选择 (BinaryTournamentSelection.m)
% =========================================================================
function selected_idx = BinaryTournamentSelection(Pop, k)
    % 从Pop中进行k次二元锦标赛选择，每次选出1个最好的个体
    selected_idx = zeros(1, k);
    nPop = numel(Pop);
    
    for i = 1:k
        % 随机选择两个个体
        p1_idx = randi(nPop);
        p2_idx = randi(nPop);
        
        % 确保选择两个不同的个体
        while p1_idx == p2_idx
            p2_idx = randi(nPop);
        end
        
        p1 = Pop(p1_idx);
        p2 = Pop(p2_idx);
        
        % 比较非支配等级（Rank），等级越低越好
        if p1.Rank < p2.Rank
            selected_idx(i) = p1_idx;
        elseif p2.Rank < p1.Rank
            selected_idx(i) = p2_idx;
        else % 如果Rank相同，比较拥挤度距离，拥挤度距离越大越好
            if p1.CrowdingDistance > p2.CrowdingDistance
                selected_idx(i) = p1_idx;
            else
                selected_idx(i) = p2_idx;
            end
        end
    end
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