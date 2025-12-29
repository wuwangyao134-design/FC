% =========================================================================
% NSGA-II 算法主函数 (MyNSGA_II.m)
% =========================================================================
function Pop = MyNSGA_II(problem, params, LastSlotArchive) % 新增 LastSlotArchive 参数
    %% 1. 初始化
    % NSGA-II的particle_template需要Rank和CrowdingDistance字段
    tic;
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'CrowdingDistance', [], 'Rank', []); 
    
    % !!! FIX: 在NSGA_II函数内部重新定义 problem.bounds.deployment
    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2);
    
    for k = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k-1)*2 + 1) = problem.area(1,1); % x_min = 0
        problem.bounds.deployment(2, (k-1)*2 + 1) = problem.area(1,2); % x_max = 100
        problem.bounds.deployment(1, (k-1)*2 + 2) = problem.area(2,1); % y_min = 0
        problem.bounds.deployment(2, (k-1)*2 + 2) = problem.area(2,2); % y_max = 100
    end

    % 调用新的初始化函数，传入 LastSlotArchive
    Pop = InitializePopulation_NSGAII(problem, params, particle_template, LastSlotArchive);
    
    % 评估初始粒子 (Objectives和Tmax将在InitializePopulation_NSGAII中评估)
    % 原始NSGA_II函数中的循环评估逻辑已经移动到InitializePopulation_NSGAII
    % 确保Pop中的所有粒子都已评估
    for i = 1:params.N
        if isempty(Pop(i).Objectives) % 避免重复评估
            EvalResults = feval(problem.objFunc, Pop(i).Position, problem);
            Pop(i).Objectives = EvalResults.Objectives;
            Pop(i).Tmax = EvalResults.Tmax;
        end
    end

    % 对初始种群进行第一次非支配排序和拥挤度距离计算
    Fronts = FindAllFronts(Pop); 
    
    % 在Pop中更新Rank和CrowdingDistance，以便后续选择和变异使用
    for f_idx = 1:numel(Fronts)
        current_front = Fronts{f_idx};
        if isempty(current_front), continue; end
        
        if isscalar(current_front)
            crowd_dist_vals = inf;
        else
            crowd_dist_vals = CalculateCrowdingDistance(vertcat(current_front.Objectives));
        end
        
        for i_mem = 1:numel(current_front)
            % 找到原始Pop中的对应索引并更新其Rank和CrowdingDistance
            % 使用精确匹配，因为在处理惩罚值或浮点数时，Objectives可能完全相同
            original_idx = find(arrayfun(@(p) isequal(p.Position.deployment, current_front(i_mem).Position.deployment) && ...
                                           isequal(p.Position.bandwidth, current_front(i_mem).Position.bandwidth) && ...
                                           isequal(p.Position.offloading, current_front(i_mem).Position.offloading), Pop), 1, 'first');
            if ~isempty(original_idx)
                Pop(original_idx).Rank = f_idx;
                Pop(original_idx).CrowdingDistance = crowd_dist_vals(i_mem);
            end
        end
    end

    %% 2. 主循环
    for gen = 1:params.T_max
        % 1. 选择父代 (Binary Tournament Selection)
        Parents = repmat(Pop(1), params.N, 1);
        selected_indices = BinaryTournamentSelection(Pop, params.N);
        Parents = Pop(selected_indices);

        % 2. 交叉 (Crossover) 和 变异 (Mutation) 生成子代 (Offspring)
        Offspring = repmat(particle_template, params.N, 1);
        for i = 1:2:params.N-1
            parent1 = Parents(i);
            parent2 = Parents(i+1);
            
            [child1_pos, child2_pos] = Crossover(parent1.Position, parent2.Position, problem, params);
            
            child1 = particle_template; child1.Position = child1_pos;
            child2 = particle_template; child2.Position = child2_pos;
            
            child1 = Mutate(child1, problem, params, parent1.CrowdingDistance,gen);
            child2 = Mutate(child2, problem, params, parent2.CrowdingDistance,gen);
            
            EvalResults_child1 = feval(problem.objFunc, child1.Position, problem);
            child1.Objectives = EvalResults_child1.Objectives;
            child1.Tmax = EvalResults_child1.Tmax;
            
            EvalResults_child2 = feval(problem.objFunc, child2.Position, problem);
            child2.Objectives = EvalResults_child2.Objectives;
            child2.Tmax = EvalResults_child2.Tmax;
            
            Offspring(i) = child1;
            Offspring(i+1) = child2;
        end
        if mod(params.N, 2) == 1
            last_parent = Parents(params.N);
            mutated_child = particle_template; % 从模板开始
            mutated_child.Position = last_parent.Position; % 继承最后一个父代的位置
            mutated_child = Mutate(mutated_child, problem, params, last_parent.CrowdingDistance,gen); % 只变异一次
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
        if isempty(R), Pop = []; break; end

        Fronts = FindAllFronts(R);
        
        % 4. 选择下一代种群
        Pop_next = repmat(particle_template, 0, 1);
        current_pop_size = 0;
        
        for f_idx = 1:numel(Fronts)
            current_front = Fronts{f_idx};
            
            if isempty(current_front), continue; end
            
            for i_mem = 1:numel(current_front)
                current_front(i_mem).Rank = f_idx;
            end

            if isscalar(current_front)
                crowd_dist_vals = inf;
            else
                crowd_dist_vals = CalculateCrowdingDistance(vertcat(current_front.Objectives));
            end
            
            for i_mem = 1:numel(current_front)
                current_front(i_mem).CrowdingDistance = crowd_dist_vals(i_mem);
            end
            
            [~, sorted_crowd_idx] = sort([current_front.CrowdingDistance], 'descend');
            current_front_sorted = current_front(sorted_crowd_idx);
            
            num_to_add = min(numel(current_front_sorted), params.N - current_pop_size);
            Pop_next = [Pop_next; current_front_sorted(1:num_to_add)];
            current_pop_size = numel(Pop_next);
            
            if current_pop_size >= params.N
                break; 
            end
        end
        Pop = Pop_next;
        if mod(gen, 100) == 0 || gen == 1 || gen == params.T_max
            fprintf('迭代: %d/%d, 种群大小: %d, 第一前沿大小: %d\n', gen, params.T_max, numel(Pop), numel(Fronts{1}));
        end
        % --- 探索过程可视化 (更新为NSGA-II的种群) ---
        % (这部分保持不变，用于展示探索过程)
        %if mod(gen, 50) == 0 || gen == 1 || gen == params.T_max
        % if gen == params.T_max
        %     figure('Name', ['MyNSGA-II 探索过程 - 迭代 ' num2str(gen)]);
        %     hold on;
        % 
        %     plot(problem.terminalProperties.positions(:,1), problem.terminalProperties.positions(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', '终端');
        % 
        %     all_fog_pos = [];
        %     for i = 1:numel(Pop)
        %         fog_pos_flat = Pop(i).Position.deployment;
        %         fog_pos = reshape(fog_pos_flat, [2, problem.nFogNodes])';
        %         all_fog_pos = [all_fog_pos; fog_pos];
        %     end
        %     plot(all_fog_pos(:,1), all_fog_pos(:,2), '.', 'Color', [0.5 0.5 0.5], 'DisplayName', '种群探索位置');
        % 
        %     if ~isempty(Fronts) && ~isempty(Fronts{1})
        %         elite_fog_pos = [];
        %         for i = 1:numel(Fronts{1})
        %             fog_pos_flat = Fronts{1}(i).Position.deployment;
        %             fog_pos = reshape(fog_pos_flat, [2, problem.nFogNodes])';
        %             elite_fog_pos = [elite_fog_pos; fog_pos];
        %         end
        %         plot(elite_fog_pos(:,1), elite_fog_pos(:,2), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', '第一前沿精英解位置');
        %     end
        % 
        %     xlim(problem.area(1,:));
        %     ylim(problem.area(2,:));
        %     title(['NSGA-II 探索过程 - 迭代 ' num2str(gen)]);
        %     xlabel('X 坐标');
        %     ylabel('Y 坐标');
        %     legend;
        %     grid on;
        %     axis equal;
        %     hold off;
        %     drawnow;
        % end
    end
    toc;
end

% =========================================================================
% 辅助函数 - 初始化种群 (InitializePopulation_NSGAII.m)
% 专为NSGA-II设计，并支持记忆机制
% =========================================================================
function Pop = InitializePopulation_NSGAII(problem, params, particle_template, LastSlotArchive)
    Pop = repmat(particle_template, params.N, 1);
    
    % Define initial deployment perturbation strength for non-inherited particles
    % Using 5% of the overall deployment area range (e.g., 0-100)
    initial_deployment_perturb_strength = 0.05 * (problem.area(2,2) - problem.area(2,1)); 
    
    num_inherited = 0;
    if ~isempty(LastSlotArchive) % If there's an archive from the previous time slot (i.e., not the first time slot)
        nMemory = floor(params.memory_ratio * params.N); % Calculate number of particles to inherit
        if nMemory > 0
            % Ensure not to inherit more elite solutions than actually available
            memory_indices = randperm(numel(LastSlotArchive), min(nMemory, numel(LastSlotArchive)));
            num_inherited = numel(memory_indices);
            
            for i = 1:num_inherited
                inherited_particle = LastSlotArchive(memory_indices(i));
                Pop(i).Position = inherited_particle.Position;
                
                % Inherited particles already have evaluated objectives, no need to re-evaluate
                Pop(i).Objectives = inherited_particle.Objectives;
                Pop(i).Tmax = inherited_particle.Tmax;
            end
            fprintf('Inherited %d elite solutions from the previous time slot for initialization.\n', num_inherited);
        end
    end
    
    % Initialize remaining particles (random initialization, but based on unified reference)
    % This loop starts from (num_inherited + 1) to avoid overwriting inherited particles
    for i = (num_inherited + 1):params.N
        % === Unified Fog Node Deployment Initialization ===
        % The first particle in the *non-inherited* set (if any)
        % will use the exact unified initial fog node position.
        % Subsequent non-inherited particles will be a small perturbation around it.
        if i == (num_inherited + 1) % This is the first particle that is *not* inherited
            deployment_pos = problem.initial_fog_deployment_flat;
        else % Remaining non-inherited particles are perturbed around the unified position
            deployment_pos = problem.initial_fog_deployment_flat + ...
                             initial_deployment_perturb_strength .* (rand(1, problem.nFogNodes * 2) - 0.5); % rand-0.5 for +/- perturbation
            
            % Ensure perturbed deployment positions stay within the global bounds
            deployment_pos = max(deployment_pos, problem.bounds.deployment(1,:));
            deployment_pos = min(deployment_pos, problem.bounds.deployment(2,:));
        end
        
        % Initialize continuous part (bandwidth)
        bandwidth = problem.bounds.bandwidth(1,:) + rand(1, problem.nTerminals) .* (problem.bounds.bandwidth(2,:) - problem.bounds.bandwidth(1,:));
        
        % Initialize discrete part (offloading)
        offloading_plan = randi([1, problem.nFogNodes], 1, problem.nTerminals);
        
        Pop(i).Position.deployment = deployment_pos;
        Pop(i).Position.bandwidth = bandwidth;
        Pop(i).Position.offloading = offloading_plan;
        
        % Evaluate initialized particle
        EvalResults = feval(problem.objFunc, Pop(i).Position, problem);
        Pop(i).Objectives = EvalResults.Objectives;
        Pop(i).Tmax = EvalResults.Tmax;
    end
end


% =========================================================================
% 辅助函数 - 交叉 (Crossover.m)
% 模拟二进制交叉 (SBX) for continuous, 单点交叉 for discrete
% =========================================================================
function [child1_pos, child2_pos] = Crossover(parent1_pos, parent2_pos, problem, params)
    child1_pos = parent1_pos; % 默认复制父代
    child2_pos = parent2_pos;

    % 1. 连续部分交叉 (Simulated Binary Crossover - SBX)
    % 部署位置
    if rand() <= params.pc
        [child1_dep, child2_dep] = SBX(parent1_pos.deployment, parent2_pos.deployment, params.mu);
        child1_pos.deployment = child1_dep;
        child2_pos.deployment = child2_dep;
    end
    % 带宽
    if rand() <= params.pc
        [child1_bw, child2_bw] = SBX(parent1_pos.bandwidth, parent2_pos.bandwidth, params.mu);
        child1_pos.bandwidth = child1_bw;
        child2_pos.bandwidth = child2_bw;
    end

    % 2. 离散部分交叉 (单点交叉 for discrete part)
    if rand() <= params.pc % 离散部分也受pc控制是否交叉
        % 选择一个交叉点，范围是 1 到 nTerminals-1
        crossover_point = randi(problem.nTerminals - 1); % 例如，对于20个终端，交叉点可在1到19之间

        % 执行单点交叉
        temp_offloading1 = [parent1_pos.offloading(1:crossover_point), parent2_pos.offloading(crossover_point+1:end)];
        temp_offloading2 = [parent2_pos.offloading(1:crossover_point), parent1_pos.offloading(crossover_point+1:end)];

        child1_pos.offloading = temp_offloading1;
        child2_pos.offloading = temp_offloading2;
    end
    
    % 边界检查 (确保在变异后仍在范围内)
    % 注意：对于连续变量，SBX的结果可能超出边界，所以这里需要边界检查。
    % 对于离散变量（卸载方案），由于是索引选择，通常不会超出边界，但需要确保值有效（1到nFogNodes）。
    child1_pos.deployment = max(child1_pos.deployment, problem.bounds.deployment(1,:));
    child1_pos.deployment = min(child1_pos.deployment, problem.bounds.deployment(2,:));
    child1_pos.bandwidth = max(child1_pos.bandwidth, problem.bounds.bandwidth(1,:));
    child1_pos.bandwidth = min(child1_pos.bandwidth, problem.bounds.bandwidth(2,:));

    child2_pos.deployment = max(child2_pos.deployment, problem.bounds.deployment(1,:));
    child2_pos.deployment = min(child2_pos.deployment, problem.bounds.deployment(2,:));
    child2_pos.bandwidth = max(child2_pos.bandwidth, problem.bounds.bandwidth(1,:));
    child2_pos.bandwidth = min(child2_pos.bandwidth, problem.bounds.bandwidth(2,:));
    
    % 确保离散变量在有效范围内 (虽然 randi 已经保证了，但交叉可能导致边缘情况，所以加固)
    child1_pos.offloading = max(1, min(problem.nFogNodes, round(child1_pos.offloading)));
    child2_pos.offloading = max(1, min(problem.nFogNodes, round(child2_pos.offloading)));
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
function child = Mutate(parent, problem, params, crowding_distance,gen) % 增加 crowding_distance 参数
    child = parent; % 默认复制父代
    % % 自适应参数（新增）
    adaptive_pm = params.pm;   % 默认使用 params.pm 作为固定值
    adaptive_mum = params.mum; % 默认使用 params.mum 作为固定值
    % 策略：拥挤度越小（越密集），变异强度越大（mum越小，扰动越大）
    % 假设拥挤度距离的范围有一个大致的尺度，这里我们进行归一化或直接映射
    % 避免 crowding_distance 为 inf 或 0
    % 计算自适应参数
    if params.adaptive_enabled
        % 基于进化进度的因子
        gen_factor = 1 - (gen / params.T_max)^2;
        if isinf(crowding_distance)
           crowd_factor = 0.5;  
        else
           crowd_factor = min(1, max(0.01, crowding_distance));
        end
        % 计算自适应参数
        adaptive_pm = params.pm_min + (params.pm_max - params.pm_min) * ...
                      (gen_factor * crowd_factor);
        adaptive_mum = params.mum_max - (params.mum_max - params.mum_min) * ...
                       (gen_factor * crowd_factor);
        
        % 确保在设定范围内
        adaptive_pm = max(params.pm_min, min(params.pm_max, adaptive_pm));
        adaptive_mum = max(params.mum_min, min(params.mum_max, adaptive_mum));
    else
        % 不启用自适应时使用固定参数
        adaptive_pm = params.pm;
        adaptive_mum = params.mum;
    end
    
    % 1. 连续部分变异 (Polynomial Mutation)
    % 部署位置
    dim = numel(child.Position.deployment);
    for d = 1:dim
        if rand() <= adaptive_pm
            child.Position.deployment(d) = PolynomialMutation(child.Position.deployment(d), ...
                problem.bounds.deployment(1,d), problem.bounds.deployment(2,d), adaptive_mum); % 使用 adaptive_mum
        end
    end
    % 带宽
    dim_bw = numel(child.Position.bandwidth);
    for d = 1:dim_bw
        if rand() <= adaptive_pm
            child.Position.bandwidth(d) = PolynomialMutation(child.Position.bandwidth(d), ...
                problem.bounds.bandwidth(1,d), problem.bounds.bandwidth(2,d), adaptive_mum); % 使用 adaptive_mum
        end
    end

    % 2. 离散部分变异 (Random Mutation)
    discrete_pm = 0; % 先给个占位符
    if params.adaptive_enabled
         % gen_factor 在前面的 if 块中已计算
       discrete_pm = params.pm_min + (params.pm_max - params.pm_min) * gen_factor; % 离散 pm 随代数自适应
    else
       discrete_pm = params.pm; % 当不启用自适应时，使用固定 pm
    end

    for d = 1:problem.nTerminals
        if rand() <= discrete_pm% 使用统一变异概率，也可以设置独立概率
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
% FindNonDominated.m - 找到第一非支配前沿
% =========================================================================
function non_dominated_pop = FindNonDominated(pop)
    objectives = vertcat(pop.Objectives);
    nPop = size(objectives, 1);
    is_dominated = false(nPop, 1);
    for i = 1:nPop
        for j = 1:nPop
            if i == j, continue; end
            if all(objectives(j,:) <= objectives(i,:)) && any(objectives(j,:) < objectives(i,:))
                is_dominated(i) = true;
                break;
            % 修复：检查是否是inf/NaN，如果是，则不参与支配判断
            elseif any(isinf(objectives(i,:))) || any(isnan(objectives(i,:)))
                is_dominated(i) = true; % 无效解被视为被支配，不进入非支配前沿
                break;
            end
        end
    end
    non_dominated_pop = pop(~is_dominated);
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
