% =========================================================================
% 纯血 MOPSO 算法主函数 (MOPSO_Pure.m)
% =========================================================================
function Archive = MOPSO(problem, params, ~)
    %% 1. 初始化
    % !!! FIX: 在HMD_MOPSO函数内部重新定义 problem.bounds.deployment
    % 确保它是一个正确的上下界矩阵，而不是依赖外部可能被误导的定义
    problem.bounds.deployment = zeros(2, problem.nFogNodes * 2); % 行数固定为2
    for k_node = 1:problem.nFogNodes
        problem.bounds.deployment(1, (k_node-1)*2 + 1) = problem.area(1,1); % x_min = 0
        problem.bounds.deployment(2, (k_node-1)*2 + 1) = problem.area(1,2); % x_max = 100
        problem.bounds.deployment(1, (k_node-1)*2 + 2) = problem.area(2,1); % y_min = 0
        problem.bounds.deployment(2, (k_node-1)*2 + 2) = problem.area(2,2); % y_max = 100
    end

    [Swarm, pbest] = InitializePopulation(problem, params);
    Archive = [];
    % stagnation_counter = 0; % 暂时不启用重启机制
    % last_archive_size = 0;   % 暂时不启用重启机制
    
    for i = 1:params.N
        EvalResults = feval(problem.objFunc, Swarm(i).Position, problem);
        Swarm(i).Objectives = EvalResults.Objectives;
        Swarm(i).Tmax = EvalResults.Tmax;
        pbest(i) = Swarm(i);
    end


    Archive = UpdateArchive(Archive, Swarm, params);
    % last_archive_size = numel(Archive); % 暂时不启用重启机制
    
    %% 2. 主循环
    for t = 1:params.T_max
        w = AdaptiveInertiaWeight(Archive, params, t, params.T_max);
        Leaders = SelectLeaders(Archive, params);
        
        % if isempty(Leaders)
        %    % 如果存档为空，从当前种群的pbest中随机选择一个作为所有粒子的全局领导者
        %    rand_idx = randi(numel(pbest));
        %    global_leader = pbest(rand_idx);
        %    % 确保 Leaders 是一个大小为 N 的数组，每个元素都是这个 global_leader
        %    Leaders = repmat(global_leader, 1, params.N);
        %    warning('存档为空，所有粒子随机选择一个pbest作为共同领导者。');
        % end
        if isempty(Leaders)
            Leaders = pbest(randperm(numel(pbest), params.N));  % 随机分配不同的pbest作为领导者
        end

        for i = 1:params.N
            Swarm(i) = UpdateParticle(Swarm(i), pbest(i), Leaders(i), problem, params, w);
            Swarm(i) = Mutate(Swarm(i), problem, params); % 重新启用Mutate
        end
        
        for i = 1:params.N
            EvalResults = feval(problem.objFunc, Swarm(i).Position, problem);
            Swarm(i).Objectives = EvalResults.Objectives;
            Swarm(i).Tmax = EvalResults.Tmax;
        end
        
        pbest = UpdatePbest(pbest, Swarm);
        Archive = UpdateArchive(Archive, Swarm, params);
        if mod(t, 100) == 0 || t == 1 || t == params.T_max
            fprintf('迭代: %d/%d, 存档大小: %d, 惯性权重: %.4f\n', t, params.T_max, numel(Archive), w);
        end
        
        % --- 探索过程可视化 ---
        % if params.T_max
        %     figure('Name', ['部署探索过程 - 迭代 ' num2str(t)]);
        %     hold on;
        % 
        %     plot(problem.terminalProperties.positions(:,1), problem.terminalProperties.positions(:,2), 'bo', 'MarkerFaceColor', 'b', 'DisplayName', '终端');
        % 
        %     all_fog_pos = [];
        %     for i = 1:params.N
        %         fog_pos_flat = Swarm(i).Position.deployment;
        %         fog_pos = reshape(fog_pos_flat, [2, problem.nFogNodes])';
        %         all_fog_pos = [all_fog_pos; fog_pos];
        %     end
        %     plot(all_fog_pos(:,1), all_fog_pos(:,2), '.', 'Color', [0.5 0.5 0.5], 'DisplayName', '粒子群探索位置');           
        %     xlim(problem.area(1,:));
        %     ylim(problem.area(2,:));
        %     title(['部署探索过程 - 迭代 ' num2str(t)]);
        %     xlabel('X 坐标');
        %     ylabel('Y 坐标');
        %     legend('show'); % 确保图例显示所有DisplayName
        %     grid on;
        %     axis equal;
        %     hold off;
        %     drawnow;
        % end
    end
end


% =========================================================================
% 辅助函数 - 初始化种群 (InitializePopulation_MOPSO.m)
% 针对纯血MOPSO，移除记忆机制
% =========================================================================
% =========================================================================
% 辅助函数 - 初始化种群 (InitializePopulation_MOPSO.m)
% 针对纯血MOPSO，移除记忆机制
% =========================================================================
function [Swarm, pbest] = InitializePopulation(problem, params) % 确保函数名是 InitializePopulation_MOPSO

    % 粒子模板，纯血MOPSO通常不直接在粒子中存储Rank和CrowdingDistance
    particle_template = struct('Position', [], 'Velocity', [], 'Objectives', [], 'Tmax', [], 'pbest', []); 
    Swarm = repmat(particle_template, params.N, 1);
    
    % 定义初始部署扰动强度
    initial_deployment_perturb_strength = 0.05 * (problem.area(2,2) - problem.area(2,1));
    
    for i = 1:params.N
        % === 修正：基于统一的初始雾节点位置进行初始化 ===
        if i == 1 % 第一个粒子使用统一的初始雾节点位置
            deployment_pos = problem.initial_fog_deployment_flat;
        else % 剩余粒子在统一位置附近进行小范围随机扰动
            deployment_pos = problem.initial_fog_deployment_flat + ...
                             initial_deployment_perturb_strength .* (rand(1, problem.nFogNodes * 2) - 0.5);
            
            % 确保扰动后的部署位置仍在有效边界内 (使用 problem.bounds.deployment 定义的全局边界)
            % 注意：MOPSO通常会在 UpdateParticle 中才进行最终边界钳制，但这里初始化时进行一次是好的
            % 这里的 problem.bounds.deployment 应该确保在 MOPSO_Pure 主函数中已经定义并传递到 problem
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
        Swarm(i).Position.deployment = deployment_pos;
        Swarm(i).Position.bandwidth = bandwidth;
        Swarm(i).Position.offloading = offloading_plan;
        
        % 初始化速度
        Swarm(i).Velocity.deployment = zeros(1, problem.nFogNodes*2);
        Swarm(i).Velocity.bandwidth = zeros(1, problem.nTerminals);
        Swarm(i).Velocity.offloading = cell(1, problem.nTerminals); % 离散部分的概率分布
        for j = 1:problem.nTerminals
            Swarm(i).Velocity.offloading{j} = ones(1, problem.nFogNodes) / problem.nFogNodes; % 初始化为均匀概率
        end
        
        % pbest 在主函数中评估后初始化
    end
    
    pbest = Swarm; % 初始时，pbest 是 Swarm 的副本 (但在主函数中会再次赋值Objectives等)
end
% =========================================================================
% 辅助函数 - 混合更新粒子 (UpdateParticle.m)
% =========================================================================
function particle = UpdateParticle(particle, pbest, leader, problem, params, w)
    r1 = rand(); 
    r2 = rand();
    
    % 连续部分更新保持不变
    particle.Velocity.deployment = w*particle.Velocity.deployment + ...
        params.c1*r1*(pbest.Position.deployment - particle.Position.deployment) + ...
        params.c2*r2*(leader.Position.deployment - particle.Position.deployment);
    particle.Position.deployment = particle.Position.deployment + particle.Velocity.deployment;
    
    particle.Velocity.bandwidth = w*particle.Velocity.bandwidth + ...
        params.c1*r1*(pbest.Position.bandwidth - particle.Position.bandwidth) + ...
        params.c2*r2*(leader.Position.bandwidth - particle.Position.bandwidth);
    particle.Position.bandwidth = particle.Position.bandwidth + particle.Velocity.bandwidth;
    
    % 边界检查保持不变
    particle.Position.deployment = max(particle.Position.deployment, problem.bounds.deployment(1,:));
    particle.Position.deployment = min(particle.Position.deployment, problem.bounds.deployment(2,:));
    particle.Position.bandwidth = max(particle.Position.bandwidth, problem.bounds.bandwidth(1,:));
    particle.Position.bandwidth = min(particle.Position.bandwidth, problem.bounds.bandwidth(2,:));
    
    % --- 改进的离散部分更新 ---
    for d = 1:problem.nTerminals
        % 初始化权重向量
        weights = zeros(1, problem.nFogNodes);
        
        % 当前位置的权重（惯性项）
        current_pos = particle.Position.offloading(d);
        weights(current_pos) = w;
        
        % 个人最优位置的权重
        pbest_pos = pbest.Position.offloading(d);
        weights(pbest_pos) = weights(pbest_pos) + params.c1 * r1;
        
        % 全局最优位置的权重
        leader_pos = leader.Position.offloading(d);
        weights(leader_pos) = weights(leader_pos) + params.c2 * r2;
        
        % 归一化权重（避免除零）
        if sum(weights) > 0
            weights = weights / sum(weights);
        else
            weights = ones(1, problem.nFogNodes) / problem.nFogNodes;
        end
        
        % 使用softmax进一步平滑概率分布（可选）
        % weights = exp(weights - max(weights)); % 减去最大值避免数值溢出
        % weights = weights / sum(weights);
        
        % 使用轮盘赌选择新的卸载决策
        cumsum_weights = cumsum(weights);
        r = rand();
        selected_fog = find(r <= cumsum_weights, 1, 'first');
        
        % 更新卸载决策和速度
        particle.Position.offloading(d) = selected_fog;
        particle.Velocity.offloading{d} = weights; % 存储概率分布作为"速度"
    end
end



% =========================================================================
% 辅助函数 - 更新个体最优 (UpdatePbest_MOPSO.m)
% 注意：我们将pbest直接集成到Swarm结构体中，简化参数传递
% =========================================================================
function pbest = UpdatePbest(pbest, Swarm)
    for i = 1:numel(Swarm)
        obj_s = Swarm(i).Objectives;
        obj_p = pbest(i).Objectives;
        
        if isempty(obj_s) || any(isinf(obj_s)) || any(isnan(obj_s)), continue; end
        
        if all(obj_s <= obj_p) && any(obj_s < obj_p)
            pbest(i) = Swarm(i);
        elseif ~(all(obj_p <= obj_s) && any(obj_p < obj_s))
            if rand() < 0.5, pbest(i) = Swarm(i); end
        end
    end
end



% =========================================================================
% 辅助函数 - 更新外部存档 (UpdateArchive_MOPSO.m)
% =========================================================================
function Archive = UpdateArchive(Archive, Swarm, params)
    temp_pop_list = [Archive; Swarm];
    
    valid_indices = true(numel(temp_pop_list), 1);
    for i = 1:numel(temp_pop_list)
        if isempty(temp_pop_list(i).Objectives) || any(isinf(temp_pop_list(i).Objectives)) || any(isnan(temp_pop_list(i).Objectives))
            valid_indices(i) = false;
        end
    end
    temp_pop = temp_pop_list(valid_indices);
    if isempty(temp_pop), Archive = []; return; end
    
    % 在 MOPSO 中，FindNonDominated 通常用于筛选存档
    Archive = FindNonDominated(temp_pop); 
    
    if numel(Archive) > params.ArchiveSize
        % TruncateArchive 函数需要确保其逻辑与MOPSO的网格机制匹配
        Archive = TruncateArchive(Archive, params.ArchiveSize, params.nGrid);
    end
end


% =========================================================================
% 辅助函数 - 选择领导者 (SelectLeaders_MOPSO.m)
% MOPSO 版本的领导者选择，通常基于网格密度
% =========================================================================
function Leaders = SelectLeaders(Archive, params)
    if isempty(Archive)
        Leaders = []; 
        return;
    end
    
    N = params.N;
    nGrid = params.nGrid;
    
    Leaders = repmat(Archive(1), N, 1); % 初始化 Leaders 结构，与粒子群大小 N 相同

    % 基于网格创建
    [Grid, ~] = CreateGrid(Archive, nGrid); % 调用MOPSO专用的CreateGrid

    occupied_grids_indices = find([Grid.Count] > 0);
    if isempty(occupied_grids_indices)
        % 如果没有被占据的网格，则随机选择存档中的任何个体作为领导者
        % 或者更简单地，直接从Archive中随机选N个领导者
        Leaders = Archive(randi(numel(Archive), 1, N));
        return;
    end
    
    % 密度反比适应度选择
    density = [Grid(occupied_grids_indices).Count];
    fitness = 1 ./ (density + 1); % 密度越低，适应度越高
    probabilities = fitness / sum(fitness); % 归一化概率
    
    for i = 1:N
        % 轮盘赌选择一个网格
        chosen_grid_idx_in_occupied = find(rand() <= cumsum(probabilities), 1, 'first');
        chosen_grid_idx = occupied_grids_indices(chosen_grid_idx_in_occupied);
        
        % 从被选中的网格中随机选择一个粒子作为领导者
        leader_member_idx = Grid(chosen_grid_idx).Members(randi([1, Grid(chosen_grid_idx).Count]));
        
        % 复制领导者信息到 Leaders 数组
        % 确保只复制 MOPSO 粒子结构中需要的字段
        Leaders(i).Position = Archive(leader_member_idx).Position;
        Leaders(i).Velocity = Archive(leader_member_idx).Velocity; % 领导者的速度也需要
        Leaders(i).Objectives = Archive(leader_member_idx).Objectives;
        Leaders(i).Tmax = Archive(leader_member_idx).Tmax;
        Leaders(i).pbest = Archive(leader_member_idx).pbest; % 领导者也应该有pbest字段
    end
end

% =========================================================================
% 辅助函数 - 自适应惯性权重 (AdaptiveInertiaWeight_MOPSO.m)
% 通常在MOPSO中惯性权重是基于迭代次数线性递减，或根据多样性调整
% 这里保留了原HMD版本中的多样性调整逻辑
function w = AdaptiveInertiaWeight(Archive, params, t, T_max)
    if isempty(Archive)
        % 如果存档为空，给予较大的惯性权重，促进探索
        w = params.w_max; 
        return;
    end
    
    % 1. 基于迭代次数的线性递减 (经典且有效)
    % 确保 w_max > w_min，通常 w_max = 0.9, w_min = 0.4
    w_linear = params.w_max - (params.w_max - params.w_min) * (t / T_max);
    
    % 2. 基于多样性的调整 (微调)
    [Grid, ~] = CreateGrid(Archive, params.nGrid);
    G_occupied = sum([Grid.Count] > 0);
    diversity_metric = G_occupied / numel(Grid); % 被占据网格的比例，范围 [0, 1]

    % 定义一个多样性对惯性权重的“影响因子”
    % 当多样性高时，影响因子可能略大于1，使 w 稍微增大，鼓励探索
    % 当多样性低时，影响因子可能略小于1，使 w 稍微减小，鼓励开发
    % 示例：可以将 diversity_metric 线性映射到 [0.8, 1.2] 或 [0.9, 1.1] 的范围
    % 这里的映射需要根据你的实验效果来调整
    diversity_factor = 1 + (diversity_metric - 0.5) * 0.4; % 示例: 0.5 -> 1, 0 -> 0.8, 1 -> 1.2
                                                          % 调整 0.4 来控制波动强度
    
    % 将线性递减和多样性因子结合
    w = w_linear * diversity_factor;

end


% =========================================================================
% 辅助函数 - TruncateArchive_MOPSO.m (HMD-MOPSO版本)
% 保持与MOPSO的外部存档截断机制一致，基于网格
% =========================================================================
function TruncatedArchive = TruncateArchive(Archive, ArchiveSize, nGrid)
    while numel(Archive) > ArchiveSize
        [Grid, ~] = CreateGrid(Archive, nGrid);
        
        grid_counts = [Grid.Count];
        % 确保 grid_counts 非空，且至少有一个正数 Count
        if isempty(grid_counts) || all(grid_counts == 0)
            % 如果所有网格都为空，或者只有0个计数，则随机移除一个
            if numel(Archive) > 0
                Archive(randi(numel(Archive))) = [];
            end
            continue; % 继续循环直到大小达标
        end
        
        % 找到最拥挤的网格
        [~, max_grid_idx] = max(grid_counts(:));
        
        members_in_grid = Grid(max_grid_idx).Members;
        % 从最拥挤的网格中随机移除一个成员
        idx_to_remove_in_members = randi([1, numel(members_in_grid)]);
        idx_to_remove_in_archive = members_in_grid(idx_to_remove_in_members);
        
        Archive(idx_to_remove_in_archive) = [];
    end
    TruncatedArchive = Archive;
end

% =========================================================================
% 辅助函数 - CreateGrid_MOPSO.m
% 保持与MOPSO的网格划分机制一致
% =========================================================================
function [Grid, GridSubIndex] = CreateGrid(pop, nGrid)
    objectives = vertcat(pop.Objectives);
    min_obj = min(objectives, [], 1);
    max_obj = max(objectives, [], 1);
    
    for j = 1:size(objectives, 2)
        if min_obj(j) == max_obj(j)
            max_obj(j) = min_obj(j) + 1; % 避免范围为0
        end
    end
    
    Grid = repmat(struct('Count', 0, 'Members', []), nGrid, nGrid); 
    GridSubIndex = zeros(numel(pop), size(objectives, 2));
    for i = 1:numel(pop)
        for j = 1:size(objectives, 2)
            edge = linspace(min_obj(j), max_obj(j), nGrid + 1);
            sub_idx = find(objectives(i,j) <= edge, 1, 'first') - 1;
            if isempty(sub_idx) || sub_idx == 0
                sub_idx = 1;
            end
            if sub_idx > nGrid
                sub_idx = nGrid;
            end
            GridSubIndex(i, j) = sub_idx;
        end
        if size(objectives, 2) >= 2
            idx1 = GridSubIndex(i, 1);
            idx2 = GridSubIndex(i, 2);
            Grid(idx1, idx2).Count = Grid(idx1, idx2).Count + 1;
            Grid(idx1, idx2).Members = [Grid(idx1, idx2).Members, i];
        else 
            idx1 = GridSubIndex(i, 1);
            Grid(idx1, 1).Count = Grid(idx1, 1).Count + 1; 
            Grid(idx1, 1).Members = [Grid(idx1, 1).Members, i];
        end
    end
end

% =========================================================================
% 辅助函数 - 变异算子 (Mutate_MOPSO.m) - 针对纯血MOPSO
% 采用标准多项式变异或高斯变异
% =========================================================================
function particle = Mutate(particle, problem, params)
    % 变异概率
    if rand() > params.pm_mopso % 使用专门的变异概率参数，确保统一
        return; % 不发生变异
    end
    
    % --- 对所有连续变量进行多项式变异 ---
    % 部署位置
    dim_dep = numel(particle.Position.deployment);
    for d = 1:dim_dep
        particle.Position.deployment(d) = PolynomialMutation(particle.Position.deployment(d), ...
            problem.bounds.deployment(1,d), problem.bounds.deployment(2,d), params.mum_mopso); % 使用MOPSO专用参数
    end
    % 带宽
    dim_bw = numel(particle.Position.bandwidth);
    for d = 1:dim_bw
        particle.Position.bandwidth(d) = PolynomialMutation(particle.Position.bandwidth(d), ...
            problem.bounds.bandwidth(1,d), problem.bounds.bandwidth(2,d), params.mum_mopso); % 使用MOPSO专用参数
    end
    
    % --- 对离散变量进行随机变异 ---
    for d = 1:problem.nTerminals
        current_choice = particle.Position.offloading(d);
        possible_new_nodes = 1:problem.nFogNodes;
        possible_new_nodes(possible_new_nodes == current_choice) = []; % 移除当前选择
            
        if ~isempty(possible_new_nodes) % 确保有其他选项
            particle.Position.offloading(d) = possible_new_nodes(randi(numel(possible_new_nodes)));
        end
    end

    % 再次进行边界检查 (确保变异后仍在范围内)
    particle.Position.deployment = max(particle.Position.deployment, problem.bounds.deployment(1,:));
    particle.Position.deployment = min(particle.Position.deployment, problem.bounds.deployment(2,:));
    particle.Position.bandwidth = max(particle.Position.bandwidth, problem.bounds.bandwidth(1,:));
    particle.Position.bandwidth = min(particle.Position.bandwidth, problem.bounds.bandwidth(2,:));
end

% 辅助函数：Polynomial Mutation (MOPSO专用，确保参数独立)
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