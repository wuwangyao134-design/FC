function Pop = MOPO(problem, params, ~)
    %% 1. 初始化
    tic;
    % 统一粒子模板
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'Rank', [], 'Distance', []); 
    
    % 算法特有参数
    N = params.N;          % 种群大小
    K = 100;               % 外部存档容量 
    Archive = [];          % 初始化外部存档 
    
    % ======================== 核心修复代码开始 ========================
    % 确保部署边界字段存在 (解决 "无法识别的字段名称 deployment" 报错)
    if ~isfield(problem.bounds, 'deployment')
        % 根据场景区域 [x_min x_max; y_min y_max] 构造边界
        x_min = problem.area(1,1); x_max = problem.area(1,2);
        y_min = problem.area(2,1); y_max = problem.area(2,2);
        
        % 构造 2 x (2*nFogNodes) 的边界矩阵，每一列对应一个坐标维度
        lower_bound_row = repmat([x_min, y_min], 1, problem.nFogNodes);
        upper_bound_row = repmat([x_max, y_max], 1, problem.nFogNodes);
        problem.bounds.deployment = [lower_bound_row; upper_bound_row];
    end
    
    % 提取边界行向量以确保后面运算维度兼容
    LB_dep = problem.bounds.deployment(1, :);
    UB_dep = problem.bounds.deployment(2, :);
    % ======================== 核心修复代码结束 ========================
    % --- 新增：定义各部分维度的长度，防止报错 ---
    len_dep = numel(LB_dep); 
    len_bw = numel(problem.bounds.bandwidth(1,:));

    % 初始化种群位置
    Pop = repmat(particle_template, N, 1);
    for i = 1:N
        % 使用修正后的 LB_dep 和 UB_dep 进行随机初始化
        Pop(i).Position.deployment = LB_dep + rand(1, problem.nFogNodes*2) .* (UB_dep - LB_dep);
        
        Pop(i).Position.bandwidth = problem.bounds.bandwidth(1,:) + rand(1, problem.nTerminals) .* ...
                                    (problem.bounds.bandwidth(2,:) - problem.bounds.bandwidth(1,:));
        Pop(i).Position.offloading = randi([1, problem.nFogNodes], 1, problem.nTerminals);
    end
    
    % 初始评估
    EvalRes = feval(problem.objFunc, Pop, problem);
    for i = 1:N
        Pop(i).Objectives = EvalRes.Objectives(i,:);
        Pop(i).Tmax = EvalRes.Tmax(i);
    end
    %% 2. 主循环
    for t = 1:params.T_max
        % 更新当前迭代的最优个体 (从存档或当前种群中通过轮盘赌选出) [cite: 610]
        if isempty(Archive)
            Xbest = Pop(randi(N)); 
        else
            Xbest = Archive(randi(numel(Archive))); 
        end
        
        Xmean = CalculateMeanPosition(Pop); % 计算种群平均位置 [cite: 281]
        
       % 遍历更新每个 parrot 的位置
        for i = 1:N
            St = randi([1, 4]); % 随机选择一种行为策略
            
            % 1. 处理连续部分 (部署与带宽) - 保留原有的 Lévy 飞行机制
            pos_cont = [Pop(i).Position.deployment, Pop(i).Position.bandwidth];
            best_cont = [Xbest.Position.deployment, Xbest.Position.bandwidth];
            
            % 原有行为逻辑应用于连续变量
            L = Levy(numel(pos_cont));
            new_cont = pos_cont + (pos_cont - best_cont) .* L; % 示例行为1
            
            % 2. 核心修改：处理离散部分 (卸载决策) - 使用标准离散算子
            % 策略：根据不同的 St 行为，执行不同强度的离散变异或交叉
            current_off = Pop(i).Position.offloading;
            best_off = Xbest.Position.offloading;
            
            switch St
                case 1 % 寻食行为 -> 高强度离散变异 (大范围搜索)
                    new_off = DiscreteMutate(current_off, problem.nFogNodes, 0.3);
                case 2 % 停留行为 -> 局部离散变异 (精细微调)
                    new_off = DiscreteMutate(current_off, problem.nFogNodes, 0.05);
                case 3 % 交流行为 -> 均匀交叉 (向优秀解学习)
                    new_off = DiscreteCrossover(current_off, best_off);
                case 4 % 恐惧陌生人行为 -> 随机重新初始化部分决策位
                    new_off = DiscreteMutate(current_off, problem.nFogNodes, 0.5);
            end
            
            % 3. 合并并更新个体位置
            Pop(i).Position.deployment = max(min(new_cont(1:len_dep), UB_dep), LB_dep);
            Pop(i).Position.bandwidth = max(min(new_cont(len_dep+1:end), ...
                                        problem.bounds.bandwidth(2,:)), problem.bounds.bandwidth(1,:));
            Pop(i).Position.offloading = new_off;
        end
        
        % 评估新位置
        EvalRes = feval(problem.objFunc, Pop, problem);
        for i = 1:N
            Pop(i).Objectives = EvalRes.Objectives(i,:);
            Pop(i).Tmax = EvalRes.Tmax(i);
        end
        
        % 3. 环境选择与存档更新 [cite: 28, 551]
        % 合并、非支配排序
        Combined = [Archive; Pop];
        valid_mask = arrayfun(@(x) ~isempty(x.Objectives), Combined);
        Combined = Combined(valid_mask);
        
        Fronts = FindAllFronts(Combined);
        if ~isempty(Fronts)
            Archive = Fronts{1}; % 仅保留第一层非支配解作为候选存档
        end
        
        % 4. 自适应椭圆分割 (AES) 筛选存档 [cite: 24, 505, 608]
        if numel(Archive) > K
            Archive = AdaptiveEllipticalSelection(Archive, K); 
        end
        
        if mod(t, 100) == 0
            fprintf('MOPO 迭代: %d/%d, 存档规模: %d\n', t, params.T_max, numel(Archive));
        end
    end
    
    % 返回最终种群 (若存档不足 N，从 Pop 补齐)
    if numel(Archive) >= params.N
        Pop = Archive(1:params.N);
    else
        Pop = [Archive; Pop(1:params.N-numel(Archive))];
    end
    toc;
end

%% ========================== MOPO 特有辅助函数 ==========================

function Archive = AdaptiveEllipticalSelection(Archive, K)
    % 实现自适应椭圆分割筛选逻辑 [cite: 505, 533]
    Objs = vertcat(Archive.Objectives);
    nD = size(Objs, 1);
    
    % 计算个体间距离以确定椭圆轴长 [cite: 506, 512]
    dist_matrix = pdist2(Objs, Objs);
    dist_matrix(logical(eye(nD))) = inf;
    
    is_removed = false(nD, 1);
    while sum(~is_removed) > K
        remaining_indices = find(~is_removed);
        % 选取密度最高（最近距离最小）的个体
        [min_dist, idx_in_rem] = min(min(dist_matrix(remaining_indices, remaining_indices)));
        actual_idx = remaining_indices(idx_in_rem);
        
        % 以该点为中心设定椭圆长轴 a = min_dist [cite: 537]
        % 简化实现：剔除该点周围半径内的其他解 [cite: 507]
        neighbors = remaining_indices(dist_matrix(actual_idx, remaining_indices) < min_dist);
        if ~isempty(neighbors)
            is_removed(neighbors(1)) = true; % 剔除冗余解
        else
            is_removed(actual_idx) = true;
        end
    end
    Archive = Archive(~is_removed);
end

function L = Levy(dim)
    % Levy 飞行步长生成 [cite: 284, 345]
    beta = 1.5;
    sigma = (gamma(1+beta) * sin(pi*beta/2) / (gamma((1+beta)/2) * beta * 2^((beta-1)/2)))^(1/beta);
    u = randn(1, dim) * sigma;
    v = randn(1, dim);
    L = u ./ abs(v).^(1/beta);
end

% ... 位置向量转换函数 (CombinePos/SeparatePos) ...
%% ========================== MOPO 核心辅助函数 ==========================

function MeanPos = CalculateMeanPosition(Pop)
    % 计算种群中所有个体的平均位置，用于交流行为更新
    N = numel(Pop);
    
    % 提取并累加所有个体的组成部分
    sum_dep = 0; sum_bw = 0; sum_off = 0;
    for i = 1:N
        sum_dep = sum_dep + Pop(i).Position.deployment;
        sum_bw = sum_bw + Pop(i).Position.bandwidth;
        sum_off = sum_off + Pop(i).Position.offloading;
    end
    
    % 计算均值
    MeanPos.deployment = sum_dep / N;
    MeanPos.bandwidth = sum_bw / N;
    MeanPos.offloading = sum_off / N; % 离散变量在行为计算中暂按连续处理
end

function pos_vector = CombinePos(PositionStruct)
    % 将结构体形式的位置转换为扁平的行向量，以便进行算术运算（如 Levy 飞行）
    pos_vector = [PositionStruct.deployment(:)', ...
                  PositionStruct.bandwidth(:)', ...
                  PositionStruct.offloading(:)'];
end

function NewPositionStruct = SeparatePos(pos_vector, problem)
    % 将更新后的行向量重新拆分回结构体，并执行边界约束检查
    
    % 1. 计算各部分长度
    len_dep = numel(problem.bounds.deployment(1,:));
    len_bw = numel(problem.bounds.bandwidth(1,:));
    len_off = problem.nTerminals;
    
    % 2. 截取数据
    dep = pos_vector(1 : len_dep);
    bw = pos_vector(len_dep + 1 : len_dep + len_bw);
    off = pos_vector(len_dep + len_bw + 1 : end);
    
    % 3. 执行边界检查与约束处理
    % 连续变量：部署位置和带宽
    NewPositionStruct.deployment = max(min(dep, problem.bounds.deployment(2,:)), problem.bounds.deployment(1,:));
    NewPositionStruct.bandwidth = max(min(bw, problem.bounds.bandwidth(2,:)), problem.bounds.bandwidth(1,:));
    
    % 离散变量：卸载决策 (四舍五入并限制在雾节点索引范围内)
    off_rounded = round(off);
    NewPositionStruct.offloading = max(min(off_rounded, problem.nFogNodes), 1);
end

function off = DiscreteMutate(off, nFogNodes, prob)
    % 标准离散重置变异
    mask = rand(size(off)) < prob;
    off(mask) = randi([1, nFogNodes], 1, sum(mask));
end

function off = DiscreteCrossover(off1, off2)
    % 标准均匀交叉
    mask = rand(size(off1)) < 0.5;
    off = off1;
    off(mask) = off2(mask);
end