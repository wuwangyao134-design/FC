function Fronts = FindAllFronts(pop)
    % 1. 预提取所有目标值，形成 N x M 矩阵 (M=2, G1和G2)
    objs = vertcat(pop.Objectives);
    nPop = size(objs, 1);
    
    % 2. 初始化
    Fronts = {};
    is_assigned = false(nPop, 1); % 追踪哪些个体已被分配到前沿
    front_idx = 1;
    
    % 3. 预计算支配关系矩阵 (这是提速的核心！)
    % DomMatrix(i,j) = 1 表示 i 支配 j
    % 逻辑：i的所有目标 <= j，且至少有一个目标 < j
    DomMatrix = false(nPop, nPop);
    for i = 1:nPop
        % 向量化比较：个体 i 与所有其他个体比较
        less_equal = all(objs(i,:) <= objs, 2);
        strictly_less = any(objs(i,:) < objs, 2);
        DomMatrix(i,:) = less_equal & strictly_less;
    end
    
    % 计算每个个体被支配的次数 (被多少人压着)
    domination_count = sum(DomMatrix, 1)';
    
    % 4. 剥洋葱式提取前沿
    while ~all(is_assigned)
        % 当前前沿是个体：尚未分配 且 被支配次数为0
        current_front_mask = (domination_count == 0) & ~is_assigned;
        
        if ~any(current_front_mask)
            break; 
        end
        
        % 存储当前层级的个体
        Fronts{front_idx} = pop(current_front_mask);
        
        % 标记为已分配
        is_assigned(current_front_mask) = true;
        
        % 关键：模拟“移除”当前前沿，更新剩余个体的被支配次数
        % 找到当前前沿支配的所有个体
        dominated_by_current = any(DomMatrix(current_front_mask, :), 1)';
        % 这些个体的被支配次数减去 1 (或根据具体支配关系减去对应次数)
        % 更稳健的做法：直接更新剩余个体的计数
        % 这里使用标准快速非支配排序逻辑的变体
        indices_in_front = find(current_front_mask);
        for idx = indices_in_front'
            targets = DomMatrix(idx, :);
            domination_count(targets) = domination_count(targets) - 1;
        end
        
        % 强制将已分配个体的计数设为无穷大，防止重复选中
        domination_count(current_front_mask) = inf;
        
        front_idx = front_idx + 1;
    end
end