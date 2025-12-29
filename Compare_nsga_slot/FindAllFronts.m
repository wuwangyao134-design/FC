function Fronts = FindAllFronts(pop)
    Fronts = {};

    nPop_original = numel(pop);
    % 使用逻辑标志数组而非索引，更高效地追踪未分配的个体
    current_pop_flags = true(nPop_original, 1); 

    front_idx = 1;
    % 只要还有未分配的个体，就继续寻找下一层前沿
    while any(current_pop_flags)
        % 使用逻辑索引直接获取当前剩余种群
        current_remaining_pop_indices = find(current_pop_flags);
        current_remaining_pop = pop(current_remaining_pop_indices);

        % 寻找当前剩余种群中的第一非支配前沿
        current_front_members = FindNonDominated(current_remaining_pop);

        % 如果找不到新的非支配前沿，则停止算法
        if isempty(current_front_members)
            break; 
        end

        % 将找到的前沿添加到结果中
        Fronts{front_idx} = current_front_members;

        % 高效地标记当前前沿的成员为"已处理"
        % 提取当前前沿所有成员的目标值，用于快速预筛选
        current_front_objectives = vertcat(current_front_members.Objectives);

        % 遍历当前剩余种群的个体
        for i = 1:numel(current_remaining_pop)
            is_in_current_front = false;

            % 先用目标值进行快速预筛选
            for j = 1:size(current_front_objectives, 1)
                % 如果目标值相同，才进行更详细的比较
                if all(current_remaining_pop(i).Objectives == current_front_objectives(j,:))
                    % 目标值相同时，进一步比较Position以确保是同一个解
                    if isequal(current_remaining_pop(i).Position.deployment, current_front_members(j).Position.deployment) && ...
                       isequal(current_remaining_pop(i).Position.bandwidth, current_front_members(j).Position.bandwidth) && ...
                       isequal(current_remaining_pop(i).Position.offloading, current_front_members(j).Position.offloading)
                        is_in_current_front = true;
                        break; % 找到匹配，跳出内层循环
                    end
                end
            end

            if is_in_current_front
                % 找到该个体在原始种群中的索引，并将其标志设为 false
                original_index_of_this_solution = current_remaining_pop_indices(i);
                current_pop_flags(original_index_of_this_solution) = false;
            end
        end

        front_idx = front_idx + 1;
    end
end

