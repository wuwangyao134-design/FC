% 这是 NSGA-II 内部也用到的非支配排序机制
function non_dominated_indices = FindNonDominatedSolutions(objectives_matrix)
    if isempty(objectives_matrix)
        non_dominated_indices = [];
        return;
    end
    
    num_solutions = size(objectives_matrix, 1);
    is_dominated = false(num_solutions, 1);
    
    for i = 1:num_solutions
        if any(isinf(objectives_matrix(i,:))) || any(isnan(objectives_matrix(i,:)))
            is_dominated(i) = true; % 无效解被视为被支配
            continue; % 跳过与此无效解的比较
        end
        for j = 1:num_solutions
            if i == j, continue; end
            
            % 检查 j 是否支配 i
            % j 支配 i 的条件：j 在所有目标上都不劣于 i，且至少在一个目标上优于 i
            dominates = true;
            is_strictly_better = false;
            for k = 1:size(objectives_matrix, 2)
                if objectives_matrix(j, k) > objectives_matrix(i, k) % 最小化问题，j 的目标值更大表示更差
                    dominates = false;
                    break;
                end
                if objectives_matrix(j, k) < objectives_matrix(i, k) % j 至少在一个目标上优于 i
                    is_strictly_better = true;
                end
            end
            
            if dominates && is_strictly_better
                is_dominated(i) = true;
                break; % i 被支配，无需再检查其他解
            end
        end
    end
    non_dominated_indices = find(~is_dominated);
end