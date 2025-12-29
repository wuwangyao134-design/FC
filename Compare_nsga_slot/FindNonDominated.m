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
