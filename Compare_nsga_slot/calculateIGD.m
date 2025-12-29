%% **新增辅助函数：计算反向世代距离 (Inverted Generational Distance, IGD)**
% 输入: obtained_front (N x M 矩阵), reference_front (P x M 矩阵，真实或近似前沿)
function igd = calculateIGD(obtained_front, reference_front)
    if isempty(obtained_front) || isempty(reference_front)
        igd = inf; % 如果任一为空，则 IGD 为无穷大
        return;
    end
    
    num_ref_points = size(reference_front, 1);
    distances = zeros(num_ref_points, 1);
    
    % 遍历参考前沿中的每个点
    for i = 1:num_ref_points
        ref_point = reference_front(i, :);
        min_dist_to_obtained = inf;
        
        % 计算该参考点到 obtained_front 中所有点的最小欧氏距离
        for j = 1:size(obtained_front, 1)
            dist = norm(ref_point - obtained_front(j, :));
            if dist < min_dist_to_obtained
                min_dist_to_obtained = dist;
            end
        end
        distances(i) = min_dist_to_obtained;
    end
    
    % 求平均距离
    igd = mean(distances);
end
