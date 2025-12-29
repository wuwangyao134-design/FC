% function hv = calculateHV(obtained_front, reference_point)
%     if isempty(obtained_front)
%         hv = 0;
%         return;
%     end
%     if size(obtained_front, 2) ~= length(reference_point)
%         error('目标维度与参考点维度不匹配。');
%     end
% 
%     % 在最小化问题中，HV通常计算的是参考点到解集之间的超体积。
%     % 对于2D问题，可以简化为计算矩形面积之和。
%     % 需要将解集相对于参考点进行翻转，以便使用现有的超体积计算方法（通常是最大化问题）
%     % 或者更直接地，对于每个解 (f1, f2)，其贡献的矩形面积是 (ref_f1 - f1) * (ref_f2 - f2)
% 
%     % 确保reference_point每个维度都大于或等于obtained_front的相应维度
%     % 如果不是，可能参考点设置不合理
%     if any(obtained_front > reference_point, 'all')
%         warning('参考点可能设置不合理，部分解超出参考点。HV计算可能不准确。');
%     end
% 
%     % 过滤掉目标值大于参考点的解，这些解不会贡献超体积（实际上是“惩罚”解）
%     % 并且 HV 值通常计算的是帕累托前沿到参考点之间的超体积
%     % 对于最小化问题，超体积是指被帕累托前沿“支配”且在参考点边界内的区域
% 
%     % Step 1: 对帕累托前沿按第一个目标进行排序
%     [~, sortIdx] = sort(obtained_front(:,1));
%     sorted_front = obtained_front(sortIdx, :);
% 
%     hv = 0;
%     num_objectives = size(sorted_front, 2);
% 
%     if num_objectives == 2 % 仅支持2D超体积计算
%         % 假设目标都是越小越好 (最小化问题)
%         % 从最右侧/最下方（最差）的解开始计算，逐步向最优解累加矩形面积
% 
%         % 添加一个虚拟点，其第一个目标为-无穷大，第二个目标为参考点第二个维度
%         % 这确保了最左侧的解也能形成一个完整的矩形区域
%         virtual_point_left = [-inf, reference_point(2)];
%         % 添加一个虚拟点，其第一个目标为参考点第一个维度，第二个目标为-无穷大
%         % 这确保了最下方的解也能形成一个完整的矩形区域
%         % 但对于2D的HV计算，通常可以这样处理:
%         % 将每个点 (f1, f2) 的贡献视为 (ref_f1 - f1) * (ref_f2 - f2)
%         % 并且需要从最右边的点开始计算。
% 
%         % HV的计算逻辑通常是：遍历排序后的解，每个解与其前一个解（或参考点）共同定义一个矩形。
%         % 在最小化问题中，我们通常将目标值反转，使其成为最大化问题，然后使用现成的HV算法，
%         % 或者采用专门为最小化设计的算法。
% 
%         % 一个简化的2D超体积计算（假定参考点在所有前沿点右上方）：
%         % 沿着f1轴（或f2轴）排序，然后计算矩形面积
% 
%         % 遍历排序后的帕累托前沿点
%         for i = 1:size(sorted_front, 1)
%             p_i = sorted_front(i,:);
% 
%             % 确保解在参考点内部，否则不贡献HV (或者贡献负值，取决于实现)
%             if any(p_i > reference_point)
%                 continue; % 跳过支配参考点的点，这些点不贡献超体积
%             end
% 
%             % 计算当前矩形的f1边长
%             rect_f1_length = reference_point(1) - p_i(1);
% 
%             % 计算当前矩形的f2边长
%             if i == size(sorted_front, 1) % 最右边的点 (f1最大)
%                 % 它的f2边长由参考点f2决定
%                 rect_f2_length = reference_point(2) - p_i(2);
%             else
%                 % 它的f2边长由下一个点（f1更小）的f2决定
%                 % 取下一个点的f2和当前点f2的最小值作为高度，以避免重叠
%                 rect_f2_length = reference_point(2) - max(p_i(2), sorted_front(i+1, 2));
%             end
% 
%             % 确保边长非负
%             rect_f1_length = max(0, rect_f1_length);
%             rect_f2_length = max(0, rect_f2_length);
% 
%             hv = hv + rect_f1_length * rect_f2_length;
%         end
%         % 以上HV计算方式可能过于简化，对更复杂的形状不准确。
%         % 推荐使用更健壮的HV库函数，例如 PlatEMO 工具箱中的 'Metric.HV' 或类似实现。
%         % 这里提供一个更通用的概念化实现，可能需要根据具体库函数进行调整。
% 
%         % 更准确的2D HV计算 (基于矩形分解)
%         hv = 0;
% 
%         % 1. 确保所有目标都小于参考点
%         % 过滤掉任何一个目标都大于参考点的解
%         filtered_for_hv = obtained_front(all(obtained_front <= reference_point, 2), :);
% 
%         if isempty(filtered_for_hv)
%             hv = 0;
%             return;
%         end
% 
%         % 2. 找出这些解中的非支配解（如果传入的already_front不是严格的非支配集）
%         non_dominated_idx = FindNonDominatedSolutions(filtered_for_hv);
%         pf = filtered_for_hv(non_dominated_idx, :);
% 
%         % 3. 按第一个目标排序（升序）
%         [~, sort_idx] = sort(pf(:,1));
%         pf = pf(sort_idx, :);
% 
%         % 4. 添加一个辅助点到前沿的左侧，用于计算第一个矩形
%         % 这个点的f1是负无穷，f2是参考点f2。或者，取f1的最小值，f2是参考点f2
%         % 为了方便计算，可以将目标值从原点平移到参考点为零点，然后求面积
% 
%         % 目标值相对于参考点的距离 (所有值都应该是正的)
%         shifted_front = reference_point - pf;
% 
%         % 按第一个目标 (f1) 降序排序
%         [~, sort_idx_shifted] = sort(shifted_front(:,1), 'descend');
%         shifted_front = shifted_front(sort_idx_shifted, :);
% 
%         % 计算矩形
%         hv = 0;
%         current_max_f2 = 0; % 追踪当前最大 f2 投影高度
%         for i = 1:size(shifted_front, 1)
%             f1_val = shifted_front(i,1);
%             f2_val = shifted_front(i,2);
% 
%             if f2_val > current_max_f2
%                 % 计算当前点与之前最大f2之间的矩形面积
%                 % 宽度由f1决定，高度由f2增量决定
%                 area_width = f1_val;
%                 area_height = f2_val - current_max_f2;
%                 hv = hv + area_width * area_height;
%                 current_max_f2 = f2_val;
%             end
%         end
% 
%     else % 对于更高维度的HV计算，需要更复杂的算法，如Zitzler的算法或专用的HV库
%         error('目前仅支持2D超体积 (HV) 计算。');
%     end
% end
function hv = calculateHV(obtained_front, reference_point)
    % 基本检查
    if isempty(obtained_front)
        hv = 0;
        return;
    end
    if size(obtained_front, 2) ~= length(reference_point)
        error('目标维度与参考点维度不匹配。');
    end
    
    % 检查参考点
    if any(obtained_front > reference_point, 'all')
        warning('参考点可能设置不合理，部分解超出参考点。HV计算可能不准确。');
    end

    % 选择计算方法（可以通过参数控制使用哪种方法）
    method = 3;  % 1: 直接计算法, 2: 坐标转换法, 3: 论文定义法
    
    switch method
        case 1
            % ================ 方法1：直接计算法 ================
            [~, sortIdx] = sort(obtained_front(:,1));
            sorted_front = obtained_front(sortIdx, :);
            
            hv = 0;
            for i = 1:size(sorted_front, 1)
                p_i = sorted_front(i,:);
                
                if any(p_i > reference_point)
                    continue;
                end
                
                rect_f1_length = reference_point(1) - p_i(1);
                
                if i == size(sorted_front, 1)
                    rect_f2_length = reference_point(2) - p_i(2);
                else
                    rect_f2_length = reference_point(2) - max(p_i(2), sorted_front(i+1, 2));
                end
                
                rect_f1_length = max(0, rect_f1_length);
                rect_f2_length = max(0, rect_f2_length);
                
                hv = hv + rect_f1_length * rect_f2_length;
            end
            
        case 2
            % ================ 方法2：坐标转换法 ================
            filtered_for_hv = obtained_front(all(obtained_front <= reference_point, 2), :);
            
            if isempty(filtered_for_hv)
                hv = 0;
                return;
            end

            non_dominated_idx = FindNonDominatedSolutions(filtered_for_hv);
            pf = filtered_for_hv(non_dominated_idx, :);
            
            shifted_front = reference_point - pf;
            [~, sort_idx_shifted] = sort(shifted_front(:,1), 'descend');
            shifted_front = shifted_front(sort_idx_shifted, :);

            hv = 0;
            current_max_f2 = 0;
            for i = 1:size(shifted_front, 1)
                f1_val = shifted_front(i,1);
                f2_val = shifted_front(i,2);
                
                if f2_val > current_max_f2
                    area_width = f1_val;
                    area_height = f2_val - current_max_f2;
                    hv = hv + area_width * area_height;
                    current_max_f2 = f2_val;
                end
            end
            
        case 3
            % ================ 方法3：论文定义法 ================
            % HV(PF,r) = Volume(∪_{p∈PF} hypercube(p,r))
            
            % 1. 过滤无效解并获取非支配解
            valid_solutions = obtained_front(all(obtained_front <= reference_point, 2), :);
            if isempty(valid_solutions)
                hv = 0;
                return;
            end
            
            non_dominated_idx = FindNonDominatedSolutions(valid_solutions);
            pf = valid_solutions(non_dominated_idx, :);
            
            % 2. 按第一个目标排序
            [~, sort_idx] = sort(pf(:,1));
            sorted_pf = pf(sort_idx, :);
            
            % 3. 计算超体积
            hv = 0;
            prev_f2 = reference_point(2);
            
            for i = 1:size(sorted_pf, 1)
                current_point = sorted_pf(i,:);
                
                % 计算当前点贡献的hypercube体积
                width = reference_point(1) - current_point(1);
                height = prev_f2 - current_point(2);
                
                if width > 0 && height > 0
                    hv = hv + width * height;
                end
                
                prev_f2 = current_point(2);
            end
    end
end