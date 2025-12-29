function spread = calculateSpread(obtained_front, reference_front)
    % 基本检查
    if isempty(obtained_front) || isempty(reference_front) || ...
       size(obtained_front,1) < 2 || size(reference_front,1) < 2
        spread = NaN;
        return;
    end

    % 1. 计算参考前沿的范围R
    R = sum(max(reference_front) - min(reference_front));
    if R < eps
        spread = NaN;
        return;
    end

    % 2. 计算边界偏差
    ext_obtained = [min(obtained_front); max(obtained_front)];
    ext_reference = [min(reference_front); max(reference_front)];
    delta_ext = sum(vecnorm(ext_obtained - ext_reference, 2, 2));

    % 3. 计算相邻解之间的距离
    [~, idx] = sort(obtained_front(:,1));
    sorted_front = obtained_front(idx,:);
    distances = zeros(size(sorted_front,1)-1, 1);

    for i = 1:size(sorted_front,1)-1
        distances(i) = norm(sorted_front(i+1,:) - sorted_front(i,:));
    end

    % 计算平均距离
    d_bar = mean(distances);

    % 4. 计算最终的spread值
    spread = (delta_ext + sum(abs(distances - d_bar))) / R;
end
