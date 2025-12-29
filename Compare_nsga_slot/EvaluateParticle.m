% =========================================================================
% 核心评估函数 (EvaluateParticle.m) - 优化版 (使用固定于每次run_idx的阴影衰落)
% =========================================================================
function Results = EvaluateParticle(position, problem)
    %% 0. 定义通用物理和系统常数
    No_dbm_per_hz = -174;
    task_cycles_per_bit = 1000; % 每个比特数据在CPU上处理所需的周期数
    eta = 0.2;
    gamma = 5e-6;
    penalty_value = 1e9; % 惩罚值，用于处理约束违反

    % --- 移除：这里不再定义和生成随机衰落 ---
    % shadow_std_dev_LoS = 3.0;
    % shadow_std_dev_NLoS = 8.29;

    %% 1. 解析决策变量和异构属性
    nTerminals = problem.nTerminals;
    nFogNodes = problem.nFogNodes;

    % 解码决策变量
    fog_positions = reshape(position.deployment, [2, nFogNodes])'; % 雾节点部署位置 (K x 2 矩阵)
    bandwidths = position.bandwidth; % 分配给终端的带宽 (1 x M 向量)
    task_sizes_vec = problem.terminalProperties.task_sizes; % 终端任务大小 (1 x M 向量)
    offloading_plan = position.offloading; % 终端卸载方案 (1 x M 向量, 整数)
    offloading_plan = round(offloading_plan); % 确保索引是整数，防止浮点数微小偏差导致报错

    % 获取问题参数 (终端属性和雾节点CPU率)
    terminal_positions = problem.terminalProperties.positions; % 终端位置 (M x 2 矩阵)
    Pt_dbm_vec = problem.terminalProperties.Pt_dbm; % 终端发射功率 (1 x M 向量)
    fc_vec = problem.terminalProperties.fc; % 终端工作频率 (1 x M 向量)
    cpu_cycle_rate_vec = problem.fogNodeProperties.cpu_cycle_rate; % 雾节点CPU周期率 (1 x K 向量)
    
    % --- 新增：从 problem 结构体中获取固定的阴影衰落值 ---
    % 这些值将在 main_compare.m 的 run_idx 循环中生成并传入
    fixed_shadow_LoS_val = problem.fixed_shadow_LoS_val;
    fixed_shadow_NLoS_val = problem.fixed_shadow_NLoS_val;
    % -----------------------------------------------------

    % 获取约束参数
    Tslot = problem.Tslot; 
    systemTotalBandwidth = problem.systemTotalBandwidth;

    % 如果总带宽超出系统总带宽，按比例缩放 (确保可行)
    if sum(bandwidths) > systemTotalBandwidth
       bandwidths = bandwidths * (systemTotalBandwidth / sum(bandwidths));
    end

    %% 2. 计算通信延迟 (T2) - 优化：大部分计算已向量化
    % 关键：计算每个终端到其选择的雾节点的距离
    % 假设 offloading_plan 包含了终端选择的雾节点索引
    d = zeros(1, nTerminals);
    for i = 1:nTerminals
        d(i) = norm(terminal_positions(i, :) - fog_positions(offloading_plan(i), :));
    end
    d = max(d, 1); % 避免距离过小导致计算问题 (d < 1 改为 max(d,1))

    % 路径损耗模型 (向量化计算)
    d1 = 5; d2 = 65;
    pLos = (d < d1) + (d >= d1) .* exp(-(d - d1) / d2);

    % --- 关键修改：在这里使用从 problem 传入的固定阴影衰落值 ---
    PL_LoS_db = 16.9.*log10(d) + 32.8 + 20.*log10(fc_vec/1e9) + fixed_shadow_LoS_val;
    PL_NLoS_db = 38.3.*log10(d) + 17.30 + 24.9.*log10(fc_vec/1e9) + fixed_shadow_NLoS_val;
    % ----------------------------------------------------------

    PL_db = pLos .* PL_LoS_db + (1 - pLos) .* PL_NLoS_db; % 注意这里是 .* 元素乘法 (PL_NLoS_db 是向量，pLos是向量，因此这里应该是 .*)

    % 信噪比 (SNR) 和信道容量 (Cn) (向量化计算)
    No_dbm = No_dbm_per_hz + 10 * log10(bandwidths);
    SNR_db = Pt_dbm_vec - PL_db - No_dbm;
    SNR_linear = 10.^(SNR_db / 10); % 元素级别 10 的次方
    Cn = bandwidths .* log2(1 + SNR_linear); % 元素级别乘法

    % 通信延迟 (T2) (向量化计算)
    T2 = task_sizes_vec ./ Cn; % 元素级别除法
    T2(Cn <= 0 | isnan(Cn) | isinf(Cn)) = inf; % 检查信道容量是否有效

    %% 3. 计算处理与排队延迟 (部分优化：内部循环逻辑不变)
    % 这部分涉及到任务排队顺序，不易完全向量化
    T_finish = zeros(1, nTerminals); % 记录每个任务的完成时间
    T_max = 0; % 记录所有任务中最晚的完成时间

    for k = 1:nFogNodes % 遍历每个雾节点
        terminals_at_k_indices = find(offloading_plan == k); % 找到分配给雾节点 k 的所有终端
        if ~isempty(terminals_at_k_indices) % 如果有任务分配到此雾节点
            arrivals_at_k = T2(terminals_at_k_indices); % 获取这些任务的到达时间 (通信延迟)

            % 将终端索引和到达时间配对，以便排序后能找到原始终端
            [sorted_arrivals, sort_order_relative] = sort(arrivals_at_k);
            sorted_indices_global = terminals_at_k_indices(sort_order_relative); % 得到排序后的原始终端索引

            cpu_cycle_rate = cpu_cycle_rate_vec(k); % 获取当前雾节点的CPU周期率
            last_finish_time_at_k = 0; % 记录该雾节点上上一个任务的完成时间

            % 计算处理时间 for current_terminal_idx
            T3_current_terminals_at_k = (task_sizes_vec(sorted_indices_global) * task_cycles_per_bit) / cpu_cycle_rate; % 向量化计算所有T3

            for j = 1:numel(sorted_indices_global) % 遍历分配到此雾节点的每个任务
                current_terminal_idx = sorted_indices_global(j);
                T3_this_terminal = T3_current_terminals_at_k(j); % 获取当前任务的处理时间
                arrival_time = sorted_arrivals(j); % 获取当前任务的到达时间 (已排序)

                start_time = max(arrival_time, last_finish_time_at_k); % 计算任务开始处理时间
                finish_time = start_time + T3_this_terminal; % 计算任务完成时间

                T_finish(current_terminal_idx) = finish_time; % 更新原始终端的完成时间
                last_finish_time_at_k = finish_time; % 更新该雾节点上一个任务的完成时间
            end
            if last_finish_time_at_k > T_max, T_max = last_finish_time_at_k; end % 更新全局最晚完成时间
        end
    end

    %% 4. 计算原始目标函数 (在施加惩罚之前)
    G1 = mean(T_finish); % 平均总延迟

    Pt_W_vec = 10.^(Pt_dbm_vec / 10) / 1e3;    
    E = (Pt_W_vec /eta + gamma.* bandwidths) .* T2;
    % G2 = sum(E); % 终端总能耗
    epsilon_m_val = 0.1e-6;
    E_comp = epsilon_m_val .* task_sizes_vec; 
    
    % 总能耗 G2 = 所有终端的 (传输能耗 + 计算能耗) 之和
    G2 = sum(E + E_comp);

    % === 新增：检查雾节点 CPU 容量约束 (保持不变，因为循环难以向量化) ===
    is_cpu_capacity_violated = false; % 初始化CPU容量违反标志
    for k = 1:nFogNodes % 遍历每个雾节点
        terminals_at_k_indices = find(offloading_plan == k); % 找到分配给当前雾节点 k 的所有终端索引
        if ~isempty(terminals_at_k_indices) % 如果有终端分配到这个雾节点
            total_cycles_demanded_at_k = sum(task_sizes_vec(terminals_at_k_indices) * task_cycles_per_bit); % 优化：sum 向量化
            max_cycles_per_slot_k = cpu_cycle_rate_vec(k) * Tslot; 
            if total_cycles_demanded_at_k > max_cycles_per_slot_k
                is_cpu_capacity_violated = true; % 只要有一个节点超载，就标记为违反
                % 移除 fprintf('超载'); 以减少运行时开销，如果需要调试，可以使用 disp 或 log
                break; 
            end
        end
    end

    %% 5. 处理约束 (惩罚法 - 约束违反度优化版)
    % 计算各项约束的具体违反程度 (Violation Degree)
    % 1. 时间约束违反度
    cv_time = max(0, (T_max - Tslot) / Tslot); 
    
    % 2. 雾节点 CPU 容量违反度
    cv_cpu = 0;
    for k = 1:nFogNodes
        terminals_at_k = find(offloading_plan == k);
        if ~isempty(terminals_at_k)
            total_demand = sum(task_sizes_vec(terminals_at_k) * task_cycles_per_bit);
            max_capacity = cpu_cycle_rate_vec(k) * Tslot;
            if total_demand > max_capacity
                cv_cpu = cv_cpu + (total_demand - max_capacity) / max_capacity;
            end
        end
    end
    
    % 3. 总带宽约束违反度
    cv_bw = max(0, (sum(bandwidths) - systemTotalBandwidth) / systemTotalBandwidth);
    
    % 总违反度
    total_cv = cv_time + cv_cpu + cv_bw;

    % 施加惩罚
    if total_cv > 0
        % 基础惩罚项 penalty_value (1e9) 确保不可行解永远劣于可行解
        % 附加惩罚项 (total_cv * 1e6) 确保不可行解之间也能“比出好坏”
        G1 = G1 + penalty_value + total_cv * 1e6;
        G2 = G2 + penalty_value + total_cv * 1e6;
    end
    %% 6. 返回结果
    Results.Objectives = [G1, G2]; % 最终的目标函数值 (可能包含惩罚)
    Results.Tmax = T_max; % 返回最晚任务完成时间 (用于调试和详细分析)
end