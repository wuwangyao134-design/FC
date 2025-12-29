function Results = EvaluateParticle(input_data, problem)
    % --- 智能识别输入类型 ---
    if isstruct(input_data) && isfield(input_data, 'Position')
        % 情况1: 输入是种群数组 (IMyNSGA_II 使用)
        pop_positions = input_data;
    elseif isstruct(input_data) && isfield(input_data, 'deployment')
        % 情况2: 输入是单个 Position 结构体 (旧版 MOPSO/NSGA_II 使用)
        % 包装成数组格式以兼容后续逻辑
        temp_pop.Position = input_data;
        pop_positions = temp_pop;
    else
        error('EvaluateParticle: 输入数据格式无法识别');
    end

    nPop = numel(pop_positions);
    
    nPop = numel(pop_positions);
    nTerminals = problem.nTerminals;
    nFogNodes = problem.nFogNodes;
    
    % --- 1. 预提取问题常量 (避免在循环内重复读取结构体) ---
    term_pos = problem.terminalProperties.positions; % M x 2
    task_sizes = problem.terminalProperties.task_sizes; % 1 x M
    Pt_dbm = problem.terminalProperties.Pt_dbm;
    fc = problem.terminalProperties.fc;
    cpu_rates = problem.fogNodeProperties.cpu_cycle_rate;
    Tslot = problem.Tslot;
    sysBW = problem.systemTotalBandwidth;
    fixed_shadow_LoS = problem.fixed_shadow_LoS_val;
    fixed_shadow_NLoS = problem.fixed_shadow_NLoS_val;
    
    % 预计算常量
    Pt_W = 10.^(Pt_dbm / 10) / 1e3;
    E_comp = 0.1e-6 .* task_sizes;
    task_cycles = task_sizes * 1000; % task_cycles_per_bit = 1000

    % 初始化结果矩阵
    all_Objectives = zeros(nPop, 2);
    all_Tmax = zeros(nPop, 1);

    % --- 2. 种群循环 (由于逻辑复杂，我们保留一层循环，但优化内部) ---
    for p = 1:nPop
        % 获取当前个体的 Position 结构体
        current_pos = pop_positions(p).Position; 
        
        % 解析 (现在层级对齐了)
        f_pos = reshape(current_pos.deployment, [2, nFogNodes])'; 
        bws = current_pos.bandwidth;
        off_plan = current_pos.offloading;
        % 带宽约束强制修正 (原地修改加速)
        sumBW = sum(bws);
        if sumBW > sysBW
            bws = bws * (sysBW / sumBW);
        end

        % --- 向量化距离计算 ---
        % 使用矩阵索引直接提取终端对应的雾节点位置，避免循环 norm
        selected_f_pos = f_pos(off_plan, :); 
        d = sqrt(sum((term_pos - selected_f_pos).^2, 2))'; % 1 x M 向量
        d = max(d, 1);

        % --- 快速路径损耗 ---
        pLos = (d < 5) + (d >= 5) .* exp(-(d - 5) / 65);
        PL_LoS = 16.9*log10(d) + 32.8 + 20*log10(fc/1e9) + fixed_shadow_LoS;
        PL_NLoS = 38.3*log10(d) + 17.3 + 24.9*log10(fc/1e9) + fixed_shadow_NLoS;
        PL = pLos .* PL_LoS + (1 - pLos) .* PL_NLoS;

        % --- 信噪比与容量 ---
        No_dbm = -174 + 10 * log10(bws);
        SNR_linear = 10.^((Pt_dbm - PL - No_dbm) / 10);
        Cn = bws .* log2(1 + SNR_linear);
        T2 = task_sizes ./ Cn;
        T2(Cn <= 0 | isnan(Cn)) = 1e6; % 替代 inf 方便计算

        % --- 优化排队逻辑 ---
        T_finish = zeros(1, nTerminals);
        node_T_max = 0;
        is_cpu_violated = false;
        
        for k = 1:nFogNodes
            idx_k = (off_plan == k);
            if any(idx_k)
                % CPU 容量检查
                if sum(task_cycles(idx_k)) > cpu_rates(k) * Tslot
                    is_cpu_violated = true;
                end
                
                % 排序与排队
                [sorted_arrivals, s_order] = sort(T2(idx_k));
                current_T3 = task_cycles(idx_k) / cpu_rates(k);
                sorted_T3 = current_T3(s_order);
                
                % 使用更快的累加逻辑替代内部循环计算完成时间
                finish_times = zeros(1, numel(sorted_arrivals));
                last_f = 0;
                for j = 1:numel(sorted_arrivals)
                    last_f = max(sorted_arrivals(j), last_f) + sorted_T3(j);
                    finish_times(j) = last_f;
                end
                
                % 写回 T_finish (逻辑索引写回)
                temp_finish = T_finish(idx_k);
                temp_finish(s_order) = finish_times;
                T_finish(idx_k) = temp_finish;
                
                if last_f > node_T_max, node_T_max = last_f; end
            end
        end

        % --- 惩罚项与目标 ---
        G1 = mean(T_finish);
        E_trans = (Pt_W / 0.2 + 5e-6 .* bws) .* T2;
        G2 = sum(E_trans + E_comp);

        if node_T_max > Tslot || is_cpu_violated
            G1 = G1 + 1e9;
            G2 = G2 + 1e9;
        end
        
        all_Objectives(p, :) = [G1, G2];
        all_Tmax(p) = node_T_max;
    end
    
    Results.Objectives = all_Objectives;
    Results.Tmax = all_Tmax;
end