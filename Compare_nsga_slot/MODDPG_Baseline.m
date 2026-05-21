function Pop = MODDPG_Baseline(problem, params, ~)
    % =====================================================================
    % 标准多目标 DDPG (Multi-Objective DDPG) 基线
    % 优化变量：带宽分配 + 雾节点位置（连续动作）
    % 100% 纯 MATLAB 实现，无需任何工具箱，运行极速且稳定
    % =====================================================================
    particle_template = struct('Position', [], 'Objectives', [], 'Tmax', [], 'Rank', [], 'CrowdingDistance', []); 
    N = params.N;
    Pop = repmat(particle_template, N, 1);
    
    nTerm = problem.nTerminals;
    nFog  = problem.nFogNodes;
    action_dim_bw   = nTerm;
    action_dim_deploy = 2 * nFog;
    action_dim_total = action_dim_bw + action_dim_deploy;
    
    % ===================== 超参数（自适应、稳定） =====================
    max_episodes     = 120;
    lr_actor         = 0.0015;
    lr_critic        = 0.0025;
    tau              = 0.005;
    batch_size       = 32;
    gamma            = 0;
    action_noise     = 0.12;
    hidden_dim       = 64;
    
    % ===================== 网络结构与初始化 =====================
    % 4维输入: [任务均值, 0.5, w1, w2]
    actor_net  = init_mlp([4, hidden_dim, action_dim_total]);
    target_act = actor_net;
    
    % Critic 输入: 4维状态 + 动作维度，输出 2维 Q 值 (时延, 能耗)
    critic_net = init_mlp([4 + action_dim_total, hidden_dim, 2]);
    target_cri = critic_net;
    
    replay_buffer = cell(0, 4);
    
    % ===================== 环境固定信息 =====================
    term_pos     = problem.terminalProperties.positions;
    task_sizes   = problem.terminalProperties.task_sizes;
    max_task     = max(task_sizes);
    norm_tasks   = task_sizes / max_task;
    deploy_fixed = problem.initial_fog_deployment_flat;
    
    area_min     = problem.area(1,1);
    area_max     = problem.area(1,2);
    deploy_min   = ones(1, action_dim_deploy) * area_min;
    deploy_max   = ones(1, action_dim_deploy) * area_max;
    
    bw_min       = problem.bounds.bandwidth(1,:);
    bw_max       = problem.bounds.bandwidth(2,:);
    total_bw     = problem.systemTotalBandwidth;
    
    % 就近卸载策略（固定，保证公平）
    fog_init = reshape(deploy_fixed, [2, nFog])';
    offload_assign = zeros(1, nTerm);
    for i = 1:nTerm
        d = vecnorm(term_pos(i,:) - fog_init, 2, 2);
        [~, idx] = min(d);
        offload_assign(i) = idx;
    end

    % ===================== 【新增】奖励记录 =====================
    reward_history = [];

    %% ===================== 训练阶段 =====================
   % ===================== 超参数修正 =====================
    % 在外层将 gamma 改为 0 (单步模型绝对不能有未来折扣)
    gamma = 0.0; 

    %% ===================== 训练阶段 =====================
    for ep = 1:max_episodes
        w1 = rand(); w2 = 1 - w1;
        state = [mean(norm_tasks), 0.5, w1, w2];
        
        % 1. 生成动作并添加探索噪声 (这里逻辑没问题)
        action_raw = predict_fwd(actor_net, state);
        action_raw = action_raw + action_noise * randn(1, action_dim_total);
        
        % 2. 解耦动作：带宽 + 部署
        bw_raw  = action_raw(1:action_dim_bw);
        dep_raw = action_raw(action_dim_bw+1:end);
        
        bw_act  = 1 ./ (1 + exp(-bw_raw)); % Sigmoid
        dep_act = tanh(dep_raw);           % Tanh
        
        % 3. 边界投影与缩放
        deploy  = deploy_fixed + dep_act * (area_max - area_min) * 0.35;
        deploy  = max(deploy, deploy_min);
        deploy  = min(deploy, deploy_max);
        
        % 4. 保命防御：防止带宽除零
        sum_bw = sum(bw_act);
        if sum_bw < 1e-6 || isnan(sum_bw)
            bw_act = ones(1, nTerm) / nTerm;
            sum_bw = 1.0;
        end
        bw_plan = (bw_act / sum_bw) * total_bw;
        bw_plan = max(bw_plan, bw_min);
        bw_plan = min(bw_plan, bw_max);
        
        % 5. 评估
        sol_tmp.offloading = offload_assign;
        sol_tmp.bandwidth  = bw_plan;
        sol_tmp.deployment = deploy;
        
        res = feval(problem.objFunc, sol_tmp, problem);
        g1 = res.Objectives(1);
        g2 = res.Objectives(2);
        
        % 🌟 修复3：连续势能惩罚，引导网络走向可行域
        if res.Tmax > problem.Tslot || g1 >= 1e8
            r1 = -5.0 * (res.Tmax / problem.Tslot); 
            r2 = -5.0 * (res.Tmax / problem.Tslot);
        else
            r1 = -(g1 / 5.0);
            r2 = -(g2 / 1000.0); % 假设能耗量级在百/千级别，做适当缩放
        end
        
        % 🌟 修复1：奖励必须包含权重，否则学不出 Pareto 前沿
        total_reward = w1 * r1 + w2 * r2;
        reward_history = [reward_history, total_reward];
        
        % 🌟 修复2：存入经验池 (注意，这里的 r 直接存标量 total_reward，因为优化目标被融合成了一个)
        replay_buffer(end+1, :) = {state, action_raw, total_reward, state};
        if size(replay_buffer, 1) > 2000
            replay_buffer(1, :) = [];
        end
        
        % 7. 神经网络参数更新
        if size(replay_buffer, 1) >= batch_size
            idx = randperm(size(replay_buffer, 1), batch_size);
            batch = replay_buffer(idx, :);
            
            s_batch  = cat(1, batch{:,1});
            a_batch  = cat(1, batch{:,2});
            r_batch  = cat(1, batch{:,3});
            ns_batch = cat(1, batch{:,4}); % 这里的 ns_batch 其实没用了，因为 gamma = 0
            
            % --- Critic 更新 (gamma=0，所以 y_target 直接等于 r_batch) ---
            y_target = r_batch;
            critic_net = bp_update(critic_net, [s_batch, a_batch], y_target, lr_critic);
            
            % --- Actor 更新 ---
            a_pred = predict_fwd(actor_net, s_batch);
            actor_net = bp_actor_update(actor_net, s_batch, a_pred, critic_net, lr_actor);
            
            % --- 目标网络软更新 ---
            target_act = soft_update(target_act, actor_net, tau);
            target_cri = soft_update(target_cri, critic_net, tau);
        end
    end

    % ===================== 【新增】绘制奖励收敛曲线 =====================
    % figure('Color','white');
    % plot(reward_history, 'b-', 'LineWidth', 1.5);
    % xlabel('Training Episode', 'FontSize', 11);
    % ylabel('Total Reward', 'FontSize', 11);
    % title('MODDPG Training Convergence Curve', 'FontSize', 12);
    % grid on; grid minor;
    % saveas(gcf, 'MODDPG_Convergence.fig');
    % saveas(gcf, 'MODDPG_Convergence.png');

    %% ===================== 推理生成解集 =====================
    for p = 1:N
        w1 = (p-1)/(N-1);
        w2 = 1 - w1;
        state = [mean(norm_tasks), 0.5, w1, w2];
        
        action_raw = predict_fwd(target_act, state);
        
        bw_raw  = action_raw(1:action_dim_bw);
        dep_raw = action_raw(action_dim_bw+1:end);
        
        bw_act  = 1 ./ (1 + exp(-bw_raw));
        dep_act = tanh(dep_raw);
        
        deploy  = deploy_fixed + dep_act * (area_max - area_min) * 0.35;
        deploy  = max(deploy, deploy_min);
        deploy  = min(deploy, deploy_max);
        
        sum_bw  = sum(bw_act);
        if sum_bw < 1e-6
            bw_act = ones(1, nTerm);
            sum_bw = nTerm;
        end
        bw_plan = (bw_act / sum_bw) * total_bw;
        bw_plan = max(bw_plan, bw_min);
        bw_plan = min(bw_plan, bw_max);
        
        sol.offloading = offload_assign;
        sol.bandwidth  = bw_plan;
        sol.deployment = deploy;
        
        res = feval(problem.objFunc, sol, problem);
        Pop(p).Position = sol;
        Pop(p).Objectives = res.Objectives;
        Pop(p).Tmax = res.Tmax;
    end
end

%% ===================== 神经网络解析梯度反向传播工具 =====================
function net = init_mlp(layers)
    net.layers = layers;
    net.w = cell(length(layers)-1, 1);
    net.b = cell(length(layers)-1, 1);
    for i = 1:length(layers)-1
        net.w{i} = randn(layers(i), layers(i+1)) * sqrt(2 / layers(i));
        net.b{i} = zeros(1, layers(i+1));
    end
end

function out = predict_fwd(net, x)
    out = x;
    for i = 1:length(net.w)-1
        out = max(0, out * net.w{i} + net.b{i}); % ReLU
    end
    out = out * net.w{end} + net.b{end}; % Linear 输出
end

function net = bp_update(net, x, y, lr)
    % 批处理正向与反向传播 (用于 Critic 的均方误差最小化)
    z = cell(length(net.w), 1);
    a = cell(length(net.w)+1, 1);
    a{1} = x;
    
    % 正向传播
    for i = 1:length(net.w)-1
        z{i} = a{i} * net.w{i} + net.b{i};
        a{i+1} = max(0, z{i});
    end
    z{end} = a{end-1} * net.w{end} + net.b{end};
    a{end} = z{end};
    
    % 反向传播
    delta = a{end} - y;
    dW = cell(length(net.w), 1);
    db = cell(length(net.w), 1);
    
    dW{end} = a{end-1}' * delta;
    db{end} = sum(delta, 1);
    
    for i = length(net.w)-1:-1:1
        delta = (delta * net.w{i+1}') .* (a{i+1} > 0);
        dW{i} = a{i}' * delta;
        db{i} = sum(delta, 1);
    end
    
    % 梯度裁剪与参数更新
    clip_v = 5;
    for i = 1:length(net.w)
        net.w{i} = net.w{i} - lr * max(min(dW{i}/size(x,1), clip_v), -clip_v);
        net.b{i} = net.b{i} - lr * max(min(db{i}/size(x,1), clip_v), -clip_v);
    end
end

function net = bp_actor_update(net, x, a_pred, critic, lr)
    % 显式策略梯度更新 (提升 Q 值目标)
    % 链式法则求导: d_Loss/d_W = d_Critic/d_Act * d_Act/d_W
    
    % 1. 估算 Critic 关于 Action 的梯度
    delta_cri = 1e-4;
    grad_act = zeros(size(a_pred));
    for i = 1:size(a_pred, 2)
        a_plus = a_pred;
        a_plus(:, i) = a_plus(:, i) + delta_cri;
        q_plus = predict_fwd(critic, [x, a_plus]);
        
        a_minus = a_pred;
        a_minus(:, i) = a_minus(:, i) - delta_cri;
        q_minus = predict_fwd(critic, [x, a_minus]);
        
        % 优化两目标之和 Q1 + Q2
        grad_act(:, i) = mean((sum(q_plus, 2) - sum(q_minus, 2)) / (2 * delta_cri), 1);
    end
    
    % 2. 将梯度反向传播给 Actor
    z = cell(length(net.w), 1);
    a = cell(length(net.w)+1, 1);
    a{1} = x;
    
    for i = 1:length(net.w)-1
        z{i} = a{i} * net.w{i} + net.b{i};
        a{i+1} = max(0, z{i});
    end
    z{end} = a{end-1} * net.w{end} + net.b{end};
    
    % 负梯度，因为我们要最大化 Q 值
    delta = -grad_act; 
    
    dW = cell(length(net.w), 1);
    db = cell(length(net.w), 1);
    dW{end} = a{end-1}' * delta;
    db{end} = sum(delta, 1);
    
    for i = length(net.w)-1:-1:1
        delta = (delta * net.w{i+1}') .* (a{i+1} > 0);
        dW{i} = a{i}' * delta;
        db{i} = sum(delta, 1);
    end
    
    % 参数更新
    clip_v = 5;
    for i = 1:length(net.w)
        net.w{i} = net.w{i} - lr * max(min(dW{i}/size(x,1), clip_v), -clip_v);
        net.b{i} = net.b{i} - lr * max(min(db{i}/size(x,1), clip_v), -clip_v);
    end
end

function t = soft_update(t, s, tau)
    for i = 1:length(t.w)
        t.w{i} = (1-tau)*t.w{i} + tau*s.w{i};
        t.b{i} = (1-tau)*t.b{i} + tau*s.b{i};
    end
end