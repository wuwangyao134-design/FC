function [Pop, agent, info] = REMODDPG_Baseline(problem, params, agent)
%REMODDPG_BASELINE Preference-conditioned multi-objective DDPG baseline.
%
%   [POP, AGENT, INFO] = REMODDPG_BASELINE(PROBLEM, PARAMS, AGENT)
%
%   This implementation treats the optimization performed in one time slot
%   as a terminal one-step decision process.  A two-output critic predicts
%   the vector reward [delay reward, energy reward], while the actor is updated
%   with the preference-weighted deterministic policy gradient.  Because a
%   transition terminates after one decision, gamma = 0 and the Bellman
%   target is exactly the observed vector reward.  A replay buffer is still
%   used to decorrelate samples and stabilize training.
%
%   The continuous actor output is decoded into all three decision groups:
%     1) bandwidth allocation;
%     2) MTD-to-FN association (continuous rank parameter -> discrete FN);
%     3) two-dimensional FN deployment over the complete service region.
%   A deterministic capacity repair is applied only when the actor-proposed
%   association would overload an FN.  Hence offloading is no longer fixed
%   according to the initial FN deployment.
%
%   Passing AGENT=[] trains a new policy and then performs inference.
%   Passing a previously returned AGENT performs inference only.  This makes
%   it possible to report offline training time separately from online
%   per-slot inference time.
%
%   To prevent critic extrapolation from collapsing the deterministic actor
%   at discrete decision boundaries, training retains a small feasible elite
%   bank for uniformly spaced preference regions.  Inference evaluates a
%   fixed-size pool composed of actor outputs, nearby elite experiences, and
%   bounded local perturbations; the reported inference time includes all of
%   these objective evaluations.
%
%   修改记录：
%   [REMODDPG修改-01 | 2026-09-22]
%   原因：S2中出现确定性策略输出集中、非支配前沿规模过小和个别轮次
%   无有效前沿的问题。统一加入偏好分箱精英经验与推理候选细化机制，
%   同一组参数用于S1和S2，不进行逐场景调参。

    if nargin < 3
        agent = [];
    end

    cfg = make_config(problem, params);
    state_dim  = state_dimension(problem);
    action_dim = 2 * problem.nTerminals + 1 + 2 * problem.nFogNodes;
    context_state = build_state(problem, [0.5, 0.5]);

    must_train = isempty(agent);
    if ~must_train
        must_train = ~isfield(agent, 'state_dim') || ...
                     agent.state_dim ~= state_dim || ...
                     agent.action_dim ~= action_dim || ...
                     ~isfield(agent, 'actor_mode') || ...
                     ~strcmp(agent.actor_mode, 'dual_preference_head_v2');
        % [REMODDPG修改-01] 版本号升级，确保旧版塌缩策略不会被继续复用。
        if ~must_train
            must_train = ~isfield(agent, 'context_state') || ...
                numel(agent.context_state) ~= numel(context_state) || ...
                any(abs(agent.context_state - context_state) > 1e-12);
        end
    end

    info = struct('WasTrained', false, ...
                  'TrainingTime', 0, ...
                  'InferenceTime', 0, ...
                  'RewardHistory', [], ...
                  'FeasibleHistory', [], ...
                  'CriticLossHistory', [], ...
                  'ObjectiveScales', [], ...
                  'TrainingEvaluations', 0, ...
                  'InferenceEvaluations', 0, ...
                  'InferenceFeasibleCandidates', 0, ...
                  'InferenceSelectedFeasible', 0);
    % [REMODDPG修改-01] 增加候选池与最终输出的可行性诊断字段，便于定位塌缩。

    if must_train
        train_timer = tic;
        [agent, train_info] = train_agent(problem, cfg, state_dim, action_dim);
        info.WasTrained = true;
        info.TrainingTime = toc(train_timer);
        info.RewardHistory = train_info.RewardHistory;
        info.FeasibleHistory = train_info.FeasibleHistory;
        info.CriticLossHistory = train_info.CriticLossHistory;
        info.ObjectiveScales = agent.objective_scales;
        info.TrainingEvaluations = train_info.TrainingEvaluations;
    else
        info.ObjectiveScales = agent.objective_scales;
    end

    infer_timer = tic;
    [Pop, n_eval, inference_stats] = infer_population( ...
        problem, params.N, agent);
    info.InferenceTime = toc(infer_timer);
    info.InferenceEvaluations = n_eval;
    info.InferenceFeasibleCandidates = ...
        inference_stats.FeasibleCandidates;
    info.InferenceSelectedFeasible = ...
        inference_stats.SelectedFeasible;
    % [REMODDPG修改-01] 推理耗时包含候选生成、真实评价和筛选的全部时间。
end

function cfg = make_config(problem, params)
    cfg.hidden = get_option(params, 'drl_hidden', [128, 64]);
    cfg.train_episodes = get_option(params, 'drl_train_episodes', 3000);
    cfg.batch_size = get_option(params, 'drl_batch_size', 64);
    cfg.buffer_capacity = get_option(params, 'drl_buffer_capacity', ...
                                     max(5000, cfg.train_episodes));
    cfg.warmup_steps = get_option(params, 'drl_warmup_steps', ...
                                  min(512, max(128, round(0.1 * cfg.train_episodes))));
    cfg.scale_samples = get_option(params, 'drl_scale_samples', 128);
    cfg.actor_lr = get_option(params, 'drl_actor_lr', 1e-4);
    cfg.critic_lr = get_option(params, 'drl_critic_lr', 3e-4);
    cfg.noise_initial = get_option(params, 'drl_noise_initial', 0.30);
    cfg.noise_final = get_option(params, 'drl_noise_final', 0.03);
    cfg.penalty_weight = get_option(params, 'drl_penalty_weight', 10.0);
    cfg.policy_delay = get_option(params, 'drl_policy_delay', 2);
    cfg.gradient_clip = get_option(params, 'drl_gradient_clip', 5.0);
    cfg.reward_clip = get_option(params, 'drl_reward_clip', 50.0);
    cfg.preference_bins = get_option(params, 'drl_preference_bins', 21);
    cfg.elites_per_bin = get_option(params, 'drl_elites_per_bin', 4);
    cfg.inference_candidates_per_preference = get_option( ...
        params, 'drl_inference_candidates_per_preference', 8);
    cfg.inference_noise = get_option( ...
        params, 'drl_inference_noise', [0.02, 0.05, 0.10]);
    % [REMODDPG修改-01] 固定偏好覆盖、精英容量和推理扰动，不针对S1/S2分别调参。
    cfg.beta1 = 0.9;
    cfg.beta2 = 0.999;
    cfg.adam_eps = 1e-8;
    cfg.problem_signature = [problem.nTerminals, problem.nFogNodes, ...
                             problem.area(:)', problem.systemTotalBandwidth];
end

function value = get_option(params, name, default_value)
    if isfield(params, name) && ~isempty(params.(name))
        value = params.(name);
    else
        value = default_value;
    end
end

function n = state_dimension(problem)
    % Seven features per MTD, three per FN, four global features, and two
    % preference weights.
    n = 7 * problem.nTerminals + 3 * problem.nFogNodes + 6;
end

function [agent, info] = train_agent(problem, cfg, state_dim, action_dim)
    % Two explicit policy heads are learned: the first is anchored to the
    % delay objective and the second to the energy objective.  Their raw
    % actions are blended by the preference vector.  This prevents the
    % preference input from being ignored by a single shared output head.
    actor = init_mlp([state_dim, cfg.hidden, 2 * action_dim], true);
    critic = init_mlp([state_dim + action_dim, cfg.hidden, 2], false);
    actor_opt = init_adam(actor);
    critic_opt = init_adam(critic);

    objective_scales = estimate_objective_scales(problem, cfg, action_dim);
    buffer = init_replay(cfg.buffer_capacity, state_dim, action_dim);
    elite_archive = init_elite_archive( ...
        cfg.preference_bins, cfg.elites_per_bin, action_dim);
    % [REMODDPG修改-01] 保存各偏好区域中经真实环境评价的可行精英动作。

    reward_history = nan(cfg.train_episodes, 1);
    feasible_history = false(cfg.train_episodes, 1);
    critic_loss_history = nan(cfg.train_episodes, 1);

    for step = 1:cfg.train_episodes
        % Stratified preference sampling guarantees that both endpoints and
        % the interior of the Pareto preference interval are visited evenly.
        preference_index = mod(step - 1, cfg.preference_bins);
        w1 = preference_index / max(cfg.preference_bins - 1, 1);
        % [REMODDPG修改-01] 用统一分箱数覆盖两个目标之间的完整偏好区间。
        preference = [w1, 1 - w1];
        state = build_state(problem, preference);

        if step <= cfg.warmup_steps
            raw_action = randn(1, action_dim);
        else
            raw_action = dual_head_action(actor, state, action_dim);
            progress = (step - cfg.warmup_steps) / ...
                       max(1, cfg.train_episodes - cfg.warmup_steps);
            sigma = cfg.noise_initial * (1 - progress) + ...
                    cfg.noise_final * progress;
            raw_action = raw_action + sigma * randn(1, action_dim);
        end

        raw_action = max(min(raw_action, 8), -8);
        solution = decode_action(raw_action, problem);
        [reward_vec, feasible] = environment_reward( ...
            solution, problem, objective_scales, cfg);

        buffer = replay_add(buffer, state, raw_action, reward_vec);
        reward_history(step) = preference * reward_vec(:);
        feasible_history(step) = feasible;

        % Retain several genuinely evaluated feasible actions for every
        % preference region.  They keep the deterministic actor inside the
        % support of successful experience when the critic is inaccurate
        % around discrete association boundaries.
        if feasible
            scalar_reward = preference * reward_vec(:);
            elite_archive = elite_add(elite_archive, ...
                preference_index + 1, raw_action, reward_vec, scalar_reward);
        end
        % [REMODDPG修改-01] 仅收集真实可行动作，防止critic外推误差主导推理。

        if buffer.count >= cfg.batch_size
            [states, actions, rewards] = replay_sample(buffer, cfg.batch_size);

            % Terminal one-step transition: y = r because gamma = 0.
            [q_pred, critic_cache] = mlp_forward(critic, [states, actions]);
            critic_loss_history(step) = mean((q_pred(:) - rewards(:)).^2);
            d_q = 2 * (q_pred - rewards) / (cfg.batch_size * 2);
            [critic_grad, ~] = mlp_backward(critic, critic_cache, d_q);
            [critic, critic_opt] = adam_step( ...
                critic, critic_grad, critic_opt, cfg.critic_lr, cfg);

            if mod(step, cfg.policy_delay) == 0
                [actor_heads, actor_cache] = mlp_forward(actor, states);
                preferences = states(:, end-1:end);
                delay_head = actor_heads(:, 1:action_dim);
                energy_head = actor_heads(:, action_dim+1:end);
                actor_actions = preferences(:, 1) .* delay_head + ...
                                preferences(:, 2) .* energy_head;
                unclipped_actions = actor_actions;
                actor_actions = max(min(actor_actions, 8), -8);
                [~, policy_critic_cache] = mlp_forward( ...
                    critic, [states, actor_actions]);

                % L_actor = -mean(w_delay*Q_delay + w_energy*Q_energy).
                d_q_actor = -preferences / cfg.batch_size;
                [~, d_critic_input] = mlp_backward( ...
                    critic, policy_critic_cache, d_q_actor);
                d_actor_output = d_critic_input(:, state_dim+1:end);
                d_actor_output = d_actor_output .* (abs(unclipped_actions) < 8);

                % Back-propagate through the explicit convex mixture of the
                % delay and energy action heads.
                d_actor_heads = zeros(size(actor_heads));
                d_actor_heads(:, 1:action_dim) = ...
                    preferences(:, 1) .* d_actor_output;
                d_actor_heads(:, action_dim+1:end) = ...
                    preferences(:, 2) .* d_actor_output;
                [actor_grad, ~] = mlp_backward( ...
                    actor, actor_cache, d_actor_heads);
                [actor, actor_opt] = adam_step( ...
                    actor, actor_grad, actor_opt, cfg.actor_lr, cfg);
            end
        end
    end

    agent.actor = actor;
    agent.objective_scales = objective_scales;
    agent.state_dim = state_dim;
    agent.action_dim = action_dim;
    agent.actor_mode = 'dual_preference_head_v2';
    agent.elite_archive = elite_archive;
    % [REMODDPG修改-01] 将精英经验随agent保存，供后续时隙直接推理使用。
    agent.problem_signature = cfg.problem_signature;
    agent.context_state = build_state(problem, [0.5, 0.5]);
    agent.config = cfg;

    info.RewardHistory = reward_history;
    info.FeasibleHistory = feasible_history;
    info.CriticLossHistory = critic_loss_history;
    info.TrainingEvaluations = cfg.scale_samples + cfg.train_episodes;
end

function scales = estimate_objective_scales(problem, cfg, action_dim)
    samples = nan(cfg.scale_samples, 2);
    for k = 1:cfg.scale_samples
        solution = decode_action(randn(1, action_dim), problem);
        result = feval(problem.objFunc, solution, problem);
        if ~isfield(result, 'RawObjectives')
            error(['REMODDPG_Baseline requires EvaluateParticle to return ', ...
                   'RawObjectives, ConstraintViolation, and IsFeasible.']);
        end
        samples(k, :) = result.RawObjectives(1, :);
    end

    scales = nan(1, 2);
    for j = 1:2
        values = samples(:, j);
        values = values(isfinite(values) & values > 0);
        if isempty(values)
            scales(j) = 1;
        else
            scales(j) = median(values);
        end
    end
    scales = max(scales, [problem.Tslot * 0.05, 1e-6]);
end

function [reward_vec, feasible] = environment_reward(solution, problem, scales, cfg)
    result = feval(problem.objFunc, solution, problem);
    raw_obj = result.RawObjectives(1, :);
    violation = result.ConstraintViolation(1);
    feasible = logical(result.IsFeasible(1));

    normalized_cost = raw_obj ./ scales;
    reward_vec = -normalized_cost - cfg.penalty_weight * violation;
    reward_vec(~isfinite(reward_vec)) = -cfg.reward_clip;
    reward_vec = max(reward_vec, -cfg.reward_clip);
end

function [Pop, n_eval, stats] = infer_population( ...
        problem, population_size, agent)
    % [REMODDPG修改-01] 以actor为中心，结合邻近偏好精英和局部扰动构建候选池。
    template = struct('Position', [], 'Objectives', [], 'Tmax', [], ...
                      'Rank', [], 'CrowdingDistance', [], ...
                      'RawObjectives', [], 'ConstraintViolation', [], ...
                      'IsFeasible', []);
    Pop = repmat(template, population_size, 1);

    if population_size == 1
        weights = 0.5;
    else
        weights = linspace(0, 1, population_size);
    end

    cfg = agent.config;
    candidates_per_preference = max(1, ...
        round(cfg.inference_candidates_per_preference));
    candidate_actions = zeros( ...
        population_size * candidates_per_preference, agent.action_dim);
    next_candidate = 1;

    for p = 1:population_size
        preference = [weights(p), 1 - weights(p)];
        state = build_state(problem, preference);
        base_action = dual_head_action( ...
            agent.actor, state, agent.action_dim);
        local_actions = inference_actions( ...
            base_action, weights(p), agent, candidates_per_preference);
        last_candidate = next_candidate + size(local_actions, 1) - 1;
        candidate_actions(next_candidate:last_candidate, :) = local_actions;
        next_candidate = last_candidate + 1;
    end
    candidate_actions = candidate_actions(1:next_candidate-1, :);

    % Do not spend evaluations on exact duplicates introduced when a nearby
    % preference selects the same elite experience.
    rounded_actions = round(candidate_actions * 1e8) / 1e8;
    [~, unique_indices] = unique(rounded_actions, 'rows', 'stable');
    candidate_actions = candidate_actions(sort(unique_indices), :);

    n_candidates = size(candidate_actions, 1);
    candidate_population = repmat(struct('Position', []), n_candidates, 1);
    for c = 1:n_candidates
        candidate_population(c).Position = decode_action( ...
            candidate_actions(c, :), problem);
    end
    result = feval(problem.objFunc, candidate_population, problem);

    normalized_cost = result.RawObjectives ./ agent.objective_scales;
    normalized_cost(~isfinite(normalized_cost)) = cfg.reward_clip;
    violation = result.ConstraintViolation(:);
    violation(~isfinite(violation)) = cfg.reward_clip;
    feasible = logical(result.IsFeasible(:));
    used = false(n_candidates, 1);

    % Select one different candidate for each requested preference.  Actual
    % objective and constraint values, rather than critic predictions, are
    % used for this inexpensive policy-guided refinement.
    for p = 1:population_size
        preference = [weights(p), 1 - weights(p)];
        score = normalized_cost * preference(:) + ...
                cfg.penalty_weight * violation;
        score(~feasible) = score(~feasible) + cfg.reward_clip;
        score(used) = inf;
        [~, selected] = min(score);
        if isempty(selected) || ~isfinite(score(selected))
            score = normalized_cost * preference(:) + ...
                    cfg.penalty_weight * violation;
            [~, selected] = min(score);
        end
        used(selected) = true;

        Pop(p).Position = candidate_population(selected).Position;
        Pop(p).Objectives = result.Objectives(selected, :);
        Pop(p).Tmax = result.Tmax(selected);
        Pop(p).RawObjectives = result.RawObjectives(selected, :);
        Pop(p).ConstraintViolation = result.ConstraintViolation(selected);
        Pop(p).IsFeasible = result.IsFeasible(selected);
    end
    n_eval = n_candidates;
    stats.FeasibleCandidates = sum(feasible);
    stats.SelectedFeasible = sum([Pop.IsFeasible]);
    % [REMODDPG修改-01] 返回真实候选评价次数和可行数量，避免隐藏推理开销。
end

function actions = inference_actions( ...
        base_action, weight, agent, requested_count)
    % [REMODDPG修改-01] 为每个偏好生成统一数量的actor/精英/扰动候选动作。
    actions = base_action;
    seeds = base_action;

    if isfield(agent, 'elite_archive') && ...
            any(agent.elite_archive.valid(:))
        archive = agent.elite_archive;
        center_bin = 1 + round(weight * (archive.n_bins - 1));
        neighboring_bins = unique(max(1, min(archive.n_bins, ...
            [center_bin, center_bin - 1, center_bin + 1])));
        for b = neighboring_bins
            valid_slots = find(archive.valid(b, :));
            if isempty(valid_slots)
                continue;
            end
            [~, best_local] = max(archive.scores(b, valid_slots));
            elite_slot = valid_slots(best_local);
            elite_action = reshape( ...
                archive.actions(b, elite_slot, :), 1, []);
            seeds(end + 1, :) = elite_action; %#ok<AGROW>
            if size(actions, 1) < requested_count
                actions(end + 1, :) = elite_action; %#ok<AGROW>
            end
        end
    end

    noise_levels = agent.config.inference_noise(:)';
    if isempty(noise_levels)
        noise_levels = 0.05;
    end
    perturbation_index = 1;
    while size(actions, 1) < requested_count
        source_index = mod(perturbation_index - 1, size(seeds, 1)) + 1;
        noise_index = mod(perturbation_index - 1, numel(noise_levels)) + 1;
        candidate = seeds(source_index, :) + ...
            noise_levels(noise_index) * randn(size(base_action));
        actions(end + 1, :) = max(min(candidate, 8), -8); %#ok<AGROW>
        perturbation_index = perturbation_index + 1;
    end
    actions = actions(1:requested_count, :);
end

function action = dual_head_action(actor, state, action_dim)
    heads = mlp_predict(actor, state);
    preference = state(:, end-1:end);
    delay_head = heads(:, 1:action_dim);
    energy_head = heads(:, action_dim+1:end);
    action = preference(:, 1) .* delay_head + ...
             preference(:, 2) .* energy_head;
    action = max(min(action, 8), -8);
end

function state = build_state(problem, preference)
    n_term = problem.nTerminals;
    n_fog = problem.nFogNodes;
    area_x0 = problem.area(1, 1);
    area_x1 = problem.area(1, 2);
    area_y0 = problem.area(2, 1);
    area_y1 = problem.area(2, 2);
    range_x = max(area_x1 - area_x0, eps);
    range_y = max(area_y1 - area_y0, eps);

    term_pos = problem.terminalProperties.positions;
    task_size = problem.terminalProperties.task_sizes(:);
    pt_dbm = problem.terminalProperties.Pt_dbm(:);
    fc = problem.terminalProperties.fc(:);
    shadow_los = problem.fixed_shadow_LoS_val(:);
    shadow_nlos = problem.fixed_shadow_NLoS_val(:);

    term_features = [ ...
        (term_pos(:, 1) - area_x0) / range_x, ...
        (term_pos(:, 2) - area_y0) / range_y, ...
        task_size / 1e6, ...
        pt_dbm / 30, ...
        fc / 6e9, ...
        tanh(shadow_los / 10), ...
        tanh(shadow_nlos / 10)];

    fog_pos = problem.initial_fog_positions_matrix;
    cpu = problem.fogNodeProperties.cpu_cycle_rate(:);
    fog_features = [ ...
        (fog_pos(:, 1) - area_x0) / range_x, ...
        (fog_pos(:, 2) - area_y0) / range_y, ...
        cpu / 5e9];

    area_size = range_x * range_y;
    global_features = [n_term / 300, n_fog / 15, ...
                       n_term / max(area_size, 1) * 1e3, ...
                       problem.systemTotalBandwidth / max(n_term * 10e6, 1)];

    state = [reshape(term_features', 1, []), ...
             reshape(fog_features', 1, []), ...
             global_features, preference];
end

function solution = decode_action(raw_action, problem)
    n_term = problem.nTerminals;
    n_fog = problem.nFogNodes;

    idx = 0;
    bw_logits = raw_action(idx + (1:n_term));
    idx = idx + n_term;
    total_bw_raw = raw_action(idx + 1);
    idx = idx + 1;
    association_raw = raw_action(idx + (1:n_term));
    idx = idx + n_term;
    deployment_raw = raw_action(idx + (1:(2 * n_fog)));

    % Deployment spans the complete rectangular service region.
    deployment_unit = stable_sigmoid(deployment_raw);
    fog_positions = zeros(n_fog, 2);
    fog_positions(:, 1) = problem.area(1, 1) + ...
        deployment_unit(1:2:end)' * diff(problem.area(1, :));
    fog_positions(:, 2) = problem.area(2, 1) + ...
        deployment_unit(2:2:end)' * diff(problem.area(2, :));

    % Each continuous association parameter selects a distance-ranked FN.
    % The selected association is subsequently repaired only if required by
    % the per-FN CPU-capacity constraint.
    term_positions = problem.terminalProperties.positions;
    dx = term_positions(:, 1) - fog_positions(:, 1)';
    dy = term_positions(:, 2) - fog_positions(:, 2)';
    distances = sqrt(dx.^2 + dy.^2);
    [~, distance_order] = sort(distances, 2, 'ascend');
    rank_index = min(n_fog, floor(stable_sigmoid(association_raw(:)) * n_fog) + 1);
    preferred = distance_order(sub2ind([n_term, n_fog], ...
                                       (1:n_term)', rank_index));
    offloading = capacity_repair(preferred, distances, problem);

    % The extra scalar controls total bandwidth utilization; the remaining
    % logits distribute that total under all lower/upper bounds.
    lower = problem.bounds.bandwidth(1, :);
    upper = problem.bounds.bandwidth(2, :);
    min_total = sum(lower);
    max_total = min(problem.systemTotalBandwidth, sum(upper));
    if min_total > max_total + 1e-9
        error('Infeasible bandwidth bounds: sum(lower) exceeds total budget.');
    end
    target_total = min_total + stable_sigmoid(total_bw_raw) * ...
                   (max_total - min_total);
    bandwidth = bounded_softmax_allocation( ...
        bw_logits, lower, upper, target_total);

    solution.offloading = offloading;
    solution.bandwidth = bandwidth;
    solution.deployment = reshape(fog_positions', 1, []);
end

function offloading = capacity_repair(preferred, distances, problem)
    n_term = problem.nTerminals;
    task_cycles = problem.terminalProperties.task_sizes(:) * 1000;
    capacity = problem.fogNodeProperties.cpu_cycle_rate(:) * problem.Tslot;
    remaining = capacity;
    offloading = zeros(1, n_term);

    % Largest-first packing avoids creating an avoidable late infeasibility.
    [~, task_order] = sort(task_cycles, 'descend');
    diagonal = hypot(diff(problem.area(1, :)), diff(problem.area(2, :)));
    diagonal = max(diagonal, 1);

    for k = 1:n_term
        i = task_order(k);
        candidates = find(remaining >= task_cycles(i) - 1e-9);
        if isempty(candidates)
            % The instance itself may be CPU-infeasible. Select the FN with
            % the largest residual capacity and expose the violation to the
            % reward instead of silently discarding the solution.
            [~, chosen] = max(remaining);
        else
            score = distances(i, candidates) / diagonal + ...
                    0.35 * (candidates' ~= preferred(i));
            [~, local_index] = min(score);
            chosen = candidates(local_index);
        end
        offloading(i) = chosen;
        remaining(chosen) = remaining(chosen) - task_cycles(i);
    end
end

function bandwidth = bounded_softmax_allocation(logits, lower, upper, target_total)
    logits = logits - max(logits);
    weights = exp(max(min(logits, 40), -40));
    weights = max(weights, eps);

    bandwidth = lower;
    headroom = upper - lower;
    remaining = max(0, target_total - sum(lower));
    active = headroom > 1e-12;

    for iteration = 1:(numel(lower) + 1)
        if remaining <= 1e-6 || ~any(active)
            break;
        end
        active_weights = weights(active);
        shares = remaining * active_weights / sum(active_weights);
        active_indices = find(active);
        additions = min(shares, headroom(active));
        bandwidth(active_indices) = bandwidth(active_indices) + additions;
        headroom(active_indices) = headroom(active_indices) - additions;
        used = sum(additions);
        remaining = max(0, remaining - used);
        active = headroom > 1e-6;
        if used <= 1e-12
            break;
        end
    end
end

function y = stable_sigmoid(x)
    x = max(min(x, 40), -40);
    y = 1 ./ (1 + exp(-x));
end

function net = init_mlp(layers, small_output)
    n_layers = numel(layers) - 1;
    net.W = cell(n_layers, 1);
    net.b = cell(n_layers, 1);
    for layer = 1:n_layers
        if layer == n_layers && small_output
            net.W{layer} = 1e-3 * randn(layers(layer), layers(layer + 1));
        else
            net.W{layer} = randn(layers(layer), layers(layer + 1)) * ...
                           sqrt(2 / layers(layer));
        end
        net.b{layer} = zeros(1, layers(layer + 1));
    end
end

function y = mlp_predict(net, x)
    y = x;
    for layer = 1:numel(net.W)-1
        y = max(0, y * net.W{layer} + net.b{layer});
    end
    y = y * net.W{end} + net.b{end};
end

function [y, cache] = mlp_forward(net, x)
    n_layers = numel(net.W);
    cache.a = cell(n_layers + 1, 1);
    cache.z = cell(n_layers, 1);
    cache.a{1} = x;
    for layer = 1:n_layers
        cache.z{layer} = cache.a{layer} * net.W{layer} + net.b{layer};
        if layer < n_layers
            cache.a{layer + 1} = max(0, cache.z{layer});
        else
            cache.a{layer + 1} = cache.z{layer};
        end
    end
    y = cache.a{end};
end

function [grad, d_input] = mlp_backward(net, cache, d_output)
    n_layers = numel(net.W);
    grad.W = cell(n_layers, 1);
    grad.b = cell(n_layers, 1);
    delta = d_output;

    for layer = n_layers:-1:1
        grad.W{layer} = cache.a{layer}' * delta;
        grad.b{layer} = sum(delta, 1);
        d_previous = delta * net.W{layer}';
        if layer > 1
            delta = d_previous .* (cache.z{layer - 1} > 0);
        end
    end
    d_input = d_previous;
end

function opt = init_adam(net)
    opt.mW = cell(size(net.W));
    opt.vW = cell(size(net.W));
    opt.mb = cell(size(net.b));
    opt.vb = cell(size(net.b));
    for layer = 1:numel(net.W)
        opt.mW{layer} = zeros(size(net.W{layer}));
        opt.vW{layer} = zeros(size(net.W{layer}));
        opt.mb{layer} = zeros(size(net.b{layer}));
        opt.vb{layer} = zeros(size(net.b{layer}));
    end
    opt.t = 0;
end

function [net, opt] = adam_step(net, grad, opt, learning_rate, cfg)
    opt.t = opt.t + 1;
    for layer = 1:numel(net.W)
        gW = max(min(grad.W{layer}, cfg.gradient_clip), -cfg.gradient_clip);
        gb = max(min(grad.b{layer}, cfg.gradient_clip), -cfg.gradient_clip);

        opt.mW{layer} = cfg.beta1 * opt.mW{layer} + (1-cfg.beta1) * gW;
        opt.vW{layer} = cfg.beta2 * opt.vW{layer} + (1-cfg.beta2) * (gW.^2);
        opt.mb{layer} = cfg.beta1 * opt.mb{layer} + (1-cfg.beta1) * gb;
        opt.vb{layer} = cfg.beta2 * opt.vb{layer} + (1-cfg.beta2) * (gb.^2);

        mW_hat = opt.mW{layer} / (1 - cfg.beta1^opt.t);
        vW_hat = opt.vW{layer} / (1 - cfg.beta2^opt.t);
        mb_hat = opt.mb{layer} / (1 - cfg.beta1^opt.t);
        vb_hat = opt.vb{layer} / (1 - cfg.beta2^opt.t);

        net.W{layer} = net.W{layer} - learning_rate * ...
            mW_hat ./ (sqrt(vW_hat) + cfg.adam_eps);
        net.b{layer} = net.b{layer} - learning_rate * ...
            mb_hat ./ (sqrt(vb_hat) + cfg.adam_eps);
    end
end

function buffer = init_replay(capacity, state_dim, action_dim)
    buffer.states = zeros(capacity, state_dim);
    buffer.actions = zeros(capacity, action_dim);
    buffer.rewards = zeros(capacity, 2);
    buffer.capacity = capacity;
    buffer.count = 0;
    buffer.next = 1;
end

function archive = init_elite_archive(n_bins, elites_per_bin, action_dim)
    % [REMODDPG修改-01] 初始化分偏好的可行精英经验库。
    archive.n_bins = n_bins;
    archive.elites_per_bin = elites_per_bin;
    archive.actions = nan(n_bins, elites_per_bin, action_dim);
    archive.reward_vectors = nan(n_bins, elites_per_bin, 2);
    archive.scores = -inf(n_bins, elites_per_bin);
    archive.valid = false(n_bins, elites_per_bin);
end

function archive = elite_add( ...
        archive, bin_index, action, reward_vector, scalar_reward)
    % [REMODDPG修改-01] 每个偏好分箱仅保留标量回报最高的固定数量动作。
    invalid_slot = find(~archive.valid(bin_index, :), 1);
    if ~isempty(invalid_slot)
        target_slot = invalid_slot;
    else
        [worst_score, target_slot] = min(archive.scores(bin_index, :));
        if scalar_reward <= worst_score
            return;
        end
    end

    archive.actions(bin_index, target_slot, :) = reshape(action, 1, 1, []);
    archive.reward_vectors(bin_index, target_slot, :) = ...
        reshape(reward_vector, 1, 1, []);
    archive.scores(bin_index, target_slot) = scalar_reward;
    archive.valid(bin_index, target_slot) = true;
end

function buffer = replay_add(buffer, state, action, reward)
    index = buffer.next;
    buffer.states(index, :) = state;
    buffer.actions(index, :) = action;
    buffer.rewards(index, :) = reward;
    buffer.count = min(buffer.count + 1, buffer.capacity);
    buffer.next = mod(index, buffer.capacity) + 1;
end

function [states, actions, rewards] = replay_sample(buffer, batch_size)
    indices = randperm(buffer.count, batch_size);
    states = buffer.states(indices, :);
    actions = buffer.actions(indices, :);
    rewards = buffer.rewards(indices, :);
end
