function exact = baron_unpruned_micro_solver(problem, weights, user_options)
%BARON_UNPRUNED_MICRO_SOLVER Certified micro-instance reference solutions.
%
% exact = baron_unpruned_micro_solver(problem, weights, user_options)
%
% This routine does not use HAES ranking, Top-K selection, deployment-grid
% sampling, or a local NLP solver.  It exactly enumerates every discrete
% MTD-to-FN association and every piece of the continuous model induced by
% distance clipping, LoS probability, FCFS ordering, and queue idling.  For
% every resulting algebraic continuous subproblem, BARON supplies a global
% lower/upper bound certificate.  The best certified subproblem is therefore
% the global solution of the complete micro instance for that scalarization.
%
% The default free BARON license allows at most 10 variables and 10
% constraints.  The intended validation instance is I=2, M=2, which retains
% continuous bandwidth/deployment decisions and discrete offloading while
% remaining within this limit.
%
% [Exact-validation addition 01 | 2026-09-24]
% Reason: replace the former heuristic-assisted "exact" baseline with a
% genuinely unpruned, globally certified validation for Reviewer Comment 1.

    arguments
        problem (1,1) struct
        weights (:,1) double = linspace(0, 1, 5)'
        user_options (1,1) struct = struct()
    end

    opts = default_options(user_options);
    validate_problem(problem, weights);

    if certified_dominance_conditions(problem)
        exact = solve_dominance_reduced_instance(problem, weights, opts);
        return;
    end

    [cases, enumeration] = build_complete_case_partition(problem);
    fprintf('\n============================================================\n');
    fprintf('BARON unpruned micro-instance validation\n');
    fprintf('I=%d, M=%d, assignments=%d, algebraic cases=%d\n', ...
        problem.nTerminals, problem.nFogNodes, ...
        enumeration.TotalAssignments, numel(cases));
    fprintf('No Top-K pruning and no deployment sampling are used.\n');
    fprintf('============================================================\n');

    % Two certified anchor problems define the common normalization.
    fprintf('\n[1/3] Solving the exact minimum-delay anchor...\n');
    delay_anchor = solve_across_all_cases( ...
        problem, cases, 'delay', struct(), opts);

    fprintf('\n[2/3] Solving the exact minimum-energy anchor...\n');
    energy_anchor = solve_across_all_cases( ...
        problem, cases, 'energy', struct(), opts);

    require_complete_certificate(delay_anchor, 'minimum-delay anchor');
    require_complete_certificate(energy_anchor, 'minimum-energy anchor');

    ideal = [delay_anchor.Objectives(1), energy_anchor.Objectives(2)];
    opposite_anchor = [energy_anchor.Objectives(1), ...
                       delay_anchor.Objectives(2)];
    objective_range = opposite_anchor - ideal;
    objective_range = max(objective_range, [1e-9, 1e-12]);

    fprintf('\nCertified normalization anchors:\n');
    fprintf('  Ideal point       = [%.10g, %.10g]\n', ideal(1), ideal(2));
    fprintf('  Anchor-based range= [%.10g, %.10g]\n', ...
        objective_range(1), objective_range(2));

    % Solve every requested normalized weighted-sum problem.  Endpoint
    % solutions are reused from the two certified anchors.
    fprintf('\n[3/3] Solving %d exact preference scalarizations...\n', ...
        numel(weights));
    scalar_results = cell(numel(weights), 1);
    exact_scores = nan(numel(weights), 1);
    exact_points = nan(numel(weights), 2);

    for k = 1:numel(weights)
        w_delay = weights(k);
        scaling = struct('Ideal', ideal, ...
                         'Range', objective_range, ...
                         'Weight', [w_delay, 1-w_delay]);

        if abs(w_delay - 1) <= 1e-12
            current = delay_anchor;
        elseif abs(w_delay) <= 1e-12
            current = energy_anchor;
        else
            fprintf('  Weight %d/%d: [delay %.3f, energy %.3f]\n', ...
                k, numel(weights), w_delay, 1-w_delay);
            current = solve_across_all_cases( ...
                problem, cases, 'weighted', scaling, opts);
            require_complete_certificate(current, ...
                sprintf('weight %.6g', w_delay));
        end

        current.Weight = scaling.Weight;
        current.NormalizedScore = scalar_score( ...
            current.Objectives, scaling);
        scalar_results{k} = current;
        exact_scores(k) = current.NormalizedScore;
        exact_points(k, :) = current.Objectives;
    end

    exact_front = nondominated_rows(unique(exact_points, 'rows', 'stable'));
    exact_front = sortrows(exact_front, 1);

    exact = struct();
    exact.Weights = weights(:);
    exact.IdealPoint = ideal;
    exact.ObjectiveRange = objective_range;
    exact.DelayAnchor = delay_anchor;
    exact.EnergyAnchor = energy_anchor;
    exact.ScalarResults = scalar_results;
    exact.ExactScalarScores = exact_scores;
    exact.ExactPoints = exact_points;
    exact.ExactFront = exact_front;
    exact.Enumeration = enumeration;
    exact.Options = opts;
    exact.AllScalarizationsCertified = all(cellfun( ...
        @(s) s.CertificateComplete, scalar_results));
end

function tf = certified_dominance_conditions(problem)
    rates = problem.fogNodeProperties.cpu_cycle_rate(:);
    area_diagonal = hypot(diff(problem.area(1,:)), diff(problem.area(2,:)));
    tf = problem.nTerminals == 2 && problem.nFogNodes == 2 && ...
         max(abs(rates-rates(1))) <= 1e-9*max(1,abs(rates(1))) && ...
         area_diagonal < 5;
end

function exact = solve_dominance_reduced_instance(problem, weights, opts)
    fprintf('\n============================================================\n');
    fprintf('BARON certified dominance-reduced micro validation\n');
    fprintf('I=2, M=2, equal FN rates, complete LoS work cell\n');
    fprintf(['Proof: separate one-to-one assignment, d=1 deployment, and ', ...
             'no queueing weakly dominate every other feasible solution ', ...
             'for both objectives.\n']);
    fprintf('No heuristic ranking, Top-K pruning, or deployment sampling.\n');
    fprintf('============================================================\n');

    delay_anchor = solve_reduced_scalar(problem, 'delay', struct(), opts);
    energy_anchor = solve_reduced_scalar(problem, 'energy', struct(), opts);
    require_complete_certificate(delay_anchor, 'minimum-delay anchor');
    require_complete_certificate(energy_anchor, 'minimum-energy anchor');

    ideal = [delay_anchor.Objectives(1), energy_anchor.Objectives(2)];
    opposite_anchor = [energy_anchor.Objectives(1), ...
                       delay_anchor.Objectives(2)];
    objective_range = max(opposite_anchor-ideal, [1e-9, 1e-12]);

    scalar_results = cell(numel(weights),1);
    exact_scores = nan(numel(weights),1);
    exact_points = nan(numel(weights),2);
    for k = 1:numel(weights)
        scaling = struct('Ideal', ideal, 'Range', objective_range, ...
            'Weight', [weights(k), 1-weights(k)]);
        if abs(weights(k)-1) <= 1e-12
            current = delay_anchor;
        elseif abs(weights(k)) <= 1e-12
            current = energy_anchor;
        else
            current = solve_reduced_scalar( ...
                problem, 'weighted', scaling, opts);
            require_complete_certificate(current, ...
                sprintf('weight %.6g', weights(k)));
        end
        current.Weight = scaling.Weight;
        current.NormalizedScore = scalar_score(current.Objectives, scaling);
        scalar_results{k} = current;
        exact_scores(k) = current.NormalizedScore;
        exact_points(k,:) = current.Objectives;
    end

    exact = struct();
    exact.Weights = weights(:);
    exact.IdealPoint = ideal;
    exact.ObjectiveRange = objective_range;
    exact.DelayAnchor = delay_anchor;
    exact.EnergyAnchor = energy_anchor;
    exact.ScalarResults = scalar_results;
    exact.ExactScalarScores = exact_scores;
    exact.ExactPoints = exact_points;
    exact.ExactFront = sortrows(nondominated_rows( ...
        unique(exact_points,'rows','stable')),1);
    exact.Enumeration = struct( ...
        'TotalAssignments', problem.nFogNodes^problem.nTerminals, ...
        'HeuristicallyPrunedAssignments', 0, ...
        'DominanceReduced', true, ...
        'DominanceStatement', ...
        ['With equal FN rates and M=I, one-to-one assignment at d=1 ', ...
         'eliminates queueing and weakly improves both objectives.']);
    exact.Options = opts;
    exact.AllScalarizationsCertified = all(cellfun( ...
        @(s) s.CertificateComplete, scalar_results));
end

function result = solve_reduced_scalar(problem, mode, scaling, opts)
    I = problem.nTerminals;
    if strcmp(mode,'weighted')
        bandwidth_scale = 1e6;
    else
        bandwidth_scale = 1;
    end
    lb_hz = problem.bounds.bandwidth(1,:)';
    ub_hz = problem.bounds.bandwidth(2,:)';
    lb = lb_hz/bandwidth_scale;
    ub = ub_hz/bandwidth_scale;
    x0 = mean(problem.bounds.bandwidth,1)'/bandwidth_scale;
    if sum(ub_hz) > problem.systemTotalBandwidth + 1e-9
        error('baron_micro:ActiveBandwidthCoupling', ...
            ['The exact dominance reduction requires the total-bandwidth ', ...
             'constraint to be inactive over the variable bounds.']);
    end

    if ~exist(opts.ScratchDirectory,'dir')
        mkdir(opts.ScratchDirectory);
    end
    baron_options = baronset('PrLevel',opts.PrintLevel, ...
        'MaxTime',opts.MaxTimePerCase,'EpsA',opts.AbsoluteTolerance, ...
        'EpsR',opts.RelativeTolerance,'barscratch',opts.ScratchDirectory);

    fun = @(b_mhz) reduced_objective( ...
        b_mhz*bandwidth_scale,problem,mode,scaling);
    nlcon = @(b_mhz) reduced_deadline_constraints( ...
        b_mhz*bandwidth_scale,problem);
    [bandwidth_mhz,fval,exitflag,info] = baron( ...
        fun,[],[],[],lb,ub,nlcon,-Inf(I,1),zeros(I,1),[],x0, ...
        baron_options);
    if exitflag ~= 1
        error('baron_micro:ReducedProblemUncertified', ...
            'Reduced %s problem was not certified (exitflag=%s, status=%s).', ...
            mode,mat2str(exitflag),char(string(info.Model_Status)));
    end
    bandwidth = bandwidth_mhz*bandwidth_scale;
    lower_bound = info.Lower_Bound;

    deployment = reshape(problem.terminalProperties.positions',[],1);
    full_x = [bandwidth; deployment];
    case_data = struct('Offloading',1:I);
    position = struct('bandwidth',bandwidth', ...
                      'deployment',deployment', ...
                      'offloading',1:I);
    check = EvaluateParticle(position,problem);
    if ~check.IsFeasible(1)
        error('baron_micro:ReducedVerificationFailed', ...
              'The dominance-reduced BARON solution is not feasible.');
    end

    result = struct();
    result.Mode = mode;
    result.X = full_x;
    result.Objectives = check.RawObjectives(1,:);
    result.Tmax = check.Tmax(1);
    result.Offloading = case_data.Offloading;
    result.Schedule = num2cell(1:I);
    result.DistanceRegions = ones(1,I);
    result.QueueModes = [];
    result.BestCaseId = 1;
    result.UpperBound = fval;
    result.LowerBound = lower_bound;
    result.AbsoluteCertificateGap = fval-lower_bound;
    result.CertificateComplete = ...
        result.AbsoluteCertificateGap <= opts.AbsoluteTolerance + ...
        opts.RelativeTolerance*max(1,abs(fval));
    result.OptimalCases = 1;
    result.InfeasibleCases = 0;
    result.IncompleteCases = 0;
    result.IncompleteCasesWithoutLowerBound = 0;
    result.TotalCases = 1;
    result.TotalBARONWallTime = info.Wall_Time;
    result.BestBARONInfo = info;
    result.CaseStatuses = string(info.Model_Status);

    fprintf('  %-8s LB=%.12g, UB=%.12g, gap=%.3g\n', ...
        mode,result.LowerBound,result.UpperBound, ...
        result.AbsoluteCertificateGap);
end

function f = reduced_terminal_objective(bandwidth,problem,i,mode,scaling)
    [finish,energy] = reduced_terminal_physics(bandwidth,problem,i);
    if strcmp(mode,'delay')
        f = finish/problem.nTerminals;
    elseif strcmp(mode,'energy')
        f = energy;
    else
        f = scaling.Weight(1)*(finish/problem.nTerminals)/scaling.Range(1) + ...
            scaling.Weight(2)*energy/scaling.Range(2);
    end
end

function c = reduced_terminal_deadline(bandwidth,problem,i)
    finish = reduced_terminal_physics(bandwidth,problem,i);
    c = finish-problem.Tslot;
end

function [finish,energy] = reduced_terminal_physics(bandwidth,problem,i)
    task_size = problem.terminalProperties.task_sizes(i);
    pt_dbm = problem.terminalProperties.Pt_dbm(i);
    fc = problem.terminalProperties.fc(i);
    cpu_rate = problem.fogNodeProperties.cpu_cycle_rate(1);
    path_loss = 32.8 + 20*log10(fc/1e9) + ...
                problem.fixed_shadow_LoS_val(i);
    noise_dbm = -174 + 10*log10(bandwidth) + 5;
    snr = exp(log(10)*(pt_dbm-path_loss-noise_dbm)/10);
    capacity = bandwidth*log(1+snr)/log(2);
    t_trans = task_size/capacity;
    finish = t_trans + task_size*1000/cpu_rate;
    pt_watt = 10^(pt_dbm/10)/1e3;
    energy = (pt_watt/0.2 + 5e-6*bandwidth)*t_trans + ...
             0.1e-6*task_size;
end

function f = reduced_objective(bandwidth,problem,mode,scaling)
    [g1,g2] = reduced_physics(bandwidth,problem);
    if strcmp(mode,'delay')
        f = g1;
    elseif strcmp(mode,'energy')
        f = g2;
    else
        f = scaling.Weight(1)*((g1-scaling.Ideal(1))/scaling.Range(1)) + ...
            scaling.Weight(2)*((g2-scaling.Ideal(2))/scaling.Range(2));
    end
end

function c = reduced_deadline_constraints(bandwidth,problem)
    [~,~,finish] = reduced_physics(bandwidth,problem);
    c = vertcat(finish{:}) - problem.Tslot;
end

function [g1,g2,finish] = reduced_physics(bandwidth,problem)
    I = problem.nTerminals;
    task_sizes = problem.terminalProperties.task_sizes;
    pt_dbm = problem.terminalProperties.Pt_dbm;
    fc = problem.terminalProperties.fc;
    cpu_rate = problem.fogNodeProperties.cpu_cycle_rate(1);
    finish = cell(1,I);
    energy = cell(1,I);
    for i = 1:I
        % d=max(raw_distance,1)=1 and p_LoS=1 by the dominance proof.
        path_loss = 32.8 + 20*log10(fc(i)/1e9) + ...
                    problem.fixed_shadow_LoS_val(i);
        noise_dbm = -174 + 10*log10(bandwidth(i)) + 5;
        snr = exp(log(10)*(pt_dbm(i)-path_loss-noise_dbm)/10);
        capacity = bandwidth(i)*log(1+snr)/log(2);
        t_trans = task_sizes(i)/capacity;
        service = task_sizes(i)*1000/cpu_rate;
        finish{i} = t_trans + service;
        pt_watt = 10^(pt_dbm(i)/10)/1e3;
        energy{i} = (pt_watt/0.2 + 5e-6*bandwidth(i))*t_trans + ...
                    0.1e-6*task_sizes(i);
    end
    g1 = sum([finish{:}])/I;
    g2 = sum([energy{:}]);
end

function opts = default_options(user)
    opts = struct('MaxTimePerCase', 120, ...
                  'AbsoluteTolerance', 1e-7, ...
                  'RelativeTolerance', 1e-7, ...
                  'PrintLevel', 0, ...
                  'ProgressEvery', 15, ...
                  'StopOnIncomplete', false, ...
                  'ScratchDirectory', fullfile(tempdir, ...
                                               'baron_micro_scratch'));
    names = fieldnames(user);
    for i = 1:numel(names)
        opts.(names{i}) = user.(names{i});
    end
end

function validate_problem(problem, weights)
    required = {'nTerminals','nFogNodes','area','Tslot', ...
                'systemTotalBandwidth','terminalProperties', ...
                'fogNodeProperties','bounds','fixed_shadow_LoS_val', ...
                'fixed_shadow_NLoS_val'};
    for k = 1:numel(required)
        if ~isfield(problem, required{k})
            error('baron_micro:MissingField', ...
                'problem.%s is required.', required{k});
        end
    end
    if problem.nTerminals ~= 2 || problem.nFogNodes ~= 2
        error('baron_micro:FreeLicenseSize', ...
            ['This certified implementation is intentionally configured for ', ...
             'I=2 and M=2 so every partition stays within the free BARON ', ...
             '10-variable/10-constraint limit.']);
    end
    if any(weights < 0 | weights > 1)
        error('baron_micro:InvalidWeights', ...
              'All delay weights must lie in [0,1].');
    end
    if sum(problem.bounds.bandwidth(1, :)) > ...
            problem.systemTotalBandwidth + 1e-9
        error('baron_micro:InfeasibleBandwidthBounds', ...
              'The sum of bandwidth lower bounds exceeds the system budget.');
    end
end

function [cases, stats] = build_complete_case_partition(problem)
    I = problem.nTerminals;
    M = problem.nFogNodes;
    total_assignments = M^I;
    assignment_ids = (0:total_assignments-1)';
    assignments = zeros(total_assignments, I);
    for i = 1:I
        assignments(:, i) = mod(floor(assignment_ids / M^(i-1)), M) + 1;
    end

    cases = cell(0, 1);
    cpu_infeasible = 0;
    task_cycles = problem.terminalProperties.task_sizes(:)' * 1000;
    cpu_capacity = problem.fogNodeProperties.cpu_cycle_rate(:)' * problem.Tslot;

    % Exact bound propagation determines which distance regions can exist.
    % This is not heuristic pruning: a region is omitted only when the
    % rectangular deployment bounds prove it empty for that terminal.
    region_options = cell(1, I);
    for i = 1:I
        tx = problem.terminalProperties.positions(i,1);
        ty = problem.terminalProperties.positions(i,2);
        corners = [problem.area(1,1), problem.area(2,1); ...
                   problem.area(1,1), problem.area(2,2); ...
                   problem.area(1,2), problem.area(2,1); ...
                   problem.area(1,2), problem.area(2,2)];
        max_distance = max(sqrt((corners(:,1)-tx).^2 + ...
                                (corners(:,2)-ty).^2));
        possible = 1; % raw distance <= 1 is always possible inside the area.
        if max_distance >= 1-1e-12
            possible(end+1) = 2; %#ok<AGROW>
        end
        if max_distance >= 5-1e-12
            possible(end+1) = 3; %#ok<AGROW>
        end
        region_options{i} = possible;
    end
    distance_patterns = cartesian_region_patterns(region_options);
    n_distance_patterns = size(distance_patterns, 1);

    next_case = 0;
    for a = 1:size(assignments, 1)
        offloading = assignments(a, :);

        exact_cpu_feasible = true;
        for node = 1:M
            exact_cpu_feasible = exact_cpu_feasible && ...
                sum(task_cycles(offloading == node)) <= cpu_capacity(node) + 1e-9;
        end
        if ~exact_cpu_feasible
            % This is exact constraint propagation, not heuristic pruning:
            % the complete continuous domain is infeasible for this assignment.
            cpu_infeasible = cpu_infeasible + 1;
            continue;
        end

        schedules = enumerate_all_node_orders(offloading, M);
        active_nodes = sum(arrayfun(@(m) any(offloading == m), 1:M));
        n_queue_events = I - active_nodes;

        for s = 1:numel(schedules)
            for queue_id = 0:(2^n_queue_events - 1)
                queue_modes = zeros(1, n_queue_events);
                for q = 1:n_queue_events
                    queue_modes(q) = bitget(queue_id, q);
                end

                for d = 1:n_distance_patterns
                    next_case = next_case + 1;
                    current = struct();
                    current.Id = next_case;
                    current.Offloading = offloading;
                    current.Schedule = schedules{s};
                    current.QueueModes = queue_modes;
                    current.DistanceRegions = distance_patterns(d, :);
                    current.NonlinearConstraintCount = ...
                        count_nonlinear_constraints(current, M);
                    cases{next_case, 1} = current; %#ok<AGROW> 
                end
            end
        end
    end

    stats = struct();
    stats.TotalAssignments = total_assignments;
    stats.CPUInfeasibleAssignments = cpu_infeasible;
    stats.PossibleDistanceRegions = region_options;
    stats.TotalAlgebraicCases = numel(cases);
end

function patterns = cartesian_region_patterns(region_options)
    patterns = zeros(1, 0);
    for i = 1:numel(region_options)
        values = region_options{i};
        expanded = zeros(size(patterns,1)*numel(values), i);
        next = 0;
        for r = 1:size(patterns,1)
            for v = 1:numel(values)
                next = next + 1;
                if i > 1
                    expanded(next,1:i-1) = patterns(r,:);
                end
                expanded(next,i) = values(v);
            end
        end
        patterns = expanded;
    end
end

function schedules = enumerate_all_node_orders(offloading, M)
    schedules = {cell(1, M)};
    for node = 1:M
        assigned = find(offloading == node);
        if isempty(assigned)
            variants = {[]};
        elseif isscalar(assigned)
            variants = {assigned};
        else
            P = perms(assigned);
            variants = cell(size(P, 1), 1);
            for r = 1:size(P, 1)
                variants{r} = P(r, :);
            end
        end

        expanded = cell(0, 1);
        for s = 1:numel(schedules)
            for v = 1:numel(variants)
                candidate = schedules{s};
                candidate{node} = variants{v};
                expanded{end+1, 1} = candidate; %#ok<AGROW>
            end
        end
        schedules = expanded;
    end
end

function n = count_nonlinear_constraints(case_data, M)
    % One inequality for regions 1/3 and two for region 2.
    n_distance = sum(1 + (case_data.DistanceRegions == 2));
    n_order = 0;
    n_deadline = 0;
    for node = 1:M
        order = case_data.Schedule{node};
        if ~isempty(order)
            n_order = n_order + max(0, numel(order)-1);
            n_deadline = n_deadline + 1;
        end
    end
    n_queue = numel(case_data.QueueModes);
    n = n_distance + n_order + n_queue + n_deadline;
end

function result = solve_across_all_cases(problem, cases, mode, scaling, opts)
    I = problem.nTerminals;
    M = problem.nFogNodes;
    lb = [problem.bounds.bandwidth(1, :), ...
          repmat([problem.area(1,1), problem.area(2,1)], 1, M)]';
    ub = [problem.bounds.bandwidth(2, :), ...
          repmat([problem.area(1,2), problem.area(2,2)], 1, M)]';
    x0 = [mean(problem.bounds.bandwidth, 1), ...
          problem.initial_fog_deployment_flat]';

    % The total-bandwidth constraint is an inequality.  The old local
    % "exact" solver incorrectly imposed equality to 225 MHz even though
    % the per-MTD upper bounds made that equality impossible.
    A = [ones(1, I), zeros(1, 2*M)];
    rl = -Inf;
    ru = problem.systemTotalBandwidth;

    if ~exist(opts.ScratchDirectory, 'dir')
        mkdir(opts.ScratchDirectory);
    end
    baron_options = baronset( ...
        'PrLevel', opts.PrintLevel, ...
        'MaxTime', opts.MaxTimePerCase, ...
        'EpsA', opts.AbsoluteTolerance, ...
        'EpsR', opts.RelativeTolerance, ...
        'barscratch', opts.ScratchDirectory);

    best_upper = Inf;
    global_lower = Inf;
    best_x = [];
    best_case = [];
    best_info = struct();
    n_optimal = 0;
    n_infeasible = 0;
    n_incomplete = 0;
    n_missing_lower_bounds = 0;
    total_wall_time = 0;
    statuses = strings(numel(cases), 1);

    for c = 1:numel(cases)
        case_data = cases{c};
        fun = @(x) scalar_objective(x, problem, case_data, mode, scaling);
        nlcon = @(x) partition_constraints(x, problem, case_data);
        cl = -Inf(case_data.NonlinearConstraintCount, 1);
        cu = zeros(case_data.NonlinearConstraintCount, 1);

        try
            [x, fval, exitflag, info] = baron( ...
                fun, A, rl, ru, lb, ub, nlcon, cl, cu, [], x0, ...
                baron_options);
        catch ME
            throwAsCaller(MException('baron_micro:BARONFailure', ...
                'BARON failed for algebraic case %d/%d. Original error: %s', ...
                c, numel(cases), ME.message));
        end

        if isfield(info, 'Wall_Time') && isfinite(info.Wall_Time)
            total_wall_time = total_wall_time + info.Wall_Time;
        end

        statuses(c) = string(info.Model_Status);
        if exitflag == 1
            n_optimal = n_optimal + 1;
            global_lower = min(global_lower, info.Lower_Bound);
            if fval < best_upper
                best_upper = fval;
                best_x = x;
                best_case = case_data;
                best_info = info;
            end
        elseif exitflag == 2
            n_infeasible = n_infeasible + 1;
        else
            n_incomplete = n_incomplete + 1;
            fprintf(['    WARNING case %d incomplete: exitflag=%s, ', ...
                     'status=%s\n'], c, mat2str(exitflag), ...
                    string(info.Model_Status));
            if opts.StopOnIncomplete
                error('baron_micro:IncompleteCase', ...
                    ['BARON did not certify algebraic case %d/%d ', ...
                     '(exitflag=%s, status=%s).'], ...
                    c, numel(cases), mat2str(exitflag), ...
                    char(string(info.Model_Status)));
            end
            if isfield(info, 'Lower_Bound') && ...
                    isscalar(info.Lower_Bound) && isfinite(info.Lower_Bound)
                global_lower = min(global_lower, info.Lower_Bound);
            else
                n_missing_lower_bounds = n_missing_lower_bounds + 1;
            end
            if ~isempty(fval) && isscalar(fval) && ...
                    isfinite(fval) && fval < best_upper
                best_upper = fval;
                best_x = x;
                best_case = case_data;
                best_info = info;
            end
        end

        if opts.ProgressEvery > 0 && ...
                (mod(c, opts.ProgressEvery) == 0 || c == numel(cases))
            fprintf('    Cases %d/%d | optimal %d | infeasible %d | incomplete %d\n', ...
                c, numel(cases), n_optimal, n_infeasible, n_incomplete);
        end
    end

    if isempty(best_x)
        error('baron_micro:NoFeasibleSolution', ...
              'No feasible solution was found for scalarization %s.', mode);
    end

    certificate_tolerance = opts.AbsoluteTolerance + ...
        opts.RelativeTolerance*max(1, abs(best_upper));
    certificate_complete = n_missing_lower_bounds == 0 && ...
        isfinite(global_lower) && ...
        best_upper-global_lower <= certificate_tolerance;

    [model_objectives, model_tmax] = numeric_case_result( ...
        best_x, problem, best_case);
    [verified_objectives, verified_tmax] = verify_against_evaluate_particle( ...
        best_x, problem, best_case, model_objectives, model_tmax);

    result = struct();
    result.Mode = mode;
    result.X = best_x;
    result.Objectives = verified_objectives;
    result.Tmax = verified_tmax;
    result.Offloading = best_case.Offloading;
    result.Schedule = best_case.Schedule;
    result.DistanceRegions = best_case.DistanceRegions;
    result.QueueModes = best_case.QueueModes;
    result.BestCaseId = best_case.Id;
    result.UpperBound = best_upper;
    result.LowerBound = global_lower;
    result.AbsoluteCertificateGap = best_upper - global_lower;
    result.CertificateComplete = certificate_complete;
    result.OptimalCases = n_optimal;
    result.InfeasibleCases = n_infeasible;
    result.IncompleteCases = n_incomplete;
    result.IncompleteCasesWithoutLowerBound = n_missing_lower_bounds;
    result.TotalCases = numel(cases);
    result.TotalBARONWallTime = total_wall_time;
    result.BestBARONInfo = best_info;
    result.CaseStatuses = statuses;

    fprintf(['    Global bounds: LB=%.12g, UB=%.12g, gap=%.3g, ', ...
             'certified=%d\n'], ...
        result.LowerBound, result.UpperBound, ...
        result.AbsoluteCertificateGap, result.CertificateComplete);
end

function f = scalar_objective(x, problem, case_data, mode, scaling)
    [g1, g2] = partition_physics(x, problem, case_data);
    switch mode
        case 'delay'
            f = g1;
        case 'energy'
            f = g2;
        case 'weighted'
            f = scaling.Weight(1) * ...
                ((g1 - scaling.Ideal(1)) / scaling.Range(1)) + ...
                scaling.Weight(2) * ...
                ((g2 - scaling.Ideal(2)) / scaling.Range(2));
        otherwise
            error('baron_micro:UnknownObjectiveMode', ...
                  'Unknown scalar objective mode: %s.', mode);
    end
end

function g = partition_constraints(x, problem, case_data)
    [~, ~, aux] = partition_physics(x, problem, case_data);
    constraint_cells = cell(0, 1);

    % Exact distance partition for d=max(raw_distance,1) and p_LoS(d).
    for i = 1:problem.nTerminals
        r2 = aux.DistanceSquared{i};
        region = case_data.DistanceRegions(i);
        if region == 1
            constraint_cells{end+1, 1} = r2 - 1; %#ok<AGROW>
        elseif region == 2
            constraint_cells{end+1, 1} = 1 - r2; %#ok<AGROW>
            constraint_cells{end+1, 1} = r2 - 25; %#ok<AGROW>
        else
            constraint_cells{end+1, 1} = 25 - r2; %#ok<AGROW>
        end
    end

    % FCFS arrival order, exact queue-state partition, and slot deadline.
    constraint_cells = [constraint_cells; aux.OrderRelations; ...
                        aux.QueueRelations; aux.DeadlineRelations];
    g = vertcat(constraint_cells{:});

    if length(g) ~= case_data.NonlinearConstraintCount
        error('baron_micro:ConstraintCountMismatch', ...
            'Expected %d nonlinear constraints but constructed %d.', ...
            case_data.NonlinearConstraintCount, length(g));
    end
end

function [G1, G2, aux] = partition_physics(x, problem, case_data)
    I = problem.nTerminals;
    M = problem.nFogNodes;
    bandwidth = x(1:I);
    task_sizes = problem.terminalProperties.task_sizes;
    task_cycles = task_sizes * 1000;
    term_pos = problem.terminalProperties.positions;
    pt_dbm = problem.terminalProperties.Pt_dbm;
    fc = problem.terminalProperties.fc;
    cpu = problem.fogNodeProperties.cpu_cycle_rate;

    t_trans = cell(1, I);
    energy = cell(1, I);
    distance_squared = cell(I, 1);

    for i = 1:I
        node = case_data.Offloading(i);
        base = I + 2*(node-1);
        dx = term_pos(i,1) - x(base+1);
        dy = term_pos(i,2) - x(base+2);
        r2 = dx.^2 + dy.^2;
        distance_squared{i} = r2;

        region = case_data.DistanceRegions(i);
        if region == 1
            d = 1;
            p_los = 1;
        elseif region == 2
            d = sqrt(r2);
            p_los = 1;
        else
            d = sqrt(r2);
            p_los = exp(-(d-5)/65);
        end

        pl_los = 16.9*log10(d) + 32.8 + ...
                 20*log10(fc(i)/1e9) + problem.fixed_shadow_LoS_val(i);
        pl_nlos = 38.3*log10(d) + 17.3 + ...
                  24.9*log10(fc(i)/1e9) + problem.fixed_shadow_NLoS_val(i);
        path_loss = p_los*pl_los + (1-p_los)*pl_nlos;

        noise_dbm = -174 + 10*log10(bandwidth(i)) + 5;
        snr_linear = exp(log(10) * ...
            (pt_dbm(i) - path_loss - noise_dbm) / 10);
        capacity = bandwidth(i) * log(1 + snr_linear) / log(2);
        t_trans{i} = task_sizes(i) / capacity;

        pt_watt = 10^(pt_dbm(i)/10) / 1e3;
        energy{i} = (pt_watt/0.2 + 5e-6*bandwidth(i)) * ...
                    t_trans{i} + 0.1e-6*task_sizes(i);
    end

    finish = cell(1, I);
    order_relations = cell(0, 1);
    queue_relations = cell(0, 1);
    deadline_relations = cell(0, 1);
    node_last = cell(0, 1);
    queue_index = 0;

    for node = 1:M
        order = case_data.Schedule{node};
        if isempty(order)
            continue;
        end

        for j = 1:(numel(order)-1)
            order_relations{end+1, 1} = ...
                t_trans{order(j)} - t_trans{order(j+1)}; %#ok<AGROW>
        end

        first = order(1);
        previous_finish = t_trans{first} + task_cycles(first)/cpu(node);
        finish{first} = previous_finish;

        for j = 2:numel(order)
            queue_index = queue_index + 1;
            current = order(j);
            arrival = t_trans{current};
            service = task_cycles(current)/cpu(node);

            if case_data.QueueModes(queue_index) == 0
                % Idle branch: the task arrives after the previous job ends.
                queue_relations{end+1, 1} = ...
                    previous_finish - arrival; %#ok<AGROW>
                current_finish = arrival + service;
            else
                % Busy branch: the task has already arrived and waits.
                queue_relations{end+1, 1} = ...
                    arrival - previous_finish; %#ok<AGROW>
                current_finish = previous_finish + service;
            end

            finish{current} = current_finish;
            previous_finish = current_finish;
        end

        node_last{end+1, 1} = previous_finish; %#ok<AGROW>
        deadline_relations{end+1, 1} = ...
            previous_finish - problem.Tslot; %#ok<AGROW>
    end

    G1 = sum([finish{:}]) / I;
    G2 = sum([energy{:}]);

    aux = struct();
    aux.DistanceSquared = distance_squared;
    aux.TransmissionTime = t_trans;
    aux.FinishTime = finish;
    aux.NodeLastFinish = node_last;
    aux.OrderRelations = order_relations;
    aux.QueueRelations = queue_relations;
    aux.DeadlineRelations = deadline_relations;
end

function [objectives, tmax] = numeric_case_result(x, problem, case_data)
    [g1, g2, aux] = partition_physics(x, problem, case_data);
    objectives = double([g1, g2]);
    tmax = max(double([aux.NodeLastFinish{:}]));
end

function [objectives, tmax] = verify_against_evaluate_particle( ...
        x, problem, case_data, model_objectives, model_tmax)
    I = problem.nTerminals;
    position = struct();
    position.bandwidth = x(1:I)';
    position.deployment = x(I+1:end)';
    position.offloading = case_data.Offloading;

    check = EvaluateParticle(position, problem);
    objectives = check.RawObjectives(1, :);
    tmax = check.Tmax(1);

    obj_scale = max(abs(objectives), [1e-9, 1e-12]);
    obj_error = max(abs(model_objectives-objectives) ./ obj_scale);
    tmax_error = abs(model_tmax-tmax) / max(abs(tmax), 1e-9);
    if obj_error > 1e-5 || tmax_error > 1e-5 || ~check.IsFeasible(1)
        error('baron_micro:PhysicsMismatch', ...
            ['Certified model and EvaluateParticle disagree. ', ...
             'Objective relative error %.3g, Tmax relative error %.3g, ', ...
             'EvaluateParticle feasible=%d.'], ...
            obj_error, tmax_error, check.IsFeasible(1));
    end
end

function require_complete_certificate(result, label)
    if ~result.CertificateComplete
        error('baron_micro:IncompleteCertificate', ...
            ['The %s was not globally certified across all cases: ', ...
             '%d incomplete of %d. Increase MaxTimePerCase or inspect ', ...
             'the returned case statuses before using the result.'], ...
            label, result.IncompleteCases, result.TotalCases);
    end
end

function score = scalar_score(objectives, scaling)
    normalized = (objectives - scaling.Ideal) ./ scaling.Range;
    score = normalized * scaling.Weight(:);
end

function nd = nondominated_rows(values)
    n = size(values, 1);
    keep = true(n, 1);
    for i = 1:n
        for j = 1:n
            if i ~= j && all(values(j,:) <= values(i,:)) && ...
                    any(values(j,:) < values(i,:))
                keep(i) = false;
                break;
            end
        end
    end
    nd = values(keep, :);
end
