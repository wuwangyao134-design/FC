function Results = EvaluateParticle(input_data, problem)
%EVALUATEPARTICLE Evaluate latency, energy, and constraint violation.
%   The legacy penalized Objectives field is retained for all existing
%   evolutionary algorithms.  RawObjectives, ConstraintViolation, and
%   IsFeasible are additionally returned for learning algorithms so that
%   constraint information is not hidden behind the 1e9 hard penalty.

    if isstruct(input_data) && isfield(input_data, 'Position')
        pop_positions = input_data;
    elseif isstruct(input_data) && isfield(input_data, 'deployment')
        temp_pop.Position = input_data;
        pop_positions = temp_pop;
    else
        error('EvaluateParticle:UnsupportedInput', ...
              'The input must be a population or a Position structure.');
    end

    n_pop = numel(pop_positions);
    n_term = problem.nTerminals;
    n_fog = problem.nFogNodes;

    term_pos = problem.terminalProperties.positions;
    task_sizes = problem.terminalProperties.task_sizes;
    pt_dbm = problem.terminalProperties.Pt_dbm;
    fc = problem.terminalProperties.fc;
    cpu_rates = problem.fogNodeProperties.cpu_cycle_rate;
    t_slot = problem.Tslot;
    system_bw = problem.systemTotalBandwidth;
    shadow_los = problem.fixed_shadow_LoS_val;
    shadow_nlos = problem.fixed_shadow_NLoS_val;

    pt_watt = 10.^(pt_dbm / 10) / 1e3;
    computation_energy = 0.1e-6 .* task_sizes;
    task_cycles = task_sizes * 1000;

    all_objectives = zeros(n_pop, 2);
    all_raw_objectives = zeros(n_pop, 2);
    all_tmax = zeros(n_pop, 1);
    all_violation = zeros(n_pop, 1);
    all_feasible = false(n_pop, 1);

    for p = 1:n_pop
        current = pop_positions(p).Position;
        fog_pos = reshape(current.deployment, [2, n_fog])';
        bandwidth = current.bandwidth;
        offloading = round(current.offloading);

        if numel(offloading) ~= n_term || any(offloading < 1 | offloading > n_fog)
            error('EvaluateParticle:InvalidOffloading', ...
                  'Every offloading index must be an integer in [1,nFogNodes].');
        end

        bandwidth = max(bandwidth, eps);
        if sum(bandwidth) > system_bw
            bandwidth = bandwidth * (system_bw / sum(bandwidth));
        end

        selected_fog = fog_pos(offloading, :);
        distance = sqrt(sum((term_pos - selected_fog).^2, 2))';
        distance = max(distance, 1);

        p_los = (distance < 5) + ...
                (distance >= 5) .* exp(-(distance - 5) / 65);
        pl_los = 16.9*log10(distance) + 32.8 + ...
                 20*log10(fc/1e9) + shadow_los;
        pl_nlos = 38.3*log10(distance) + 17.3 + ...
                  24.9*log10(fc/1e9) + shadow_nlos;
        path_loss = p_los .* pl_los + (1-p_los) .* pl_nlos;

        noise_dbm = -174 + 10*log10(bandwidth)+5;
        snr_linear = 10.^((pt_dbm - path_loss - noise_dbm) / 10);
        capacity = bandwidth .* log2(1 + snr_linear);
        transmission_time = task_sizes ./ capacity;
        transmission_time(capacity <= 0 | ~isfinite(capacity)) = 1e6;

        finish_time = zeros(1, n_term);
        node_tmax = 0;
        cpu_violation = 0;

        for node = 1:n_fog
            assigned = (offloading == node);
            if ~any(assigned)
                continue;
            end

            node_cycles = sum(task_cycles(assigned));
            node_capacity = cpu_rates(node) * t_slot;
            cpu_violation = cpu_violation + ...
                max(0, node_cycles - node_capacity) / max(node_capacity, eps);

            [arrivals, order] = sort(transmission_time(assigned));
            service = task_cycles(assigned) / cpu_rates(node);
            service = service(order);
            completions = zeros(1, numel(arrivals));
            last_finish = 0;
            for j = 1:numel(arrivals)
                last_finish = max(arrivals(j), last_finish) + service(j);
                completions(j) = last_finish;
            end

            local_finish = finish_time(assigned);
            local_finish(order) = completions;
            finish_time(assigned) = local_finish;
            node_tmax = max(node_tmax, last_finish);
        end

        raw_delay = mean(finish_time);
        transmission_energy = ...
            (pt_watt / 0.2 + 5e-6 .* bandwidth) .* transmission_time;
        raw_energy = sum(transmission_energy + computation_energy);

        time_violation = max(0, node_tmax / t_slot - 1);
        total_violation = time_violation + cpu_violation;
        feasible = isfinite(raw_delay) && isfinite(raw_energy) && ...
                   total_violation <= 1e-10;

        penalized = [raw_delay, raw_energy];
        if ~feasible
            penalized = penalized + 1e9;
        end

        all_objectives(p, :) = penalized;
        all_raw_objectives(p, :) = [raw_delay, raw_energy];
        all_tmax(p) = node_tmax;
        all_violation(p) = total_violation;
        all_feasible(p) = feasible;
    end

    Results.Objectives = all_objectives;
    Results.RawObjectives = all_raw_objectives;
    Results.Tmax = all_tmax;
    Results.ConstraintViolation = all_violation;
    Results.IsFeasible = all_feasible;
end
