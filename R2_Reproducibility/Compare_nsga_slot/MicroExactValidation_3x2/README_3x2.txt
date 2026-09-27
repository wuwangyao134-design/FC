3-MTD / 2-FN CERTIFIED MICRO-INSTANCE VALIDATION
=================================================

Purpose
-------
This is a second, independent exact-validation experiment. Unlike the
2-MTD/2-FN correctness instance, I=3>M=2 forces at least two MTDs to share
an FN. Association, heterogeneous FN rates, FCFS ordering, queue states,
and continuous bandwidth allocation therefore remain active.

Exact reference construction
----------------------------
1. All 2^3=8 MTD-to-FN associations are included.
2. Every FCFS arrival order for each association is included.
3. Every algebraic branch of max(arrival, previous finish) is included.
4. BARON globally certifies every algebraic case for the delay endpoint.
5. The separable energy endpoint is globally certified by the strictly
   increasing derivative-sign function and unique-root bisection.
6. No heuristic ranking, Top-K case pruning, or deployment sampling is used.

The work-cell size is 0.5 m by 0.5 m. Its diagonal is below the model's
1 m distance floor, so d=max(raw distance,1 m)=1 everywhere. Deployment
coordinates are consequently physically neutral for this micro instance;
this exact bound propagation does not remove any attainable objective value.
The instance still contains unavoidable queue contention because I>M.

Files
-----
Run_Micro_Exact_Validation_3x2.m
baron_unpruned_micro_solver_3x2.m

Recommended execution
---------------------
First run the smoke test:

    clear functions;
    rehash;
    validation3x2 = Run_Micro_Exact_Validation_3x2('quick');

Only if all three scalarizations are certified and the quick run finishes,
start the formal 30-run experiment:

    clear functions;
    rehash;
    validation3x2 = Run_Micro_Exact_Validation_3x2('formal');

The requested delay weights are fixed at [0, 1], representing the exact
energy and delay endpoints. The complementary 2-MTD/2-FN experiment already
provides the BARON-certified balanced preference w=0.5. The 3-MTD/2-FN
experiment is intentionally used to add unavoidable association and queue
contention without reporting any uncertified intermediate scalarization.

Output
------
Results are stored under:

    Compare_nsga_slot\MicroExactValidation_3x2_Output\<timestamp>

The QUICK mode is a smoke test only and must not be reported in the paper.
