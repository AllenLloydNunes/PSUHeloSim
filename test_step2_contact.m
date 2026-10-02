%% test_step2_contact.m
% STEP 2 (PLAN.md Sec.9): normal spring-damper contact, friction OFF (mu = 0).
%
% PASS CRITERIA:
%   2A  drop settles to a static equilibrium
%   2B  at equilibrium: sum(N) == m*g, and net force AND net moment are zero
%   2C  initialised at equilibrium, the model stays there (no drift)
%   2D  a laterally symmetric drop produces ZERO roll and ZERO yaw
%   2E  energy conserved with c = 0 in continuous contact (no liftoff)
%   2F  energy conserved with c = 0 through a full bounce (contact on/off)
%   plus the PLAN.md Sec.7 sanity checks on every run (no adhesion, etc.)
%
% NOTE ON "STAYS LEVEL": PLAN.md Sec.9 said the drop should "settle to 5 cm
% static deflection; stay level". That is only half right, and the geometry is
% why. The feet sit at body x = +2, +2, -3, so sum(x_i) = +1 m: the leg
% triangle's centroid is NOT under the CG. Equal leg loads would therefore
% leave a net pitch moment of 3270 N*m. The true equilibrium loads unevenly
%   front legs 2943 N (delta = 0.045 m), aft leg 3924 N (delta = 0.060 m)
% giving a small nose-up pitch of ~0.003 rad. The MEAN deflection is still
% 0.05 m, so k is correct as specified.
%
% The geometry IS laterally symmetric (y = -1.5, +1.5, 0 about a CG at y = 0),
% so ROLL and YAW must remain identically zero -- and that is enforced here to
% 1e-10. Pitch is legitimately nonzero and is checked against the analytic
% small-angle prediction instead.

clear; clc;
fprintf('================================================================\n');
fprintf(' STEP 2 -- spring-damper contact verification (friction OFF)\n');
fprintf('================================================================\n');

P = landing_params();
P.contact_on = true;
P.mu = 0;                       % step 2: no friction
dt = P.dt_physics;

step2_pass = true;

%% EXACT static equilibrium, derived independently of the model code
% Unknowns: zD and theta (phi = psi = 0 follows from lateral symmetry).
%
%   delta_i = zD - sin(th)*x_i + cos(th)*z_i        (exact; row 3 of R_b2i)
%   N_i     = k*delta_i
%   R1 = sum(N_i) - m*g                        = 0  (vertical force balance)
%   R2 = sum(N_i*(z_i*sin(th) + x_i*cos(th)))  = 0  (pitch moment balance)
%
% THE z_i*sin(th) TERM IN R2 MATTERS. The feet sit 1.8 m BELOW the CG, so a
% pitch tilt gives the vertical normal forces a moment arm about y. A first,
% cruder version of this cross-check used only sum(x_i*N_i) = 0 and disagreed
% with the simulation by ~11 N in leg load and 5e-5 rad in pitch. The
% SIMULATION WAS RIGHT and the hand-check was wrong -- check 2B independently
% confirms the settled state satisfies force and moment balance to ~1e-6.
% The tolerances below are kept tight enough to have caught that error.
%
% Solved here by Newton iteration on the hand-written residual, so this stays
% an independent check rather than a circular one against contact_forces.m.

xv = P.p_foot_b(1,:);
zv = P.p_foot_b(3,:);
deltaOf = @(q) q(1) - sin(q(2))*xv + cos(q(2))*zv;
resid   = @(q) [ sum(P.k*deltaOf(q)) - P.m*P.g; ...
                 sum( P.k*deltaOf(q) .* (zv*sin(q(2)) + xv*cos(q(2))) ) ];

q = [-1.75; 0.003];                     % crude starting guess
converged = false;
for it = 1:100
    Rq = resid(q);
    if norm(Rq) < 1e-10, converged = true; break; end
    J = zeros(2,2);  hstep = 1e-7;
    for j = 1:2
        qp = q;  qp(j) = qp(j) + hstep;
        J(:,j) = (resid(qp) - Rq)/hstep;
    end
    q = q - J\Rq;
end
if ~converged
    error('test_step2:equilibrium', 'Newton solve for equilibrium did not converge.');
end

zD_pred    = q(1);
theta_pred = q(2);
d_pred     = deltaOf(q);
N_pred     = P.k*d_pred;

fprintf('\n  Exact equilibrium (independent Newton solve, %d iterations):\n', it);
fprintf('    leg loads       %8.2f  %8.2f  %8.2f N\n', N_pred);
fprintf('    leg deltas      %8.5f  %8.5f  %8.5f m\n', d_pred);
fprintf('    sum of loads    %8.2f N   (m*g = %.2f N)\n', sum(N_pred), P.m*P.g);
fprintf('    mean delta      %8.5f m   (design target 0.05000)\n', mean(d_pred));
fprintf('    pitch           %8.6f rad (%.4f deg)\n', theta_pred, theta_pred*180/pi);
fprintf('    CG zD           %8.6f m   (altitude %.6f m)\n', zD_pred, -zD_pred);

%% 2A -- drop and settle -------------------------------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 2A  drop from 0.2 m above first contact, settle under damping\n');
fprintf('----------------------------------------------------------------\n');

x0 = zeros(12,1);
x0(12) = -2.0;                  % feet 0.2 m above ground (contact at zD = -1.8)
tSettle = 8.0;
S_settle = run_landing_sim(x0, tSettle, dt, P);

R = check_sanity(S_settle, '2A settle');
step2_pass = step2_pass && R.pass;

xEq = S_settle.x(end,:).';
fprintf('    settled: zD = %.6f m, theta = %.6f rad (%.4f deg)\n', ...
        xEq(12), xEq(8), xEq(8)*180/pi);
fprintf('    settled leg loads : %8.2f  %8.2f  %8.2f N\n', S_settle.N(end,:));
fprintf('    settled deltas    : %8.5f  %8.5f  %8.5f m\n', S_settle.delta(end,:));
fprintf('    mean delta        : %8.5f m  (design target 0.05000)\n', mean(S_settle.delta(end,:)));

% Cross-check the settled state against the independent Newton solution
errN     = max(abs(S_settle.N(end,:) - N_pred));
errTheta = abs(xEq(8)  - theta_pred);
errZD    = abs(xEq(12) - zD_pred);
okPred   = (errN < 1e-2) && (errTheta < 1e-7) && (errZD < 1e-7);
step2_pass = step2_pass && okPred;
fprintf('    [%s] matches exact equilibrium: dN = %.3e N (tol 1e-2), ', ...
        tern(okPred,'PASS','FAIL'), errN);
fprintf('dtheta = %.3e rad (tol 1e-7), dzD = %.3e m (tol 1e-7)\n', errTheta, errZD);

%% 2B -- equilibrium conditions at the settled state ---------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 2B  static equilibrium: force balance, moment balance\n');
fprintf('----------------------------------------------------------------\n');

xEqStatic = xEq;
xEqStatic(1:6) = 0;             % zero all velocities and rates

[xdotEq, auxEq] = landing_eom(0, xEqStatic, P, @friction_physics);
% With omega = 0: xdot(1:3) = F_b/m and xdot(4:6) = Iinv*M_b, so the residual
% force and moment recover directly from the state derivative.
netF  = P.m * xdotEq(1:3);      % residual force  [N]
netM  = P.I * xdotEq(4:6);      % residual moment [N*m]

sumN     = sum(auxEq.N);
errWeight = abs(sumN - P.m*P.g);

okF  = norm(netF) < 1e-3;
okM  = norm(netM) < 1e-2;
okW  = errWeight  < 1e-3;
step2_pass = step2_pass && okF && okM && okW;

fprintf('    sum of normal forces : %.6f N  vs  m*g = %.6f N\n', sumN, P.m*P.g);
fprintf('    [%s] force balance : |sum(N) - m*g| = %.3e N   (tol 1e-3)\n', tern(okW,'PASS','FAIL'), errWeight);
fprintf('    [%s] net force     : |F| = %.3e N              (tol 1e-3)\n', tern(okF,'PASS','FAIL'), norm(netF));
fprintf('    [%s] net moment    : |M| = %.3e N*m            (tol 1e-2)\n', tern(okM,'PASS','FAIL'), norm(netM));

%% 2C -- hold at equilibrium ---------------------------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 2C  initialised at equilibrium, must not drift\n');
fprintf('----------------------------------------------------------------\n');

S_hold = run_landing_sim(xEqStatic, 2.0, dt, P);
R = check_sanity(S_hold, '2C hold');
step2_pass = step2_pass && R.pass;

driftPos = max(abs(S_hold.x(:,10:12) - xEqStatic(10:12).'), [], 'all');
driftAtt = max(abs(S_hold.x(:,7:9)  - xEqStatic(7:9).'),   [], 'all');
okHold = (driftPos < 1e-4) && (driftAtt < 1e-4);
step2_pass = step2_pass && okHold;
fprintf('    [%s] position drift %.3e m, attitude drift %.3e rad (tol 1e-4)\n', ...
        tern(okHold,'PASS','FAIL'), driftPos, driftAtt);

%% 2D -- lateral symmetry: zero roll, zero yaw ---------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 2D  symmetric drop must produce ZERO roll and ZERO yaw\n');
fprintf('----------------------------------------------------------------\n');

maxRoll = max(abs(S_settle.x(:,7)));
maxYaw  = max(abs(S_settle.x(:,9)));
maxP    = max(abs(S_settle.x(:,4)));
maxR    = max(abs(S_settle.x(:,6)));
okSym   = (maxRoll < 1e-10) && (maxYaw < 1e-10) && (maxP < 1e-10) && (maxR < 1e-10);
step2_pass = step2_pass && okSym;

fprintf('    max |roll|  = %.3e rad     max |p| = %.3e rad/s\n', maxRoll, maxP);
fprintf('    max |yaw|   = %.3e rad     max |r| = %.3e rad/s\n', maxYaw,  maxR);
fprintf('    max |pitch| = %.3e rad  (legitimately nonzero, see header)\n', max(abs(S_settle.x(:,8))));
fprintf('    [%s] roll and yaw identically zero (tol 1e-10)\n', tern(okSym,'PASS','FAIL'));

%% 2E -- energy conservation, continuous contact -------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 2E  energy conservation, c = 0, legs never leave the ground\n');
fprintf('----------------------------------------------------------------\n');

Pc = P;  Pc.c = 0;  Pc.mu = 0;
x0e = xEqStatic;
x0e(12) = xEqStatic(12) - 0.01;     % raise 1 cm -> oscillates without liftoff

S_econt = run_landing_sim(x0e, 3.0, dt, Pc);
R = check_sanity(S_econt, '2E energy (continuous contact)', ...
                 struct('checkEnergy', true, 'energyTol', 1e-8, ...
                        'requireAllContact', true));
step2_pass = step2_pass && R.pass;

%% 2F -- energy conservation through a bounce ----------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 2F  energy conservation, c = 0, full bounce (contact on/off)\n');
fprintf('----------------------------------------------------------------\n');

x0b = zeros(12,1);
x0b(12) = -2.0;
S_ebounce = run_landing_sim(x0b, 3.0, dt, Pc);
R = check_sanity(S_ebounce, '2F energy (with bounce)', ...
                 struct('checkEnergy', true, 'energyTol', 1e-3));
step2_pass = step2_pass && R.pass;
fprintf('    (looser tolerance than 2E is expected: RK4 loses order stepping\n');
fprintf('     across the contact discontinuity at touchdown/liftoff)\n');

%% Summary ---------------------------------------------------------------
fprintf('\n================================================================\n');
if step2_pass
    fprintf(' STEP 2: PASS -- all contact and conservation criteria met\n');
else
    fprintf(' STEP 2: FAIL\n');
end
fprintf('================================================================\n');

% Persist results. run_all_tests.m reads this rather than the workspace,
% because each test script clears the workspace on entry to stay standalone.
save('step2_result.mat', 'step2_pass', 'S_settle', 'S_econt', 'S_ebounce', 'xEqStatic');

function out = tern(c,a,b)
if c, out = a; else, out = b; end
end
