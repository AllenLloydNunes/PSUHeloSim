%% test_step1_freefall.m
% STEP 1 (PLAN.md Sec.9): EOM + gravity only, contact DISABLED.
%
% PASS CRITERION: the simulated inertial trajectory must match the analytic
% free-fall parabola
%       P_i(t) = P_i(0) + V_i(0)*t + 0.5*[0;0;g]*t^2
%
% This is a MACHINE-PRECISION test, not an engineering-tolerance one. RK4
% integrates constant-acceleration dynamics exactly, so any error above
% roundoff indicates a real defect in the EOM, the rotation matrix, or the
% gravity transform.
%
% Case 1c (tilted) is the important one: with no moments acting, the attitude
% must stay frozen AND the inertial acceleration must still be exactly
% [0;0;g] despite gravity being resolved in a rotated body frame. That
% simultaneously validates rot_b2i and the R'*[0;0;mg] gravity transform.

clear; clc;
fprintf('================================================================\n');
fprintf(' STEP 1 -- free-fall verification (contact disabled)\n');
fprintf('================================================================\n');

P = landing_params();
P.contact_on = false;
dt   = P.dt_physics;
tEnd = 2.0;

TOL_POS  = 1e-9;    % m      -- generous vs expected ~1e-13 roundoff
TOL_ATT  = 1e-12;   % rad    -- attitude must be frozen
TOL_RATE = 1e-12;   % rad/s  -- rates must stay identically zero

cases = struct('name', {}, 'x0', {});

% 1a -- level, released from rest at 5 m altitude
x0 = zeros(12,1);
x0(12) = -5.0;
cases(1) = struct('name', '1a  level, at rest', 'x0', x0);

% 1b -- level, translating (also exercises the position kinematics in x,y)
x0 = zeros(12,1);
x0(1) = 3.0;  x0(2) = 1.0;
x0(12) = -5.0;
cases(2) = struct('name', '1b  level, translating', 'x0', x0);

% 1c -- tilted and translating: validates rot_b2i + gravity transform
x0 = zeros(12,1);
x0(1:3) = [3.0; 1.0; -2.0];
x0(7:9) = [20*pi/180; -15*pi/180; 30*pi/180];
x0(12)  = -5.0;
cases(3) = struct('name', '1c  tilted, translating', 'x0', x0);

step1_pass = true;

for iC = 1:numel(cases)
    x0 = cases(iC).x0;
    S  = run_landing_sim(x0, tEnd, dt, P);

    % Analytic inertial trajectory
    R0   = rot_b2i(x0(7), x0(8), x0(9));
    V_i0 = R0 * x0(1:3);
    P_i0 = x0(10:12);
    tt   = S.t;
    Pana = P_i0.' + V_i0.'.*tt + 0.5*[0 0 P.g].*(tt.^2);

    errPos  = max(abs(S.x(:,10:12) - Pana), [], 'all');
    errAtt  = max(abs(S.x(:,7:9) - x0(7:9).'), [], 'all');
    errRate = max(abs(S.x(:,4:6)), [], 'all');

    ok = (errPos < TOL_POS) && (errAtt < TOL_ATT) && (errRate < TOL_RATE);
    step1_pass = step1_pass && ok;

    fprintf('\n  %s\n', cases(iC).name);
    fprintf('    max |position error| vs analytic : %.3e m      (tol %.1e)\n', errPos,  TOL_POS);
    fprintf('    max |attitude drift|             : %.3e rad    (tol %.1e)\n', errAtt,  TOL_ATT);
    fprintf('    max |body rate|                  : %.3e rad/s  (tol %.1e)\n', errRate, TOL_RATE);
    % Altitude is measured from the GROUND PLANE (z = 0), not from the start
    % point. Contact is disabled in this test, so the body falls straight
    % through z = 0 and a negative final altitude is expected and correct.
    fprintf('    start altitude                   : %+.4f m (above ground plane)\n', S.alt(1));
    fprintf('    final altitude                   : %+.4f m (no floor: contact disabled)\n', S.alt(end));
    fprintf('    net drop over %.1f s               : %.4f m\n', tEnd, S.alt(1) - S.alt(end));
    if ok
        fprintf('    [PASS]\n');
    else
        fprintf('    [FAIL]\n');
    end
end

fprintf('\n================================================================\n');
if step1_pass
    fprintf(' STEP 1: PASS -- free-fall matches analytic to machine precision\n');
else
    fprintf(' STEP 1: FAIL\n');
end
fprintf('================================================================\n');

% Persist the result. run_all_tests.m reads this rather than the workspace,
% because each test script clears the workspace on entry to stay standalone.
save('step1_result.mat', 'step1_pass');
