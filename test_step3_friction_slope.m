%% test_step3_friction_slope.m
% STEP 3: friction ON, on a 5-degree incline (terrain rises toward North).
%
%   3A  settle and HOLD on the slope      -- statics and the creep artifact
%   3B  slide to a stop, uphill & downhill -- kinetic friction vs analytic
%   3C  integrator stability vs dt, v_bl   -- linear prediction vs seeded runs
%
% ANALYTIC REFERENCES (quasi-static, no normal acceleration):
%   holding requires tan(alpha) < mu  (0.0875 < 0.6: holds)
%   sum of normal forces      = m g cos(alpha)
%   total friction            = m g sin(alpha), pointing UPHILL
%   regularized friction cannot hold perfectly still: to generate m g sin(a)
%   inside the boundary layer it needs a sliding speed, so it CREEPS downhill:
%        v_creep = v_bl * tan(alpha) / mu          (exact in steady state)
%   sliding deceleration: g*(mu cos a - sin a) downhill, g*(mu cos a + sin a) uphill
%
% STABILITY FINDING THAT SHAPES THIS TEST:
% PLAN.md Sec.6 estimated the friction pole from TRANSLATION only (-589 rad/s)
% and concluded v_bl = 0.01 is stable at dt = 0.001, and v_bl = 0.05 stable at
% dt = 0.01. Both conclusions are wrong. The feet sit 1.8 m below the CG, so
% friction also resists ROTATION, and the roll inertia is small: the effective
% mass seen by lateral friction is ~171 kg, not 1000 kg. Linearizing the model
% gives a lateral/roll friction pole of ~-3450 rad/s at v_bl = 0.01, which is
% RK4-UNSTABLE at dt = 0.001.
%
% A SYMMETRIC TEST HIDES THIS. On this slope roll, yaw, v, p, r are exactly
% zero by symmetry, so the unstable lateral mode is never seeded and a
% dt = 0.001 run looks perfect. 3C therefore seeds it deliberately.
%
% RESOLUTION: v_bl was raised to the UH-60 model's value, 0.1 ft/s = 0.03048
% m/s (LG_VBL, InitLG.m). The pole scales as 1/v_bl, so it drops ~3x and
% dt = 0.001 becomes stable. The VALIDATION runs (3A, 3B) therefore run at
% P.dt_physics = 0.001 -- but only after the linearization GATE confirms that
% dt is stable for the configured v_bl. If someone tightens v_bl again, the
% gate fails loudly instead of the tests silently passing on a symmetric run.
% 3C keeps the old v_bl = 0.01 as a documented unstable case.

clear; clc;
fprintf('================================================================\n');
fprintf(' STEP 3 -- friction on a 5 deg incline\n');
fprintf('================================================================\n');

P = landing_params();                 % mu = 0.6, v_bl = 0.03048 (UH-60 value)
SLOPE_DEG = 5;
P.terrain = make_terrain(SLOPE_DEG, 0);
a  = SLOPE_DEG*pi/180;
n  = P.terrain.n_down;
tu = P.terrain.t_up;
eL = cross(n, tu);                    % lateral unit vector in the plane (= East here)
g = P.g;  m = P.m;  mu = P.mu;

DT_VAL = P.dt_physics;                % validation timestep = planned physics dt
                                      % (0.001), GATED on linear stability below
Rk     = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;           % RK4 amplification
rhoOf  = @(lam, dt) max(abs(Rk(dt*lam(abs(lam) > 1))));   % ignore neutral modes

step3_pass = true;

%% Analytic predictions --------------------------------------------------
v_creep_pred = P.v_bl*tan(a)/mu;
sumN_pred    = m*g*cos(a);
fric_pred    = m*g*sin(a);
fprintf('\n  Analytic predictions (slope %g deg, mu %.2f, v_bl %.3f m/s):\n', SLOPE_DEG, mu, P.v_bl);
fprintf('    holds?            tan(a) = %.4f < mu = %.2f  -> YES\n', tan(a), mu);
fprintf('    sum of N          %.4f N   (m g cos a)\n', sumN_pred);
fprintf('    total friction    %.4f N   (m g sin a, uphill)\n', fric_pred);
fprintf('    creep speed       %.6e m/s downhill (v_bl tan a / mu)\n', v_creep_pred);

%% 3A -- settle and hold on the slope ------------------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 3A  settle and hold on the slope (dt = %g)\n', DT_VAL);
fprintf('----------------------------------------------------------------\n');

x0 = zeros(12,1);
x0(8)     = a;               % body aligned with slope, nose uphill (body z = n)
x0(10:12) = -1.75*n;         % feet ~5 cm into the ground along the normal
S3A = run_landing_sim(x0, 10.0, DT_VAL, P);
R = check_sanity(S3A, '3A settle on slope');
step3_pass = step3_pass && R.pass;

xS = S3A.x(end,:).';
RS = rot_b2i(xS(7), xS(8), xS(9));
[xdotS, auxS] = landing_eom(0, xS, P);

% Linearization validity: all legs loaded, every foot inside the boundary layer
footSpd = zeros(1, P.nLegs);
for i = 1:P.nLegs
    Vf = RS*( xS(1:3) + cross(xS(4:6), P.p_foot_b(:,i)) );
    footSpd(i) = norm(Vf - (n.'*Vf)*n);
end
linOK = all(auxS.contact) && all(auxS.N > 0) && all(footSpd < P.v_bl);
if ~linOK
    error('test_step3:linearization', 'Settled state is not a valid linearization point.');
end

% GATE: the validation timestep must be linearly stable, or 3A/3B are invalid
[~, lamS] = linearize_eom(xS, P);
rhoVal = rhoOf(lamS, DT_VAL);
okGate = rhoVal < 1;
step3_pass = step3_pass && okGate;
fprintf('    [%s] validation dt = %g is linearly stable: max|R(h*lambda)| = %.4f\n', ...
        tern(okGate,'PASS','FAIL'), DT_VAL, rhoVal);

Vi      = RS*xS(1:3);
v_creep = -tu.'*Vi;
Ftot    = sum(auxS.F_t, 2);
sumN    = sum(auxS.N);
Fbal    = Ftot - sumN*n + [0; 0; m*g];          % inertial: contact + gravity
netF    = m*norm(xdotS(1:3));
netM    = norm(P.I*xdotS(4:6));
holdRatio = max( sqrt(sum(auxS.F_t.^2,1)) ./ (mu*auxS.N) );

errCreep = abs(v_creep - v_creep_pred)/v_creep_pred;
okCreep  = errCreep < 1e-3 && abs(n.'*Vi) < 1e-8;
okN      = abs(sumN - sumN_pred)     < 1e-2;
okFric   = abs(tu.'*Ftot - fric_pred) < 1e-2 && abs(eL.'*Ftot) < 1e-9;
okBal    = norm(Fbal) < 1e-2 && netF < 1e-2 && netM < 1e-2;
okHold   = holdRatio < 1;
maxLat   = max(abs(S3A.x(:,[2 4 6 7 9])), [], 'all');
okSym    = maxLat < 1e-10;
step3_pass = step3_pass && okCreep && okN && okFric && okBal && okHold && okSym;

fprintf('    settled pitch %.5f deg  (slope %g deg, offset %.4f deg from uneven leg loads)\n', ...
        xS(8)*180/pi, SLOPE_DEG, (xS(8)-a)*180/pi);
fprintf('    leg loads %8.2f %8.2f %8.2f N\n', auxS.N);
fprintf('    [%s] creep speed      %.6e m/s vs %.6e analytic (rel err %.2e, tol 1e-3)\n', ...
        tern(okCreep,'PASS','FAIL'), v_creep, v_creep_pred, errCreep);
fprintf('    [%s] sum of N         %.4f vs %.4f N  (|err| %.2e, tol 1e-2)\n', ...
        tern(okN,'PASS','FAIL'), sumN, sumN_pred, abs(sumN-sumN_pred));
fprintf('    [%s] friction uphill  %.4f vs %.4f N  (|err| %.2e, tol 1e-2)\n', ...
        tern(okFric,'PASS','FAIL'), tu.'*Ftot, fric_pred, abs(tu.'*Ftot-fric_pred));
fprintf('    [%s] force/moment balance: |F_inertial| %.2e N, |F_body| %.2e N, |M| %.2e N*m (tol 1e-2)\n', ...
        tern(okBal,'PASS','FAIL'), norm(Fbal), netF, netM);
fprintf('    [%s] holding, not sliding: max |F_t|/(mu N) = %.4f  (must be < 1)\n', ...
        tern(okHold,'PASS','FAIL'), holdRatio);
fprintf('    [%s] symmetry: max |v, p, r, roll, yaw| = %.2e  (tol 1e-10)\n', ...
        tern(okSym,'PASS','FAIL'), maxLat);

%% 3B -- slide to a stop --------------------------------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 3B  slide to a stop from 3 m/s (dt = %g)\n', DT_VAL);
fprintf('----------------------------------------------------------------\n');
fprintf('    Compared: distance until speed first drops below 0.10 m/s (>3x v_bl,\n');
fprintf('    so the boundary layer plays no part). Tolerance 2%%: the analytic result\n');
fprintf('    assumes sum(N) = m g cos a throughout, which the pitch/heave transient\n');
fprintf('    from sudden friction onset violates slightly.\n');

V0 = 3.0;  VTHR = 0.10;  TOL_STOP = 0.02;
if VTHR < 3*P.v_bl
    error('test_step3:vthr', 'Stop threshold must sit well above v_bl or the boundary layer contaminates the comparison.');
end
slides = struct('name', {'downhill (tail-first)', 'uphill (nose-first)'}, ...
                'dir',  {-tu, tu}, ...
                'acc',  {g*(mu*cos(a) - sin(a)), g*(mu*cos(a) + sin(a))});
S3B = cell(1, numel(slides));
for c = 1:numel(slides)
    d  = slides(c).dir;
    x0s = xS;
    x0s(1:3) = RS.' * (V0*d);          % prescribe inertial velocity along the slope
    x0s(4:6) = 0;
    Ss = run_landing_sim(x0s, 1.5, DT_VAL, P);
    S3B{c} = Ss;
    R = check_sanity(Ss, ['3B ' slides(c).name]);
    step3_pass = step3_pass && R.pass;

    Vtr = inertialVel(Ss);
    spd = sqrt(sum((Vtr - (Vtr*n)*n.').^2, 2));
    kStop = find(spd < VTHR, 1);
    if isempty(kStop)
        error('test_step3:noStop', 'Slide "%s" never slowed below %.2f m/s.', slides(c).name, VTHR);
    end
    dSim  = d.' * (Ss.x(kStop,10:12).' - Ss.x(1,10:12).');
    dPred = (V0^2 - VTHR^2)/(2*slides(c).acc);
    tPred = (V0 - VTHR)/slides(c).acc;
    errD  = abs(dSim - dPred)/dPred;
    symB  = max(abs(Ss.x(:,[2 4 6 7 9])), [], 'all');
    ok    = errD < TOL_STOP && symB < 1e-10;
    step3_pass = step3_pass && ok;
    fprintf('    [%s] %-22s decel %.4f m/s^2 | distance %.4f vs %.4f m (err %.2f%%) | time %.4f vs %.4f s | sym %.1e\n', ...
            tern(ok,'PASS','FAIL'), slides(c).name, slides(c).acc, dSim, dPred, 100*errD, ...
            Ss.t(kStop), tPred, symB);
end

%% 3C -- integrator stability ----------------------------------------------
fprintf('\n----------------------------------------------------------------\n');
fprintf(' 3C  integrator stability: linear prediction vs SEEDED runs\n');
fprintf('----------------------------------------------------------------\n');
fprintf('    Each run starts from the 3A equilibrium with a 1e-4 m/s lateral velocity\n');
fprintf('    seed. A stable integration of this dissipative system must shed it; an\n');
fprintf('    unstable one amplifies it until friction saturates (bounded chatter).\n');
fprintf('    Criterion: observed (seed grew?) must match predicted (max|R| > 1?).\n\n');

SEED = 1e-4;
vUH = P.v_bl;                          % UH-60 value, 0.1 ft/s
cfg = struct('vbl', {vUH, vUH, 0.01, 0.01, 0.15}, ...
             'dt',  {0.001, 0.01, 0.001, 0.0005, 0.01}, ...
             'label', {'physics ref (UH-60 v_bl)', ...
                       'control run B (UH-60)', ...
                       'old tight v_bl @ 0.001', ...
                       'old tight v_bl @ 0.0005', ...
                       'loosened until stable'});
res = struct('vbl',{},'dt',{},'label',{},'lamMax',{},'rho',{},'predStable',{}, ...
             'growth',{},'obsStable',{},'creep',{},'creepPred',{},'t',{},'vAlong',{},'vLat',{});

fprintf('    %-26s %5s %7s | %9s %9s %6s | %8s %6s | %11s %11s | %s\n', ...
        'case','v_bl','dt','lam_max','max|R|','pred','growth','obs','creep','analytic','match');
fprintf('    %s\n', repmat('-',1,126));
allMatch = true;
for c = 1:numel(cfg)
    Pc = P;  Pc.v_bl = cfg(c).vbl;
    [~, lamC] = linearize_eom(xS, Pc);          % valid: foot speeds < every v_bl here
    rho  = rhoOf(lamC, cfg(c).dt);
    pred = rho < 1;

    xs = xS;  xs(2) = xs(2) + SEED;
    Sc = run_landing_sim(xs, 2.0, cfg(c).dt, Pc);
    Vtr = inertialVel(Sc);
    late = Sc.t >= 1.5;
    vLat = Vtr*eL;   vAl = -(Vtr*tu);
    growth = max(abs(vLat(late)))/SEED;
    obs    = growth < 1;
    creep  = mean(vAl(late));
    cPred  = cfg(c).vbl*tan(a)/mu;
    match  = (pred == obs) && all(isfinite(Sc.x(:)));
    allMatch = allMatch && match;

    res(c) = struct('vbl',cfg(c).vbl,'dt',cfg(c).dt,'label',cfg(c).label, ...
        'lamMax',max(abs(lamC)),'rho',rho,'predStable',pred,'growth',growth, ...
        'obsStable',obs,'creep',creep,'creepPred',cPred,'t',Sc.t,'vAlong',vAl,'vLat',vLat);
    fprintf('    %-26s %5.2f %7.4f | %9.1f %9.3g %6s | %8.3g %6s | %11.4e %11.4e | %s\n', ...
        cfg(c).label, cfg(c).vbl, cfg(c).dt, max(abs(lamC)), rho, stab(pred), ...
        growth, stab(obs), creep, cPred, tern(match,'PASS','FAIL'));
end
step3_pass = step3_pass && allMatch;

% The trap, demonstrated: the old unstable configuration WITHOUT a seed
Pc = P;  Pc.v_bl = 0.01;
Sun = run_landing_sim(xS, 2.0, 0.001, Pc);
Vun = inertialVel(Sun);  lateU = Sun.t >= 1.5;
fprintf('\n    Old tight v_bl (0.01) at dt 0.001 with NO seed:\n');
fprintf('      creep %.6e vs %.6e analytic, max lateral velocity %.1e\n', ...
        mean(-(Vun(lateU,:)*tu)), Pc.v_bl*tan(a)/mu, max(abs(Vun*eL)));
fprintf('      -> lateral motion stays EXACTLY zero: symmetry never excites the unstable\n');
fprintf('         mode, so an unseeded test cannot detect the instability. (Creep is still\n');
fprintf('         converging from the v_bl = 0.0305 starting equilibrium.)\n');

% Stability boundary in v_bl for each dt (bisection on the linear prediction)
fprintf('\n    Smallest stable v_bl (bisection on max|R| = 1), and the creep it implies:\n');
vbStar = zeros(1,2);  dts = [0.001 0.01];
for q = 1:2
    lo = 0.006;  hi = 1.0;                % lo unstable, hi stable at both dts;
                                          % lo > creep speed keeps linearization valid
    for it = 1:40
        mid = 0.5*(lo+hi);  Pm = P;  Pm.v_bl = mid;
        [~, lm] = linearize_eom(xS, Pm);
        if rhoOf(lm, dts(q)) < 1, hi = mid; else, lo = mid; end
    end
    vbStar(q) = hi;
    fprintf('      dt = %-6g  v_bl >= %.4f m/s  ->  creep >= %.2f mm/s\n', ...
            dts(q), hi, 1000*hi*tan(a)/mu);
end

%% Comparison figure for 3C ------------------------------------------------
outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
fig = figure('Visible','off','Position',[100 100 1000 760]);
cols = lines(numel(res));
subplot(2,1,1); hold on; grid on;
for c = 1:numel(res)
    plot(res(c).t, abs(res(c).vAlong), 'Color', cols(c,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('%s (v_{bl} %.2f, dt %g)', res(c).label, res(c).vbl, res(c).dt));
    yline(res(c).creepPred, '--', 'Color', cols(c,:), 'HandleVisibility','off');
end
set(gca,'YScale','log'); ylabel('downhill creep speed [m/s]');
title('solid: simulated   dashed: analytic v_{bl} tan\alpha / \mu', 'FontWeight','normal');
legend('Location','eastoutside');
subplot(2,1,2); hold on; grid on;
for c = 1:numel(res)
    plot(res(c).t, abs(res(c).vLat), 'Color', cols(c,:), 'LineWidth', 1.2);
end
yline(SEED, 'k--', 'seed', 'LabelHorizontalAlignment','left');
set(gca,'YScale','log'); ylabel('|lateral CG velocity| [m/s]'); xlabel('time [s]');
title('stable integration sheds the seed; unstable integration amplifies it', 'FontWeight','normal');
sgtitle(sprintf(['Test 3C: integrator stability on a %g deg slope, seeded with a %.0e m/s lateral velocity\n' ...
                 'm = %g kg, k = %.3g N/m, c = %g N s/m, mu = %g'], SLOPE_DEG, SEED, m, P.k, P.c, mu), ...
        'Interpreter','tex', 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, 'test3C_integrator_stability_seeded_runs.png'), 'Resolution', 150);
close(fig);
fprintf('\n  figure written: test3C_integrator_stability_seeded_runs.png\n');

%% Summary ------------------------------------------------------------------
fprintf('\n================================================================\n');
if step3_pass
    fprintf(' STEP 3: PASS -- friction statics, sliding, and stability all verified\n');
else
    fprintf(' STEP 3: FAIL\n');
end
fprintf('================================================================\n');

save('step3_result.mat', 'step3_pass', 'S3A', 'S3B', 'res', 'vbStar', 'xS', 'SLOPE_DEG');

%% local functions ------------------------------------------------------------
function Vi = inertialVel(S)
% Inertial CG velocity history, (n+1) x 3
nS = size(S.x, 1);
Vi = zeros(nS, 3);
for k = 1:nS
    xx = S.x(k,:);
    Vi(k,:) = ( rot_b2i(xx(7), xx(8), xx(9)) * xx(1:3).' ).';
end
end

function s = stab(b)
if b, s = 'stable'; else, s = 'UNSTAB'; end
end

function out = tern(c,a,b)
if c, out = a; else, out = b; end
end
