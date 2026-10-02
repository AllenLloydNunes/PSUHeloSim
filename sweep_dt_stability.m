%% sweep_dt_stability.m
% Find the largest stable timestep for the CURRENT model (mu = 0.6,
% v_bl = 0.1 ft/s = 0.03048 m/s) holding on a 5 deg incline, by increasing dt
% in small steps and running the model at each one.
%
% Two sweeps, because the answer depends on what gets excited:
%   LATERAL seed      1e-4 m/s sideways velocity. Excites the lateral/roll
%                     friction mode -- the stiffest one. This is the TRUE limit.
%   LONGITUDINAL seed 1e-4 m/s along-slope velocity only. Keeps the run exactly
%                     laterally symmetric, so the lateral mode is never excited.
%                     This is the limit a symmetric test would (wrongly) report.
%
% Classification (same physical criterion as test_step3 3C): a stable
% integration of this dissipative system must shed the seed. Growth = largest
% deviation over the last 0.5 s divided by the seed; stable iff growth < 1.
%
% Each sweep: coarse pass over a dt grid, then bisection by SIMULATION between
% the last stable and first unstable dt. Both are compared against the linear
% prediction (RK4 amplification |R(h*lambda)| = 1) from linearize_eom.
%
% Run from SimpleLandingModel:  matlab -batch "sweep_dt_stability"

clear; clc;
fprintf('================================================================\n');
fprintf(' dt stability sweep -- current model on a 5 deg incline\n');
fprintf('================================================================\n');

P = landing_params();
SLOPE_DEG = 5;
P.terrain = make_terrain(SLOPE_DEG, 0);
a  = SLOPE_DEG*pi/180;
n  = P.terrain.n_down;
tu = P.terrain.t_up;
eL = cross(n, tu);
v_creep = P.v_bl*tan(a)/P.mu;
fprintf('  mu = %.2f, v_bl = %.5f m/s, analytic creep = %.4f mm/s\n', P.mu, P.v_bl, 1000*v_creep);

%% Equilibrium on the slope (settled at the known-stable dt = 0.001)
x0 = zeros(12,1);  x0(8) = a;  x0(10:12) = -1.75*n;
S0 = run_landing_sim(x0, 10.0, 0.001, P);
xS = S0.x(end,:).';

%% Linear prediction, split by symmetry
% At a laterally symmetric equilibrium the Jacobian decouples into a
% longitudinal block and a lateral block.
iLon = [1 3 5 8 10 12];          % u w q theta xN zD
iLat = [2 4 6 7 9 11];           % v p r phi psi yE
[J, lamAll] = linearize_eom(xS, P);
coupling = max(norm(J(iLon,iLat)), norm(J(iLat,iLon)));
lamLon = eig(J(iLon,iLon));
lamLat = eig(J(iLat,iLat));

Rk    = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;
rhoOf = @(lam, dt) max(abs(Rk(dt*lam(abs(lam) > 1))));
dtLim = @(lam) fzero(@(h) rhoOf(lam,h) - 1, [1e-4 0.05]);

dtLinAll = dtLim(lamAll);
dtLinLat = dtLim(lamLat);
dtLinLon = dtLim(lamLon);
fprintf('\n  Linear prediction (RK4, |R(h*lambda)| = 1):\n');
fprintf('    block coupling at symmetric equilibrium: %.2e (should be ~0)\n', coupling);
fprintf('    lateral/roll block      fastest |lambda| = %7.1f rad/s  ->  dt_max = %.4f ms\n', ...
        max(abs(lamLat)), 1000*dtLinLat);
fprintf('    longitudinal/pitch block fastest |lambda| = %7.1f rad/s  ->  dt_max = %.4f ms\n', ...
        max(abs(lamLon)), 1000*dtLinLon);
fprintf('    full model                                        ->  dt_max = %.4f ms\n', 1000*dtLinAll);

%% Sweeps
SEED  = 1e-4;
dtGrid = [0.0010:0.0001:0.0040, 0.0045:0.0005:0.0120];
modes = struct('name', {'LATERAL seed (true limit)', 'LONGITUDINAL seed (symmetric test)'}, ...
               'dir',  {eL, tu}, 'dtLin', {dtLinLat, dtLinLon});
sweep = struct();

for mIdx = 1:2
    fprintf('\n----------------------------------------------------------------\n');
    fprintf(' %s\n', modes(mIdx).name);
    fprintf('----------------------------------------------------------------\n');
    fprintf('    %9s | %10s | %9s | %12s\n', 'dt [ms]', 'growth', 'result', 'creep [mm/s]');
    g = nan(size(dtGrid));  cr = nan(size(dtGrid));
    for k = 1:numel(dtGrid)
        [g(k), cr(k)] = runCase(xS, P, dtGrid(k), SEED, modes(mIdx).dir, tu, eL, v_creep, mIdx);
        fprintf('    %9.2f | %10.3g | %9s | %12.4f\n', 1000*dtGrid(k), g(k), ...
                tern(g(k) < 1, 'stable', 'UNSTABLE'), 1000*cr(k));
    end

    % Bracket the FIRST instability on the grid, then bisect by simulation
    kU = find(~(g < 1), 1);
    if isempty(kU) || kU == 1
        error('sweep:bracket', 'Could not bracket the stability boundary on the grid.');
    end
    lo = dtGrid(kU-1);  hi = dtGrid(kU);
    for it = 1:14
        mid = 0.5*(lo + hi);
        gm = runCase(xS, P, mid, SEED, modes(mIdx).dir, tu, eL, v_creep, mIdx, 3.0);
        if gm < 1, lo = mid; else, hi = mid; end
    end
    dtSim = 0.5*(lo + hi);
    errPct = 100*(dtSim - modes(mIdx).dtLin)/modes(mIdx).dtLin;
    fprintf('\n    simulated boundary : dt = %.4f ms  (bisection bracket +/- %.1e us; the\n', 1000*dtSim, 1e6*(hi-lo)/2);
    fprintf('                         growth criterion itself resolves ~1 us near the boundary)\n');
    fprintf('    linear prediction  : dt = %.4f ms  (difference %+.2f%%)\n', 1000*modes(mIdx).dtLin, errPct);

    sweep(mIdx).name = modes(mIdx).name;
    sweep(mIdx).growth = g;  sweep(mIdx).creep = cr;
    sweep(mIdx).dtSim = dtSim;  sweep(mIdx).dtLin = modes(mIdx).dtLin;
end

%% Summary
fprintf('\n================================================================\n');
fprintf(' SUMMARY\n');
fprintf('================================================================\n');
fprintf('  Largest stable dt, real (laterally seeded) : %.3f ms   (linear %.3f ms)\n', ...
        1000*sweep(1).dtSim, 1000*sweep(1).dtLin);
fprintf('  Largest stable dt, symmetric test only     : %.3f ms   (linear %.3f ms)\n', ...
        1000*sweep(2).dtSim, 1000*sweep(2).dtLin);
fprintf('  Margin at the physics reference dt = 1 ms  : %.2fx\n', sweep(1).dtSim/0.001);
fprintf('  Target ML dt = 10 ms is %.1fx beyond the real limit\n', 0.01/sweep(1).dtSim);

%% Figure
outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
fig = figure('Visible','off','Position',[100 100 1050 800]);
cols = [0.00 0.45 0.74; 0.85 0.33 0.10];
subplot(2,1,1); hold on; grid on;
for mIdx = 1:2
    semilogy(1000*dtGrid, sweep(mIdx).growth, 'o-', 'Color', cols(mIdx,:), ...
             'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', sweep(mIdx).name);
    xline(1000*sweep(mIdx).dtLin, '--', 'Color', cols(mIdx,:), 'LineWidth', 1.0, ...
          'HandleVisibility','off');
end
set(gca, 'YScale', 'log');
yline(1, 'k-', 'growth = 1 (stability boundary)', 'LabelHorizontalAlignment','left', ...
      'HandleVisibility','off');
xline(1, 'k:', 'physics ref 1 ms', 'LabelOrientation','horizontal', 'HandleVisibility','off');
xline(10, 'k:', 'ML target 10 ms', 'LabelOrientation','horizontal', ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
ylabel('seed growth after 2 s'); legend('Location','southeast');
title('solid: simulated   dashed: linear prediction |R(h\lambda)| = 1', 'FontWeight','normal');
subplot(2,1,2); hold on; grid on;
for mIdx = 1:2
    plot(1000*dtGrid, 1000*sweep(mIdx).creep, 'o-', 'Color', cols(mIdx,:), ...
         'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName', sweep(mIdx).name);
end
yline(1000*v_creep, 'k--', sprintf('analytic creep %.2f mm/s', 1000*v_creep), ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xlabel('timestep dt [ms]'); ylabel('mean downhill creep [mm/s]'); legend('Location','northwest');
title('creep is exact while stable; wrong once the integration goes unstable', 'FontWeight','normal');
sgtitle(sprintf(['Timestep sweep: largest stable dt on a %g deg incline, RK4\n' ...
                 'mu = %g, v_{bl} = %.5f m/s (0.1 ft/s), m = %g kg, k = %.3g N/m, c = %g N s/m'], ...
                SLOPE_DEG, P.mu, P.v_bl, P.m, P.k, P.c), 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, 'sweep_dt_stability.png'), 'Resolution', 150);
close(fig);
fprintf('\n  figure written: figures/sweep_dt_stability.png\n');

save('sweep_dt_result.mat', 'sweep', 'dtGrid', 'lamLat', 'lamLon', 'lamAll', 'xS', 'P');

%% local functions -------------------------------------------------------------
function [growth, creep] = runCase(xS, P, dt, seed, dirSeed, tu, eL, v_creep, mode, tEnd)
% Seeded run from the slope equilibrium. Mode 1 measures lateral velocity;
% mode 2 measures along-slope velocity minus the analytic creep.
if nargin < 10, tEnd = 2.0; end
nSteps = round(tEnd/dt);
tEnd   = nSteps*dt;                      % exact multiple of dt -> no warning
R0 = rot_b2i(xS(7), xS(8), xS(9));
x0 = xS;
x0(1:3) = x0(1:3) + R0.' * (seed*dirSeed);
S = run_landing_sim(x0, tEnd, dt, P);
late = S.t >= tEnd - 0.5;
Vi = zeros(nnz(late), 3);  idx = find(late);
for k = 1:numel(idx)
    xx = S.x(idx(k),:);
    Vi(k,:) = ( rot_b2i(xx(7), xx(8), xx(9)) * xx(1:3).' ).';
end
if ~all(isfinite(S.x(:)))
    growth = Inf;  creep = NaN;  return
end
vAlong = -(Vi*tu);
creep  = mean(vAlong);
if mode == 1
    growth = max(abs(Vi*eL))/seed;
else
    growth = max(abs(vAlong - v_creep))/seed;
end
end

function out = tern(c,a,b)
if c, out = a; else, out = b; end
end
