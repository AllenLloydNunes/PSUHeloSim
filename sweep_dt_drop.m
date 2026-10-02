%% sweep_dt_drop.m
% Drop-test timestep sweep: find where a realistic landing goes unstable, and
% how accuracy degrades with dt before that.
%
% SCENARIO
%   5 deg incline rising toward SOUTH, so the vehicle, facing North, lands
%   NOSE-FIRST going DOWNHILL (gravity keeps pushing while friction stops it).
%   Released with the body aligned to the slope plus a 3 deg ROLL offset, so one
%   front wheel touches first. That asymmetry is essential: a laterally
%   symmetric drop never excites the lateral/roll friction mode and would report
%   a falsely large stable dt (~10.5 ms, see sweep_dt_stability.m).
%   Initial velocity 3 m/s along the slope; lowest foot ~0.12 m above ground.
%   Run 4 s: fall, touchdown, slide, stop, ~3 s of holding.
%
% REFERENCE
%   dt = 0.5 ms, checked against dt = 0.25 ms to measure its own error.
%
% Run from SimpleLandingModel:  matlab -batch "sweep_dt_drop"

clear; clc;
fprintf('================================================================\n');
fprintf(' Drop-test dt sweep: nose-first downhill landing, 5 deg, 3 deg roll\n');
fprintf('================================================================\n');

%% Scenario ----------------------------------------------------------------
P = landing_params();
SLOPE_DEG = 5;  ROLL0_DEG = 3;  V0 = 3.0;  TEND = 4.0;
P.terrain = make_terrain(SLOPE_DEG, 180);        % rises toward South
a  = SLOPE_DEG*pi/180;
n  = P.terrain.n_down;
dn = -P.terrain.t_up;                            % downhill unit vector (= North, down-slope)
v_creep = P.v_bl*tan(a)/P.mu;

x0 = zeros(12,1);
x0(7) = ROLL0_DEG*pi/180;                        % roll offset -> one wheel first
x0(8) = -a;                                      % nose down, aligned with slope
x0(10:12) = -2.0*n;                              % aligned feet 0.2 m above ground
R0 = rot_b2i(x0(7), x0(8), x0(9));
x0(1:3) = R0.' * (V0*dn);                        % 3 m/s along the slope, downhill

fprintf('  mu = %.2f, v_bl = %.5f m/s, analytic hold creep = %.4f mm/s\n', P.mu, P.v_bl, 1000*v_creep);

%% Reference runs ----------------------------------------------------------
REF  = runDrop(x0, P, 0.0005, TEND, n, dn);
REF2 = runDrop(x0, P, 0.00025, TEND, n, dn);
fprintf('\n  Reference (dt 0.5 ms) vs half-step (0.25 ms) -- the reference''s own error:\n');
refErr = compareRuns(REF, REF2);
printMetrics('0.5 ms', REF);
printMetrics('0.25 ms', REF2);
fprintf('    differences: stop dist %.2e m, peak N %.2e N, final pos %.2e m, final att %.2e deg\n', ...
        refErr.dStop, refErr.dPeakN, refErr.dFinalPos, refErr.dFinalAtt);

%% Linear prediction at the landed, held state -----------------------------
xF = REF.S.x(end,:).';
RF = rot_b2i(xF(7), xF(8), xF(9));
[~, auxF] = landing_eom(0, xF, P);
footSpd = zeros(1,P.nLegs);
for i = 1:P.nLegs
    Vf = RF*( xF(1:3) + cross(xF(4:6), P.p_foot_b(:,i)) );
    footSpd(i) = norm(Vf - (n.'*Vf)*n);
end
if ~(all(auxF.contact) && all(auxF.N > 0) && all(footSpd < P.v_bl))
    error('sweep_dt_drop:lin', 'Final reference state is not a valid linearization point.');
end
[~, lamF] = linearize_eom(xF, P);
Rk    = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;
rhoOf = @(lam, dt) max(abs(Rk(dt*lam(abs(lam) > 1))));
dtLin = fzero(@(h) rhoOf(lamF,h) - 1, [1e-4 0.05]);
fprintf('\n  Linearized at the landed state: fastest |lambda| = %.1f rad/s -> predicted dt_max = %.4f ms\n', ...
        max(abs(lamF)), 1000*dtLin);
fprintf('  (final heading %.2f deg, roll %.3f deg -- the landing left it %s)\n', ...
        xF(9)*180/pi, xF(7)*180/pi, tern(abs(xF(9)) > 1e-6, 'yawed, i.e. asymmetric', 'symmetric'));

%% Sweep -------------------------------------------------------------------
dtGrid = [0.0010:0.0001:0.0040, 0.0045:0.0005:0.0120];
fprintf('\n  %7s | %9s %9s %9s %9s %9s | %10s %10s %8s\n', 'dt[ms]', 'dStop[m]', 'dTstop[s]', ...
        'dPeakN[%]', 'dPos[m]', 'dAtt[deg]', 'wLate/ref', 'decay', 'creep');
fprintf('  %s\n', repmat('-',1,104));
R = struct('dt',{},'dStop',{},'dTstop',{},'dPeakNpct',{},'dFinalPos',{},'dFinalAtt',{}, ...
           'wRatio',{},'decay',{},'creep',{},'finite',{});
for k = 1:numel(dtGrid)
    D = runDrop(x0, P, dtGrid(k), TEND, n, dn);
    e = compareRuns(REF, D);
    R(k).dt = dtGrid(k);           R(k).dStop = e.dStop;       R(k).dTstop = e.dTstop;
    R(k).dPeakNpct = e.dPeakNpct;  R(k).dFinalPos = e.dFinalPos;
    R(k).dFinalAtt = e.dFinalAtt;  R(k).wRatio = D.wLate/REF.wLate;
    R(k).decay = D.decay;          R(k).creep = D.creep;       R(k).finite = D.finite;
    fprintf('  %7.2f | %9.2e %9.2e %9.2f %9.2e %9.2e | %10.3g %10.3g %8.3f\n', ...
        1000*dtGrid(k), e.dStop, e.dTstop, e.dPeakNpct, e.dFinalPos, e.dFinalAtt, ...
        R(k).wRatio, R(k).decay, 1000*R(k).creep);
end
fprintf('\n  (reference late |omega| = %.3e rad/s, decay = %.3f, creep = %.4f mm/s)\n', ...
        REF.wLate, REF.decay, 1000*REF.creep);

%% Is the reference's late creep converged? -------------------------------
% After 4 s the landed vehicle is still rotating slowly (the landing left it
% yawed), so its creep is not yet the analytic steady value. Run the
% reference hold much longer to confirm it converges there.
LONG = runDrop(x0, P, 0.0005, 20.0, n, dn);
fprintf('\n  Creep convergence of the reference (dt 0.5 ms):\n');
fprintf('    at 4 s : %.4f mm/s   late |omega| %.2e rad/s\n', 1000*REF.creep,  REF.wLate);
fprintf('    at 20 s: %.4f mm/s   late |omega| %.2e rad/s   analytic %.4f mm/s  (err %.2e)\n', ...
        1000*LONG.creep, LONG.wLate, 1000*v_creep, abs(LONG.creep - v_creep)/v_creep);

%% Classification ------------------------------------------------------------
% Stable integration: the late residual motion is the physical landing
% transient dying out at its natural rate -- same level as the reference
% (wRatio ~ 1) and still shrinking (decay ~ 0.51 per 0.5 s, identical to the
% reference). Unstable integration: a sustained limit cycle (decay ~ 1) at
% several times the reference level. The thresholds below sit between those
% two clusters as observed in the coarse sweep; they are CROSS-CHECKED against
% the independent linear prediction rather than tuned to it.
isStable = @(D) D.finite && (D.decay < 0.75) && (D.wLate/REF.wLate < 2);
stab = arrayfun(@(r) r.finite && r.decay < 0.75 && r.wRatio < 2, R);
kU = find(~stab, 1);
if isempty(kU) || kU == 1 || any(~stab(1:kU-1))
    error('sweep_dt_drop:bracket', 'Stable/unstable split is not a single clean boundary on the grid.');
end
fprintf('\n  Coarse grid: stable up to %.2f ms, unstable from %.2f ms; %d/%d beyond also unstable\n', ...
        1000*dtGrid(kU-1), 1000*dtGrid(kU), nnz(~stab(kU:end)), numel(stab)-kU+1);

lo = dtGrid(kU-1);  hi = dtGrid(kU);
for it = 1:12
    mid = 0.5*(lo+hi);
    if isStable(runDrop(x0, P, mid, TEND, n, dn)), lo = mid; else, hi = mid; end
end
dtSim = 0.5*(lo+hi);
fprintf('  Bisection by simulation: boundary at dt = %.4f ms (bracket +/- %.1e us)\n', 1000*dtSim, 1e6*(hi-lo)/2);
fprintf('  Linear prediction      : dt = %.4f ms  (difference %+.2f%%)\n', 1000*dtLin, 100*(dtSim-dtLin)/dtLin);

%% Summary -------------------------------------------------------------------
st = R(stab);  un = R(~stab);
fprintf('\n================================================================\n');
fprintf(' SUMMARY -- drop test, nose-first downhill, 5 deg, 3 deg roll\n');
fprintf('================================================================\n');
fprintf('  Largest stable dt          : %.3f ms  (linear prediction %.3f ms)\n', 1000*dtSim, 1000*dtLin);
fprintf('  Accuracy while stable (dt <= %.1f ms), worst case vs 0.5 ms reference:\n', 1000*dtGrid(kU-1));
fprintf('    stop distance %.2e m of %.3f m | peak load %.2f%% | final position %.2e m | attitude %.2e deg\n', ...
        max([st.dStop]), REF.stopDist, max([st.dPeakNpct]), max([st.dFinalPos]), max([st.dFinalAtt]));
fprintf('  Beyond the limit (up to %.0f ms): landing still accurate, hold is not:\n', 1000*dtGrid(end));
fprintf('    stop distance %.2e m | peak load %.2f%% | final position up to %.3f m | creep up to %.1f mm/s\n', ...
        max([un.dStop]), max([un.dPeakNpct]), max([un.dFinalPos]), 1000*max([un.creep]));
k10 = find(abs(dtGrid - 0.01) < 1e-12);
fprintf('  At dt = 10 ms (ML target): stop dist err %.2e m, peak load err %.2f%%, final pos err %.3f m, creep %.1f mm/s\n', ...
        R(k10).dStop, R(k10).dPeakNpct, R(k10).dFinalPos, 1000*R(k10).creep);

%% Figures -------------------------------------------------------------------
outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
ms = 1000*dtGrid;

% (1) stability + accuracy vs dt
fig = figure('Visible','off','Position',[100 100 1050 900]);
subplot(3,1,1); hold on; grid on;
semilogy(ms, [R.wRatio], 'o-', 'LineWidth', 1.3, 'MarkerSize', 4);
set(gca,'YScale','log');
yline(2, 'k--', 'threshold', 'LabelHorizontalAlignment','left');
xline(1000*dtLin, 'r--', sprintf('linear limit %.3f ms', 1000*dtLin), 'LabelOrientation','horizontal');
xline(1, 'k:', 'ref 1 ms', 'LabelOrientation','horizontal');
xline(10, 'k:', 'ML 10 ms', 'LabelOrientation','horizontal', 'LabelHorizontalAlignment','left');
ylabel({'late |\omega| relative','to 0.5 ms reference'});
title('stability: residual motion during the hold phase (last 0.5 s)', 'FontWeight','normal');
subplot(3,1,2); hold on; grid on;
semilogy(ms, max([R.dStop],1e-7), 'o-', 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName','stop distance error [m]');
semilogy(ms, max([R.dFinalPos],1e-7), 's-', 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName','final position error [m]');
set(gca,'YScale','log');
xline(1000*dtLin, 'r--', 'HandleVisibility','off');
ylabel('error vs reference [m]'); legend('Location','northwest');
title('the slide is accurate at every dt; the final position is not, once the hold goes unstable', 'FontWeight','normal');
subplot(3,1,3); hold on; grid on;
plot(ms, [R.dPeakNpct], 'o-', 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName','peak leg load error [%]');
plot(ms, 1000*[R.creep], 's-', 'LineWidth', 1.3, 'MarkerSize', 4, 'DisplayName','late creep [mm/s]');
yline(1000*REF.creep, 'k--', sprintf('reference creep %.2f mm/s', 1000*REF.creep), ...
      'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xline(1000*dtLin, 'r--', 'HandleVisibility','off');
xlabel('timestep dt [ms]'); legend('Location','northwest');
title('touchdown load is resolved at every dt; hold creep is not', 'FontWeight','normal');
sgtitle(sprintf(['Drop-test timestep sweep: nose-first downhill landing at %g m/s on a %g deg incline, %g deg roll\n' ...
                 'mu = %g, v_{bl} = %.5f m/s (0.1 ft/s), RK4, reference dt = 0.5 ms'], ...
                V0, SLOPE_DEG, ROLL0_DEG, P.mu, P.v_bl), 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, 'sweep_dt_drop_accuracy_stability.png'), 'Resolution', 150);
close(fig);

% (2) time histories at representative dt
showDt = [0.0005, 0.0024, 0.0025, 0.01];
cols = lines(numel(showDt));
fig = figure('Visible','off','Position',[100 100 1050 800]);
for s = 1:numel(showDt)
    D = runDrop(x0, P, showDt(s), TEND, n, dn);
    nS = size(D.S.x,1);  spd = zeros(nS,1);  w = zeros(nS,1);
    for k = 1:nS
        xx = D.S.x(k,:);  Vi = rot_b2i(xx(7),xx(8),xx(9))*xx(1:3).';
        spd(k) = norm(Vi - (n.'*Vi)*n);  w(k) = norm(xx(4:6));
    end
    lbl = sprintf('dt = %g ms', 1000*showDt(s));
    subplot(2,1,1); hold on; plot(D.S.t, spd, 'Color', cols(s,:), 'LineWidth', 1.2, 'DisplayName', lbl);
    subplot(2,1,2); hold on; plot(D.S.t, w,   'Color', cols(s,:), 'LineWidth', 1.2, 'DisplayName', lbl);
end
subplot(2,1,1); grid on; set(gca,'YScale','log'); ylim([1e-3 5]);
yline(P.v_bl, 'k--', 'v_{bl}', 'HandleVisibility','off');
ylabel('in-plane CG speed [m/s]'); legend('Location','northeast');
title('touchdown ~0.16 s, slide, stop ~0.56 s, then hold: all dt agree until the hold', 'FontWeight','normal');
subplot(2,1,2); grid on; set(gca,'YScale','log');
ylabel('|\omega| [rad/s]'); xlabel('time [s]'); legend('Location','northeast');
title('during the hold, stable dt lets the landing transient die out; unstable dt sustains chatter', 'FontWeight','normal');
sgtitle(sprintf('Drop test time histories either side of the stability limit (%.3f ms)', 1000*dtSim), ...
        'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, 'sweep_dt_drop_time_histories.png'), 'Resolution', 150);
close(fig);
fprintf('\n  figures written: sweep_dt_drop_accuracy_stability.png, sweep_dt_drop_time_histories.png\n');

save('sweep_dt_drop_result.mat', 'R', 'REF', 'REF2', 'LONG', 'dtGrid', 'lamF', 'dtLin', 'dtSim', 'x0', 'P');

%% local functions -----------------------------------------------------------
function D = runDrop(x0, P, dt, tEnd, n, dn)
nSteps = round(tEnd/dt);  tEnd = nSteps*dt;
S = run_landing_sim(x0, tEnd, dt, P);
D.S = S;  D.dt = dt;
D.finite = all(isfinite(S.x(:)));
nS = size(S.x,1);
Vi = zeros(nS,3);  w = zeros(nS,1);
for k = 1:nS
    xx = S.x(k,:);
    Vi(k,:) = ( rot_b2i(xx(7),xx(8),xx(9)) * xx(1:3).' ).';
    w(k) = norm(xx(4:6));
end
vt  = Vi - (Vi*n)*n.';                       % in-plane CG velocity
spd = sqrt(sum(vt.^2,2));
anyC = any(S.contact,2);
kTouch = find(anyC, 1);
kStop  = kTouch - 1 + find(spd(kTouch:end) < 0.1, 1);
D.tTouch = S.t(kTouch);
if isempty(kStop), D.tStop = NaN; D.stopPos = [NaN;NaN;NaN];
else, D.tStop = S.t(kStop); D.stopPos = S.x(kStop,10:12).'; end
D.stopDist = dn.' * (D.stopPos - S.x(1,10:12).');
D.peakN  = max(S.N(:));
D.finalPos = S.x(end,10:12).';
D.finalAtt = S.x(end,7:9).'*180/pi;
w1 = S.t >= tEnd-1.0 & S.t < tEnd-0.5;
w2 = S.t >= tEnd-0.5;
D.wLate = max(w(w2));
D.decay = max(w(w2)) / max(w(w1));           % < 1: disturbances still shrinking
D.creep = mean(Vi(S.t >= tEnd-1.0,:)*dn);
end

function e = compareRuns(A, B)
e.dStop     = abs(B.stopDist - A.stopDist);
e.dTstop    = abs(B.tStop - A.tStop);
e.dPeakN    = abs(B.peakN - A.peakN);
e.dPeakNpct = 100*e.dPeakN/A.peakN;
e.dFinalPos = norm(B.finalPos - A.finalPos);
e.dFinalAtt = max(abs(B.finalAtt - A.finalAtt));
end

function printMetrics(lbl, D)
fprintf('    %-8s touch %.4f s | stop %.4f s after %.4f m | peak N %.1f N | final att [%.3f %.3f %.3f] deg\n', ...
        lbl, D.tTouch, D.tStop, D.stopDist, D.peakN, D.finalAtt);
end

function out = tern(c,a,b)
if c, out = a; else, out = b; end
end
