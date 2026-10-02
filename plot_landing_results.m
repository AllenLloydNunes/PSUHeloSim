function plot_landing_results(S, tag, desc, outDir)
%PLOT_LANDING_RESULTS  Standard time-history plots, saved to PNG.
%
%   plot_landing_results(S, tag, desc, outDir)
%
%     S      run struct from run_landing_sim
%     tag    short filename-safe identifier, e.g. 'test2A_damped_drop'
%     desc   human-readable description used as the figure title, e.g.
%            'Test 2A: damped drop from 0.2 m, settles to static equilibrium'
%            (optional -- defaults to tag with underscores replaced by spaces)
%     outDir output directory (optional, defaults to ./figures)
%
%   Figures are created invisible and written to disk so this works under
%   'matlab -batch' where no display is available.
%
%   Produces:
%     <tag>_position.png     altitude, north, east position vs t
%     <tag>_attitude.png     roll, pitch, yaw vs t
%     <tag>_gear_loads.png   per-leg normal force and penetration vs t
%     <tag>_energy.png       energy components and total drift vs t

if nargin < 2
    error('plot_landing_results:notEnoughArgs', ...
        ['plot_landing_results is a FUNCTION -- it needs a run struct to plot.\n' ...
         'Typing the bare name calls it with no arguments, which is this error.\n\n' ...
         'To regenerate every standard figure:\n' ...
         '    run_all_tests\n\n' ...
         'To plot a saved run from a previous test:\n' ...
         '    r = load(''step2_result.mat'');   %% S_settle, S_econt, S_ebounce\n' ...
         '    plot_landing_results(r.S_settle, ''my_tag'', ''my description'')\n\n' ...
         'To run your own case and plot it:\n' ...
         '    P  = landing_params();\n' ...
         '    x0 = zeros(12,1);  x0(12) = -2.0;\n' ...
         '    S  = run_landing_sim(x0, 5.0, P.dt_physics, P);\n' ...
         '    plot_landing_results(S, ''my_drop'', ''My drop test'')']);
end

if nargin < 3 || isempty(desc)
    desc = strrep(tag, '_', ' ');
end
if nargin < 4 || isempty(outDir)
    outDir = fullfile(fileparts(mfilename('fullpath')), 'figures');
end
if ~exist(outDir, 'dir'), mkdir(outDir); end

t = S.t;
P = S.P;

% Run conditions, appended to every title so a stray PNG is self-describing
cond = sprintf('m = %g kg, k = %.3g N/m, c = %g N s/m, mu = %g, v_bl = %g m/s, slope = %g deg, dt = %g s', ...
               P.m, P.k, P.c, P.mu, P.v_bl, P.terrain.slope_deg, S.dt);

%% Position -------------------------------------------------------------
fig = figure('Visible','off','Position',[100 100 1050 720]);
subplot(3,1,1);
plot(t, S.alt, 'LineWidth', 1.2); grid on;
ylabel('altitude  -z_D  [m]');
subplot(3,1,2);
plot(t, S.x(:,10), 'LineWidth', 1.2); grid on; ylabel('north  x_N  [m]');
subplot(3,1,3);
plot(t, S.x(:,11), 'LineWidth', 1.2); grid on;
ylabel('east  y_E  [m]'); xlabel('time [s]');
sgtitle(sprintf('%s\nCG position vs time\n%s', desc, cond), ...
        'Interpreter','none', 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, [tag '_position.png']), 'Resolution', 150);
close(fig);

%% Attitude -------------------------------------------------------------
fig = figure('Visible','off','Position',[100 100 1050 720]);
subplot(3,1,1);
plot(t, S.x(:,7)*180/pi, 'LineWidth', 1.2); grid on;
ylabel('roll  \phi  [deg]');
subplot(3,1,2);
plot(t, S.x(:,8)*180/pi, 'LineWidth', 1.2); grid on;
ylabel('pitch  \theta  [deg]');
subplot(3,1,3);
plot(t, S.x(:,9)*180/pi, 'LineWidth', 1.2); grid on;
ylabel('yaw  \psi  [deg]'); xlabel('time [s]');
sgtitle(sprintf('%s\nAttitude vs time\n%s', desc, cond), ...
        'Interpreter','none', 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, [tag '_attitude.png']), 'Resolution', 150);
close(fig);

%% Gear loads -----------------------------------------------------------
%  With friction active a third panel shows per-leg friction magnitude against
%  the Coulomb limit mu*N (dashed): below the dashed line = holding/creeping
%  inside the boundary layer, on it = fully sliding.
hasFric = P.mu > 0 && isfield(S,'Ft');
nPan = 2 + hasFric;
fig = figure('Visible','off','Position',[100 100 1050 300+240*nPan]);
subplot(nPan,1,1);
plot(t, S.N, 'LineWidth', 1.2); grid on;
ylabel('normal force  N  [N]'); legend(P.legNames, 'Location','best');
if hasFric
    subplot(nPan,1,2);
    FtMag = squeeze(sqrt(sum(S.Ft.^2, 2)));
    hF = plot(t, FtMag, 'LineWidth', 1.2); grid on; hold on;
    for i = 1:P.nLegs
        plot(t, P.mu*S.N(:,i), '--', 'Color', hF(i).Color, 'LineWidth', 0.8);
    end
    ylabel('friction  |F_t|  [N]');
    legend(hF, P.legNames, 'Location','best');
    title('solid: friction magnitude   dashed: Coulomb limit \mu N', 'FontWeight','normal');
end
subplot(nPan,1,nPan);
plot(t, S.delta*1000, 'LineWidth', 1.2); grid on;
yline(0, 'k--', 'ground plane', 'LabelHorizontalAlignment','left');
ylabel('penetration  \delta  [mm]'); xlabel('time [s]');
legend(P.legNames, 'Location','best');
sgtitle(sprintf('%s\nPer-leg normal force%s and ground penetration\n%s', desc, ...
        repmat(', friction', 1, double(hasFric)), cond), ...
        'Interpreter','none', 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, [tag '_gear_loads.png']), 'Resolution', 150);
close(fig);

%% Energy ---------------------------------------------------------------
nE = numel(t);
KEt = zeros(nE,1); KEr = zeros(nE,1); PEg = zeros(nE,1); PEs = zeros(nE,1);
for n = 1:nE
    E = landing_energy(S.x(n,:).', P);
    KEt(n) = E.KE_trans; KEr(n) = E.KE_rot;
    PEg(n) = E.PE_grav;  PEs(n) = E.PE_spring;
end

if P.c == 0 && P.mu == 0
    note = 'no dissipation (c = 0, mu = 0): total energy MUST be conserved';
else
    note = 'damper active: total energy is expected to DECREASE';
end

fig = figure('Visible','off','Position',[100 100 1050 660]);
subplot(2,1,1);
plot(t, KEt, t, KEr, t, PEg, t, PEs, 'LineWidth', 1.2); grid on;
legend({'KE translational','KE rotational','PE gravity','PE spring'}, 'Location','best');
ylabel('energy [J]');
subplot(2,1,2);
plot(t, S.E - S.E(1), 'LineWidth', 1.2); grid on;
ylabel('total E - E(0)  [J]'); xlabel('time [s]');
sgtitle(sprintf('%s\nEnergy budget -- %s\n%s', desc, note, cond), ...
        'Interpreter','none', 'FontWeight','bold', 'FontSize', 11);
exportgraphics(fig, fullfile(outDir, [tag '_energy.png']), 'Resolution', 150);
close(fig);

fprintf('  figures written: %s_{position,attitude,gear_loads,energy}.png\n', tag);
end
