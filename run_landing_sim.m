function S = run_landing_sim(x0, tEnd, dt, P, frictionFcn)
%RUN_LANDING_SIM  Fixed-step RK4 driver with full time-history logging.
%
%   S = run_landing_sim(x0, tEnd, dt, P, frictionFcn)
%
%   Inputs:
%     x0          12x1 initial state (see landing_eom.m for layout)
%     tEnd        simulation end time [s]
%     dt          fixed timestep [s]
%     P           params struct from landing_params()
%     frictionFcn function handle, default @friction_physics
%
%   Output struct S:
%     S.t        (n+1)x1        time [s]
%     S.x        (n+1)x12       full state history
%     S.alt      (n+1)x1        altitude = -zD [m]  (convenience)
%     S.N        (n+1)xnLegs    per-leg normal force [N]
%     S.delta    (n+1)xnLegs    per-leg penetration [m]
%     S.contact  (n+1)xnLegs    per-leg contact flag
%     S.Ft       (n+1)x3xnLegs  per-leg friction force, INERTIAL [N]
%     S.E        (n+1)x1        total mechanical energy [J]
%     S.dt, S.P                 run configuration
%
%   aux quantities are logged by evaluating the EOM AT the logged state, so
%   they correspond exactly to S.x(n,:) rather than to an RK4 interior stage.

if nargin < 5 || isempty(frictionFcn)
    frictionFcn = @friction_physics;
end

x0 = x0(:);
if numel(x0) ~= 12
    error('run_landing_sim:badState', 'x0 must be a 12-element state vector.');
end

nSteps = round(tEnd/dt);
if abs(nSteps*dt - tEnd) > 1e-12
    warning('run_landing_sim:dtMismatch', ...
            'tEnd=%g is not an integer multiple of dt=%g; using %d steps (tEnd=%g).', ...
            tEnd, dt, nSteps, nSteps*dt);
end

S.t       = zeros(nSteps+1, 1);
S.x       = zeros(nSteps+1, 12);
S.N       = zeros(nSteps+1, P.nLegs);
S.delta   = zeros(nSteps+1, P.nLegs);
S.contact = false(nSteps+1, P.nLegs);
S.Ft      = zeros(nSteps+1, 3, P.nLegs);    % per-leg friction, INERTIAL [N]
S.E       = zeros(nSteps+1, 1);

f = @(tt, xx) landing_eom(tt, xx, P, frictionFcn);

x = x0;
for n = 0:nSteps
    t = n*dt;

    [~, aux] = landing_eom(t, x, P, frictionFcn);
    Ecur     = landing_energy(x, P);

    S.t(n+1)         = t;
    S.x(n+1,:)       = x.';
    S.N(n+1,:)       = aux.N;
    S.delta(n+1,:)   = aux.delta;
    S.contact(n+1,:) = aux.contact;
    S.Ft(n+1,:,:)    = reshape(aux.F_t, [1 3 P.nLegs]);
    S.E(n+1)         = Ecur.total;

    if n < nSteps
        x = rk4_step(f, t, x, dt);
    end
end

S.alt = -S.x(:,12);
S.dt  = dt;
S.P   = P;
end
