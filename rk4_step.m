function x_new = rk4_step(f, t, x, dt)
%RK4_STEP  One fixed-step classical Runge-Kutta 4 step.
%
%   x_new = rk4_step(f, t, x, dt)   with  f = @(t,x) xdot
%
%   Hand-written rather than ode45 on purpose: this study is ABOUT the
%   timestep, so dt must be explicit and fixed, not solver-chosen.
%
%   Note: RK4 integrates constant-acceleration dynamics EXACTLY, which is what
%   makes the step-1 free-fall test a machine-precision check rather than an
%   engineering-tolerance one.

k1 = f(t,        x);
k2 = f(t + dt/2, x + (dt/2)*k1);
k3 = f(t + dt/2, x + (dt/2)*k2);
k4 = f(t + dt,   x +  dt   *k3);

x_new = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end
