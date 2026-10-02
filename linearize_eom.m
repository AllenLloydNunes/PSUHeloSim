function [J, lam] = linearize_eom(x, P, frictionFcn, h)
%LINEARIZE_EOM  Numerical Jacobian of landing_eom and its eigenvalues.
%
%   [J, lam] = linearize_eom(x, P, frictionFcn, h)
%
%   Central finite differences, step h (default 1e-7), on all 12 states.
%
%   VALIDITY: the model contains switches -- contact on/off, max(N,0), and the
%   friction boundary layer max(|v_t|, v_bl). The Jacobian is only meaningful
%   if x is away from all of them: every relevant leg firmly in contact with
%   N > 0, and every foot's tangential speed strictly inside the boundary layer
%   (|v_t| < v_bl). The caller must check this; test_step3 asserts it.
%
%   Used to predict fixed-step integrator stability: for RK4 the step is
%   linearly stable iff |R(h*lambda)| < 1 for every stiff eigenvalue, with
%       R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24.

if nargin < 3 || isempty(frictionFcn), frictionFcn = @friction_physics; end
if nargin < 4 || isempty(h),           h = 1e-7;                        end

f  = @(xx) landing_eom(0, xx, P, frictionFcn);
x  = x(:);
nx = numel(x);
J  = zeros(nx);
for j = 1:nx
    xp = x;  xp(j) = xp(j) + h;
    xm = x;  xm(j) = xm(j) - h;
    J(:,j) = ( f(xp) - f(xm) ) / (2*h);
end
lam = eig(J);
end
