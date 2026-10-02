function F_t = friction_physics(v_t, N, P)
%FRICTION_PHYSICS  Regularized Coulomb friction at one contact patch.
%
%   F_t = friction_physics(v_t, N, P)
%
%   Inputs:
%     v_t : 3x1 tangential velocity of the foot, INERTIAL frame, already
%           projected onto the terrain tangent plane [m/s]
%     N   : scalar normal force magnitude [N], already clamped >= 0
%     P   : params struct (uses P.mu, P.v_bl)
%   Output:
%     F_t : 3x1 friction force on the AIRCRAFT, inertial frame, lying in the
%           tangent plane (parallel to v_t) [N]
%
%   Magnitude saturates at mu*N and opposes v_t. Below |v_t| = v_bl the law
%   becomes linear-viscous, which is what makes it integrable -- and what sets
%   the stiffness (c_eff = mu*N/v_bl) that drives the timestep requirement.
%   See PLAN.md Sec.6.
%
%   CONSEQUENCE ON A SLOPE: regularized friction cannot hold perfectly still.
%   To resist the downhill gravity component it needs a nonzero sliding speed,
%   so the vehicle CREEPS downhill at exactly
%         v_creep = v_bl * tan(alpha) / mu
%   This is the regularization artifact the study is about; step 3 tests it.
%
%   This signature is the SWAPPABLE INTERFACE. friction_ml.m must match it
%   exactly so that changing models is a one-line change.

if P.mu == 0 || N <= 0
    F_t = zeros(size(v_t));
    return
end

speed = norm(v_t);
F_t   = -P.mu * N * v_t / max(speed, P.v_bl);
end
