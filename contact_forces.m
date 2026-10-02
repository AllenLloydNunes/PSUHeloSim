function [F_b, M_b, aux] = contact_forces(x, R_b2i, P, frictionFcn)
%CONTACT_FORCES  Per-leg ground contact -> total force/moment at the CG.
%
%   [F_b, M_b, aux] = contact_forces(x, R_b2i, P, frictionFcn)
%
%   Outputs F_b, M_b in the BODY frame, ready to sum with gravity in the EOM.
%
%   TERRAIN: a plane with unit normal n = P.terrain.n_down (pointing INTO the
%   ground) through P.terrain.p0. Everything is measured relative to that
%   plane, not relative to inertial z:
%     penetration      delta    = n . (P_foot - p0)
%     penetration rate deltadot = n . V_foot
%     normal force     acts along -n (out of the ground)
%     friction         acts in the tangent plane, opposing the foot's
%                      tangential velocity  v_t = V_foot - (n.V_foot) n
%   For flat ground n = [0;0;1] and all of this reduces EXACTLY to the
%   original z-based formulation.
%
%   FRAME DISCIPLINE -- the direct lesson from the UH-60 LandingGear.m bug
%   (LANDING_GEAR_BUG_REPORT.md). On a slope these rules are what matter:
%     (1) the normal and friction forces are built together as ONE inertial
%         vector, then rotated into BODY frame BEFORE summing. Never
%         concatenate a terrain-frame component with a body-frame component.
%     (2) moments are taken about the CONTACT POINT. For this rigid-leg model
%         the contact point IS p_foot_b.

V_b     = x(1:3);
omega_b = x(4:6);
P_cg_i  = x(10:12);
n       = P.terrain.n_down;
p0      = P.terrain.p0;

F_b = zeros(3,1);
M_b = zeros(3,1);

N       = zeros(1, P.nLegs);
delta   = zeros(1, P.nLegs);
contact = false(1, P.nLegs);
F_t     = zeros(3, P.nLegs);      % friction, INERTIAL, lies in the tangent plane

for i = 1:P.nLegs
    pf_b = P.p_foot_b(:,i);

    % (1) foot position, inertial
    P_foot_i = P_cg_i + R_b2i*pf_b;

    % (2) penetration along the terrain normal (> 0 means below ground)
    delta(i) = n.' * (P_foot_i - p0);

    if delta(i) <= 0
        continue                        % airborne leg: zero force, zero moment
    end
    contact(i) = true;

    % (3) foot velocity, inertial (velocity of a body-fixed point)
    V_foot_i = R_b2i*( V_b + cross(omega_b, pf_b) );
    deltadot = n.' * V_foot_i;

    % (4) normal force magnitude. CLAMP AT ZERO: ground pushes, never pulls.
    N(i) = max( P.k*delta(i) + P.c*deltadot, 0 );

    % (5) friction -- SWAPPABLE MODEL (physics or ML).
    %     Input is the foot's velocity projected onto the tangent plane.
    v_t      = V_foot_i - deltadot*n;
    F_t(:,i) = frictionFcn(v_t, N(i), P);

    % (6) total foot force on the aircraft, INERTIAL, built as ONE vector:
    %     friction in the tangent plane + normal force out of the ground (-n)
    F_foot_i = F_t(:,i) - N(i)*n;

    % (7) RULE 1: convert to body frame BEFORE summing
    F_foot_b = R_b2i.' * F_foot_i;

    % (8) RULE 2: moment arm is the contact point
    M_foot_b = cross( pf_b, F_foot_b );

    F_b = F_b + F_foot_b;
    M_b = M_b + M_foot_b;
end

aux = struct('N', N, 'delta', delta, 'contact', contact, 'F_t', F_t);
end
