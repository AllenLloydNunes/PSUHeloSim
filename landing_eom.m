function [xdot, aux] = landing_eom(t, x, P, frictionFcn)
%LANDING_EOM  12-state rigid-body equations of motion, xdot = f(F, M).
%
%   [xdot, aux] = landing_eom(t, x, P, frictionFcn)
%
%   STATE (12):
%     x( 1: 3) = [u; v; w]          body-frame velocity            [m/s]
%     x( 4: 6) = [p; q; r]          body-frame angular rate        [rad/s]
%     x( 7: 9) = [phi; theta; psi]  Euler angles, 3-2-1            [rad]
%     x(10:12) = [xN; yE; zD]       inertial NED position          [m]
%
%   ALTITUDE = -x(12).
%
%   't' is accepted for integrator compatibility; the dynamics are autonomous.

%#ok<INUSL>

if nargin < 4 || isempty(frictionFcn)
    frictionFcn = @friction_physics;
end

V_b     = x(1:3);
omega_b = x(4:6);
phi     = x(7);
theta   = x(8);
psi     = x(9);

% Single source of truth for the transform (see rot_b2i.m)
R_b2i = rot_b2i(phi, theta, psi);

%% Gravity, expressed in body frame.
%  In NED, gravity points along +z inertial. Deriving it from R_b2i (rather
%  than writing the closed form separately) guarantees it can never drift out
%  of sync with the rotation matrix.
F_grav_b = R_b2i' * [0; 0; P.m*P.g];

%% Contact
if P.contact_on
    [F_contact_b, M_contact_b, aux] = contact_forces(x, R_b2i, P, frictionFcn);
else
    F_contact_b = zeros(3,1);
    M_contact_b = zeros(3,1);
    aux = struct('N',       zeros(1,P.nLegs), ...
                 'delta',   zeros(1,P.nLegs), ...
                 'contact', false(1,P.nLegs), ...
                 'F_t',     zeros(3,P.nLegs));
end

F_b = F_grav_b + F_contact_b;
M_b = M_contact_b;              % gravity acts at the CG -> contributes no moment

%% Translational: m*(Vdot_b + omega x V_b) = F_b
Vdot_b = F_b/P.m - cross(omega_b, V_b);

%% Rotational (Euler): I*omegadot + omega x (I*omega) = M_b
omegadot_b = P.Iinv * ( M_b - cross(omega_b, P.I*omega_b) );

%% Euler kinematics (3-2-1): body rates -> Euler angle rates
sP = sin(phi);   cP = cos(phi);
cT = cos(theta); tT = tan(theta);
p = omega_b(1);  q = omega_b(2);  r = omega_b(3);

phidot   = p + sP*tT*q + cP*tT*r;
thetadot =        cP*q -    sP*r;
psidot   =   (sP/cT)*q + (cP/cT)*r;

%% Position kinematics
Pdot_i = R_b2i * V_b;

xdot = [ Vdot_b; omegadot_b; phidot; thetadot; psidot; Pdot_i ];
end
