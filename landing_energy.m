function E = landing_energy(x, P)
%LANDING_ENERGY  Total mechanical energy of the system at one state.
%
%   E = landing_energy(x, P)  ->  struct with KE_trans, KE_rot, PE_grav,
%                                 PE_spring, total   [J]
%
%   CONSERVATION IS ONLY EXPECTED WHEN P.c == 0 AND P.mu == 0. With damping or
%   friction active this still returns a valid instantaneous energy, but it
%   will (correctly) decrease.
%
%   Gravity datum is the ground plane (z = 0), so PE_grav = m*g*altitude
%   = m*g*(-zD).
%
%   The spring is conservative: power delivered to the body is
%     F_foot . V_foot = -k*delta*deltadot = -d/dt( 0.5*k*delta^2 )
%   and the force/moment pair applied at the CG does exactly that same work,
%   since  F.V_cg + (r x F).omega = F.(V_cg + omega x r) = F.V_foot.

V_b     = x(1:3);
omega_b = x(4:6);
P_cg_i  = x(10:12);

R_b2i = rot_b2i(x(7), x(8), x(9));

E.KE_trans = 0.5 * P.m * (V_b' * V_b);
E.KE_rot   = 0.5 * (omega_b' * P.I * omega_b);
E.PE_grav  = P.m * P.g * (-P_cg_i(3));      % altitude = -zD

E.PE_spring = 0;
if P.contact_on
    for i = 1:P.nLegs
        P_foot_i = P_cg_i + R_b2i*P.p_foot_b(:,i);
        d = P.terrain.n_down.' * (P_foot_i - P.terrain.p0);   % penetration along terrain normal
        if d > 0
            E.PE_spring = E.PE_spring + 0.5*P.k*d^2;
        end
    end
end

E.total = E.KE_trans + E.KE_rot + E.PE_grav + E.PE_spring;
end
