function P = landing_params()
%LANDING_PARAMS  All parameters for the simplified landing model.
%
%   SI units throughout (m, kg, N, rad, s) -- deliberately NOT the ft/slug/lb
%   of the parent UH-60 sim, to avoid unit-mixing errors. This model is about
%   method, not matching UH-60 numbers.
%
%   FRAMES: NED inertial (x North, y East, z DOWN), ground plane at z = 0.
%           Body: x forward, y right, z down, origin at CG.
%           ALTITUDE = -zD.  A foot is penetrating when its inertial z > 0.

%% Gravity
P.g = 9.81;                     % m/s^2

%% Body -- rectangular box
P.m = 1000;                     % kg 
P.a = 10;                       % m, length (x, fore-aft)
P.b =  2;                       % m, width  (y, lateral)
P.h =  2;                       % m, height (z, vertical)

Ixx = P.m*(P.b^2 + P.h^2)/12;   %  666.67 kg*m^2
Iyy = P.m*(P.a^2 + P.h^2)/12;   % 8666.67
Izz = P.m*(P.a^2 + P.b^2)/12;   % 8666.67
P.I    = diag([Ixx Iyy Izz]);
P.Iinv = diag(1./[Ixx Iyy Izz]);

%% Landing gear -- foot positions in BODY frame [m]
%  Two forward legs splayed outboard, one aft (helicopter-like tricycle).
%  +z is BELOW the CG (NED body frame).
%
%  NOTE ON GEOMETRY: sum(x_i) = 2 + 2 - 3 = +1 m, so the leg triangle's
%  centroid is NOT under the CG. Static equilibrium therefore has UNEQUAL
%  leg loads and a small nose-up pitch -- see test_step2_contact.m. This is
%  physically correct (real rotorcraft do this) and is retained deliberately.
P.p_foot_b = [ +2.0, +2.0, -3.0; ...   % x: fore/aft
               -1.5, +1.5,  0.0; ...   % y: lateral (outboard splay)
               +1.8, +1.8, +1.8 ];     % z: below CG
P.nLegs    = size(P.p_foot_b, 2);
P.legNames = {'front-left','front-right','aft'};

%% Leg spring / damper -- identical on all legs to start
%  k sized for ~5 cm MEAN static deflection: 3*k*d = m*g -> k = m*g/(3*0.05)
P.k = 6.54e4;                   % N/m    per leg
P.c = 2800;                     % N*s/m  per leg  (zeta ~ 0.30 in heave)

%% Friction -- ACTIVE from step 3. Step-2 tests override mu = 0 explicitly.
P.mu   = 0.6;                   % [-]  Coulomb coefficient
P.v_bl = 0.1*0.3048;            % m/s, regularization boundary layer = 0.1 ft/s,
                                %      matching the UH-60 model (LG_VBL, InitLG.m).
                                %      Was 0.01 m/s, which is RK4-unstable at
                                %      dt = 0.001 (see PLAN.md Sec.6 correction).

%% Ground -- planar, possibly inclined. Default is flat.
%  make_terrain(slope_deg, uphill_heading_deg). Flat ground reproduces the
%  original z = 0 plane exactly, so steps 1-2 are unaffected.
P.terrain = make_terrain(0, 0);

%% Integration timesteps (PLAN.md Sec.6)
P.dt_physics = 0.001;           % s
P.dt_ml      = 0.01;            % s

%% Contact master switch -- step 1 disables this for pure free-fall
P.contact_on = true;
end
