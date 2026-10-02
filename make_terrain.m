function T = make_terrain(slope_deg, uphill_heading_deg)
%MAKE_TERRAIN  Define a planar (possibly inclined) ground surface.
%
%   T = make_terrain(slope_deg, uphill_heading_deg)
%
%     slope_deg           incline angle alpha of the plane [deg], 0 = flat
%     uphill_heading_deg  compass direction the ground RISES toward [deg],
%                         0 = North (+x), 90 = East (+y). Default 0.
%
%   The plane passes through the inertial origin. Its altitude rises along
%   the uphill heading:  altitude(x,y) = tan(alpha) * (uhat . [x;y]).
%
%   Returns struct T:
%     T.n_down   3x1 unit normal pointing INTO the ground (NED inertial).
%                Flat ground -> [0;0;1], i.e. straight down.
%     T.p0       3x1 a point on the plane (the origin).
%     T.t_up     3x1 unit vector along the plane pointing UPHILL.
%     T.slope_deg, T.uphill_heading_deg   the inputs, for labelling.
%
%   Everything contact-related uses ONLY n_down and p0:
%     penetration  delta = n_down . (P_foot - p0)     (> 0 means below ground)
%   For flat ground this reduces exactly to delta = z_foot, so the flat-ground
%   model is recovered bit-for-bit.

if nargin < 2 || isempty(uphill_heading_deg), uphill_heading_deg = 0; end

a   = slope_deg          * pi/180;
chi = uphill_heading_deg * pi/180;

% n_down is proportional to [tan(a)cos(chi); tan(a)sin(chi); 1]; normalized:
T.n_down = [ sin(a)*cos(chi); sin(a)*sin(chi); cos(a) ];

% Uphill direction in the plane: horizontal component uhat, rising by tan(a),
% i.e. zD DECREASING (NED). Normalized:
T.t_up   = [ cos(a)*cos(chi); cos(a)*sin(chi); -sin(a) ];

T.p0 = [0; 0; 0];
T.slope_deg          = slope_deg;
T.uphill_heading_deg = uphill_heading_deg;

% Guard: t_up must lie in the plane and both vectors must be unit length
assert(abs(dot(T.n_down, T.t_up)) < 1e-14, 'make_terrain: t_up not in plane');
assert(abs(norm(T.n_down) - 1)    < 1e-14, 'make_terrain: n_down not unit');
assert(abs(norm(T.t_up)   - 1)    < 1e-14, 'make_terrain: t_up not unit');
end
