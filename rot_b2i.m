function R = rot_b2i(phi, theta, psi)
%ROT_B2I  Rotation matrix from BODY frame to INERTIAL (NED) frame.
%
%   Euler sequence 3-2-1 (yaw psi -> pitch theta -> roll phi).
%
%   FRAMES (see PLAN.md Sec.1):
%     Inertial: NED -- x North, y East, z DOWN. Ground plane at z = 0.
%               ALTITUDE = -z.
%     Body:     x forward, y right, z down, origin at CG.
%
%   This is the SINGLE definition of the rotation matrix for the whole model.
%   Every other file must call this rather than rebuilding it inline -- a
%   duplicated/inconsistent transform is exactly the bug class that produced
%   the UH-60 LandingGear.m frame error.
%
%   Inputs:  phi, theta, psi   [rad]
%   Output:  R (3x3) such that  v_inertial = R * v_body

sP = sin(phi);   cP = cos(phi);
sT = sin(theta); cT = cos(theta);
sS = sin(psi);   cS = cos(psi);

R = [ cT*cS,  sP*sT*cS - cP*sS,  cP*sT*cS + sP*sS; ...
      cT*sS,  sP*sT*sS + cP*cS,  cP*sT*sS - sP*cS; ...
      -sT,    sP*cT,             cP*cT           ];
end
