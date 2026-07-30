function params = eaple_cartesian_to_polar(position)
    %EAPLE_CARTESIAN_TO_POLAR Convert an APLE/GIRD 3x1 Cartesian estimate to [r;omega;phi].
    % Position is metres. omega is wrapped to [0,2*pi), and phi is constrained to
    % [0,pi/2), which requires a nonnegative z coordinate for this paper's domain.
    arguments
    position (3,1) double {mustBeFinite}
end
r=norm(position); if r<=eps, error('eaple_cartesian_to_polar:ZeroRange','Position must have positive range.'); end
if position(3)<0, error('eaple_cartesian_to_polar:NegativeZ','NF_Loc polar domain requires z>=0.'); end
params=[r;mod(atan2(position(2),position(1)),2*pi);min(pi/2,acos(max(-1,min(1,position(3)/r))) )];
end
