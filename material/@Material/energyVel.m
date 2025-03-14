function [ce] = energyVel(obj, ek)
% ENERGYVEL - Compute the energy velocities.
% Each column of the return value "ce" is the energy velocity vector
% associated to the three bulk waves as obtained by wavespeeds().
% 
% Note: Applicable only to homogeneous plane wave in nondissipative media.
% 
% Literature: 
% K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive
% Testing of Materials: Theoretical Foundations (translated from German),
% 1st ed. Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
%
% Arguments:
% - obj:   Material object 
% - ek:    propagation direction of phases [3x1] (default: [1; 0; 0])
%
% Return values:
% - ce:    [3x3]-array where each column is an energy velocity vector
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, CNRS, France

if isDissipative(obj)
    warning('This function is only implemented for homogeneous waves in nondissipative media.')
end
if nargin < 2, ek = [1; 0; 0]; end
validateattributes(ek, {'numeric'}, {'vector', 'numel', 3});

ek = ek(:)/norm(ek(:)); % normalized column vector
[cp, eu] = wavespeeds(obj,ek);
svec = 1./cp.'.*ek; % slowness vectors for unit frequency w = 1

c = obj.c; % 4th-order stiffness tensor
rho = obj.rho; % mass density
ce = nan(3,3); % initialize
for i = 1:3
    sveci4 = shiftdim(svec(:,i),-3); % slowness vector at forth tensor dimension
    eui3 = shiftdim(eu(:,i),-2); % polarization vector at third tensor timension
    ce(:,i) = sum(sum(sum(1/rho*eu(:,i).*c.*sveci4.*eui3,1),4),3).';
end
end
