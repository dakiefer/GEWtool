function [ces, ea] = rayCurve(obj, alpha, erot)
% RAYCURVE - Compute energy velocities for a set of wave vector orientations.
%
% The wave vectors k (or slowness vectors) are in a plane orthogonal to the
% provided vector "erot". The wave vectors are rotatet around "erot" to the
% angles "alpha". For each of the orientations, the energy velocity vectors for
% each of the three bulk waves is computed.
% 
% Note: Applicable only to homogeneous plane wave in nondissipative media.
% 
% Literature: 
% [1] K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive
% Testing of Materials: Theoretical Foundations (translated from German),
% 1st ed. Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
% [2] D. Royer and T. Valier-Brasier, Elastic Waves in Solids 1: Propagation.
% ISTE and John Wiley & Sons, 2022.
%
% Arguments:
% - obj:    Material object 
% - alpha:  list of angles at which propagation direction of phases [3x1] (default: [1; 0; 0])
% - erot:   vector around which to rotate the wave vector k
%
% Return values:
% - ce:    [3x3xN]-array where ce(:,i,j) is the energy velocity vector for mode i.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, CNRS, France

if nargin == 3
    erot = erot(:)/norm(erot); % normalize
elseif nargin < 3
    erot = [0;0;1]; % default is rotation in x-y-plane
elseif nargin < 2
    error('GEWTOOL:rayCurve:wrongNumberOfArgs', 'Specify the rotation angles as second argument (vector).');
end
if ~isvector(alpha)
    error('GEWTOOL:rayCurve:wrongArg', 'The set of angles alpha should be specified as a vector.');
end

ebasis = null(erot.'); % svd to compute orthogonal vectors to erot
ea = ebasis(:,1); % first basis vector spanning the plane normal to erot
eb = ebasis(:,2); % second basis vector spanning the plane normal to erot

N = length(alpha);
ces = zeros(3,3,N); % initialize for 3 energy velocity vectors at every angle
for i = 1:N
    % Rodrigues' rotation formula (last term is zero due to e0.erot = 0): 
    eki = cos(alpha(i))*ea + sin(alpha(i))*eb; %  + (1 - cos(p0))*sum(e0.*erot)*erot;
    ces(:,:,i) = obj.energyVel(eki);
end

end