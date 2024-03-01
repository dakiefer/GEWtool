function plotRayCurve(obj, erot, varargin)
% plotRayCurve - plot the energy velocity curve when rotating around axis erot. 
% Iteratively calls ENERGYVEL for different propagation directions and
% plots the result into a polar plot.
%
% Usage:
% plotRayCurve(mat) % rotate around z-axis
% plotRayCurve(mat, erot) % rotate around erot [3x1]
% plotRayCurve(mat, [], varargin) % pass plotting arguments (e.g. 'k.')
% 
% See also: wavespeeds, Material.energyVel, Material.
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, CNRS, France

if nargin >= 2, erot = erot(:)/norm(erot); end % normalize
% if nargin < 3, varargin{1} = '.'; end % default plotting style
if nargin < 2 | isempty(erot)
    erot = [0;0;1]; % default is rotation in x-y-plane
end

alpha = linspace(0, 2*pi, 1000); % default angles
[ce, ~] = obj.rayCurve(alpha, erot);
cex = squeeze(ce(1,:,:)).';
cey = squeeze(ce(2,:,:)).';
cez = squeeze(ce(3,:,:)).';
plot3(cex, cey, cez, varargin{:})
% legend({'sl', 'st1', 'st2'})
% title(sprintf('%s: slowness curves around rotation axis [%g, %g, %g]\nangle 0: [%g, %g, %g]',...
%     obj.name, erot(1), erot(2), erot(3), e0(1), e0(2), e0(3)))
end
