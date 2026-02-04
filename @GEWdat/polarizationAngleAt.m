function [gamma] = polarizationAngleAt(dat, at)
%polarizationAngle - Angle between the displacement vector u and the shear-horizontal direction ey. 
% 
% For mode n with displacement vector un(z), the polarization angle gamma(z) at 
% coordinate z is defined by
%
% cos(gamma) = abs{ey.un}/||un||
% 
% where ey denotes the unit directional vector in transverse direction (to k).
% 
% Usage: 
% > gamma = polarizationAngleAt(dat,'top');
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%

if ~isscalar(dat) % compute recursively 
    compute = @(datObj) polarizationAngleAt(datObj, at); % function to apply
    gamma = arrayfun(compute,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

% reduced polarization was computed
if length(dat.gew.udof) < 3 
    if all(dat.gew.udof == dat.gew.udofLamb)
        gamma = pi/2*ones(size(dat.k));
    elseif all(dat.gew.udof == dat.gew.udofSH)
        gamma = zeros(size(dat.k));
    else
        error('GEWdat:polarizationAngleAt','Unknown polarization. Use this function on fullyCoupled waves.'); 
    end
    return;
end 

% fullyCoupled waves were computed
if isnumeric(at)
    [l, n, ~] = layerNodeCompOf(dat.gew.geom, at); 
elseif ischar(at)
    switch at
        case {'top','t','outer','o'}
            n = dat.gew.geom.N(end);  % node index
            l = dat.gew.geom.nLay;    % layer index
        case {'bottom','b','inner','i'}
            n = 1; % node index
            l = 1; % layer index 
        otherwise 
            error('GEWTOOL:polarizationAngleAt',...
                'Second argument should be one of "top" or "bottom"');
    end
end

uu = displacement(dat); % displacements of all modes
u = uu{l}(:,:,n,:);     % displacement at surface of mode n
eyDotU = abs( u(:,:,:,dat.gew.udofSH) );  % abs(ey.u) (ey is SH-component)
Unorm = vecnorm(u,2,4); % 2-norm along 4th-dimension
gamma = real(acos(eyDotU./Unorm)); % complex due to numerical inpresition

end
