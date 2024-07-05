function exc = excitabilityLUS(gew, dat, at)
% excitabilityLUS - laser-ultrasound excitability-detectability
%   An approximation for how well the waves can be excited and detected with a
%   laser ultrasonic excitation and interferometric detection. The
%   excitability/detectability depends on the modal displacements on the waveguide's
%   surface where the measurement is performed:
%   - excitability  ~ ux (shear tractions are induced by small-diamiter laser in 
%                   nonablasive regime)
%   - detectability ~ uy (normal displacements are detected) 
%   excitabilityLUS() returns the product ux*uy normlized such that an
%   excitability of 1 is obtained 1e2 above the median excitability. To mantain
%   the same normalization with different discretiations (i.e., obtained number
%   of modes) it is important to either (1) specify the desired number of modes
%   in computeW(), or (2) restrict to the solutions in a given frequency range.
%   Otherwise, more and more modes are considered with increasing discretization
%   order. Ploting decibel values in [-40, 0] usually gives nice pictures.
%
%   Usage: 
%   exc = excitabilityLUS(gew, dat, at):   provide the waveguide description
%   "gew", the dispersion data "dat" and select the surface by chosing "at": 
%       - for the last layer's surface:  'top' | 'outer' | 't' | 'o', or 
%       - for the first layer's surface: 'bottom' | 'inner' | 'b' | 'i'.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) excitabilityLUS(gewObj, datObj, at); % function to apply
    exc = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

n = size(gew.op.M,1);
numModes = size(dat.k);
if any(numModes == n) || any(numModes == 2*n)
    warning('GEWTOOL:excitabilityLUS:restrictSolution',...
        ['It seems that you have not restricted the number of modes. Note' ...
        ' that the normalization of the excitation depends on the number of' ...
        ' modes and, hence, changes with matrix size.']);
end

switch at
    case {'top','t','outer','o'}
        n = gew.geom.N(end);  % node index
        l = gew.geom.nLay;    % layer index
    case {'bottom','b','inner','i'}
        n = 1; % node index
        l = 1; % layer index 
    otherwise 
        error('GEWTOOL:excitabilityLUS',...
            'Third argument should be one of "top" or "bottom"');
end

dat = normalizeReal(gew, dat);
v = velocity(dat);
vx = v{l}(:,:,n,1);
uy = dat.u{l}(:,:,n,2);
exc = abs(vx).*abs(uy); % excitability ~vx, detectability ~ uy
excUnit = median(exc(isfinite(exc)))*1e2;
exc = exc./excUnit; % normalized excitability (don't use max as there might be singularities where cg->0)
end
