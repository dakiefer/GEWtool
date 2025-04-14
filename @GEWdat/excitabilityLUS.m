function exc = excitabilityLUS(dat, at)
% excitabilityLUS - laser-ultrasound excitability-detectability
%   An approximation for how well the waves can be excited and detected with a
%   laser ultrasonic excitation and interferometric detection. The
%   excitability/detectability depends on the modal displacements on the waveguide's
%   surface where the measurement is performed:
%   - excitability  ~ vx (shear tractions are induced by small-diameter laser in 
%                   nonablasive regime)
%   - detectability ~ uz (normal displacements are detected) 
%   excitabilityLUS() returns the product vx*uz normlized such that an
%   excitability of 1 is obtained 1e2 above the median excitability. To mantain
%   the same normalization with different discretiations (i.e., obtained number
%   of modes) it is important to either (i) specify the desired number of modes
%   in computeW(), or (ii) restrict to the solutions in a given frequency range.
%   Otherwise, more and more modes are considered with increasing discretization
%   order. Plotting decibel values in [-40, 0] usually gives nice pictures.
%
%   Usage: 
%   > exc = excitabilityLUS(dat, at); 
%
%   Arguments:
%   Provide the dispersion data "dat" and select the surface by chosing "at" as: 
%   - for the last layer's top surface:     'top' | 'outer' | 't' | 'o', or 
%   - for the first layer's bottom surface: 'bottom' | 'inner' | 'b' | 'i'.
% 
% 2022-2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    compute = @(datObj) excitabilityLUS(datObj, at); % function to apply
    exc = arrayfun(compute,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

n = size(dat.gew.op.M,1);
numModes = size(dat.k);
if any(numModes == n) || any(numModes == 2*n)
    warning('GEWTOOL:excitabilityLUS:restrictSolution',...
        ['It seems that you have not restricted the number of modes. Note' ...
        ' that the normalization of the excitation depends on the number of' ...
        ' modes and, hence, changes with matrix size. When computing solutions,'...
        ' provide the number of desired modes "nModes", e.g.,\n'...
        '>> dat = computeW(gew, k, nModes);']);
end

switch at
    case {'top','t','outer','o'}
        n = dat.gew.geom.N(end);  % node index
        l = dat.gew.geom.nLay;    % layer index
    case {'bottom','b','inner','i'}
        n = 1; % node index
        l = 1; % layer index 
    otherwise 
        error('GEWTOOL:excitabilityLUS',...
            'Third argument should be one of "top" or "bottom"');
end

dat = normalizeReal(dat);
ux = dat.u{l}(:,:,n,dat.gew.udof == dat.gew.udofAxial);
uz = dat.u{l}(:,:,n,dat.gew.udof == 3); % udof = 3 is always normal/radial to the plate/cylinder
exc = abs(dat.k).*abs(ux).*abs(uz); % excitability ~vx, detectability ~ uy
excUnit = median(exc(isfinite(exc)))*1e2;
exc = exc./excUnit; % normalized excitability (don't use max as there might be singularities where cg->0)
end
