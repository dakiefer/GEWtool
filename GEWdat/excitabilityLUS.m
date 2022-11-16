function exc = excitabilityLUS(gew, dat, at)
% excitabilityLUS - laser-ultrasound excitability-detectability
%   An approximation for how well the waves can be excited and detected with a
%   laser ultrasonic excitation and interferometric detection. The
%   excitability/detectability depends on the modal displacements on the waveguide's
%   surface where the measurement is performed:
%   excitability  ~ ux (shear tractions are induced by laser)
%   detectability ~ uy (normal displacements are detected) 
%   excitabilityLUS() returns the product ux*uy normlized such that an
%   excitability of 1 is obtained 1e2 above the median excitability. Ploting
%   decibel values in [-40, 0] usually gives nice pictures. 
%
%   Usage: 
%   exc = excitabilityLUS(gew, dat, at):   provide the waveguide description
%   "gew", the dispersion data "dat" and the surface where to compute exc: 
%       - at: one of 'top', 'bottom'.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris

    switch at
        case 'top'
            n = gew.geom.N(end);  % node index
            l = gew.geom.nLay;    % layer index
        case 'bottom'
            n = 1; % node index
            l = 1; % layer index 
        otherwise 
            error('GEWtool:excitabilityLUS',...
                'Third argument should be one of "top" or "bottom"');
    end

    dat = normalizeReal(gew, dat);
    ux = dat.u{l}(:,:,n,1);
    uy = dat.u{l}(:,:,n,2);
    exc = abs(ux).*abs(uy); % excitability ~ ux, detectability ~ uy
    excUnit = median(exc(:))*1e2;
    exc = exc./excUnit; % normalized excitability (don't use max as there might be singularities where cg->0)
end
