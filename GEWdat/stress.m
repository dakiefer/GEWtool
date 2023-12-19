function [T] = stress(gew, dat)
%STRESS compute the stress T given the displacements u and wavenumber k.
% - u: displacements
% - k: wavenumbers
% - wguide: waveguide object
% 

if isa(gew,"Cylinder")
    warning('Cylinders do not support this function yet. The results will be wrong.');
end
T = cell(gew.geom.nLay, 1); % allocate for each layer
for i = 1:gew.geom.nLay
    lay = gew.lay(i);
    D1 = 1/lay.h*lay.D1;
    c = lay.mat.c; % stiffness tensor
    udof = 1:gew.geom.Nudof(i); % Lamb and/or SH
    
    uu = permute(dat.u{i}, [1, 2, 5, 3, 4]); % additional dimension for mult. with diff. mat.
    ud = sum(shiftdim(D1, -2).*uu, 4);       % du/dy
    uu = permute(uu, [1,2,4,3,6,5]);         % additional dimension for contraction with c
    ud = permute(ud, [1,2,3,4,6,5]);         % additional dimension for contraction with c
    cx = shiftdim(c(udof,udof,udof,1), -3);
    cy = shiftdim(c(udof,udof,udof,2), -3);
    T{i} = sum(1i*dat.k.*cx.*uu + cy.*ud, 6);   % contraction 
end

% NOTE: the implementation that avoids calling strain() is a bit faster.

end
