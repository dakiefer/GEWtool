function [T] = stress(wguide, dat)
%STRESS compute the stress T given the displacements u and wavenumber k.
% - u: displacements
% - k: wavenumbers
% - wguide: waveguide object
% 

T = cell(wguide.geom.nLay, 1); % allocate for each layer
for i = 1:wguide.geom.nLay
    lay = wguide.lay(i);
    D1 = 1/lay.h*lay.D1;
    c = lay.mat.c; % stiffness tensor
    udof = 1:wguide.geom.Nudof(i); % Lamb and/or SH
    
    uu = permute(dat.u{i}, [1, 2, 5, 3, 4]); % additional dimension for mult. with diff. mat.
    ud = sum(shiftdim(D1, -2).*uu, 4);       % du/dy
    uu = permute(uu, [1,2,4,3,6,5]);         % additional dimension for contraction with c
    ud = permute(ud, [1,2,3,4,6,5]);         % additional dimension for contraction with c
    cx = shiftdim(c(udof,udof,udof,1), -3);
    cy = shiftdim(c(udof,udof,udof,2), -3);
    T{i} = sum(1i*dat.k.*cx.*uu + cy.*ud, 6);   % contraction 
end

end
