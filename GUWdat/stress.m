function [T] = stress(wguide, dat)
%STRESS compute the stress T given the displacements u and wavenumber k.
% - u: displacements
% - k: wavenumbers
% - wguide: waveguide object
% 

warning('re-use strain instead?')
D1 = 1/wguide.geom.h*wguide.lay.D1; % TODO: multilayer
c = wguide.mat.c; % stiffness tensor
udof = 1:wguide.geom.Nudof; % Lamb and/or SH

uu = permute(dat.u, [1, 2, 5, 3, 4]); % assume this to be done already? 
ud = sum(shiftdim(D1, -2).*uu, 4);
uu = permute(uu, [1,2,4,3,6,5]);
ud = permute(ud, [1,2,3,4,6,5]);
cx = shiftdim(c(udof,udof,udof,1), -3);
cy = shiftdim(c(udof,udof,udof,2), -3);
T = sum(1i*dat.k.*cx.*uu + cy.*ud, 6);

end
