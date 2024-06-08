function [S] = strain(gew,dat)
% strain - Strain tensor S.
% 
% The strain S is the symmetric part of the displacement gradient F = grad(u), i.e.,
% 
% S = 1/2(F + F^T)
% 
% In a plate (Cartesian coordinates) the displacement gradient can be computed as
% 
% F = ex ik u +  ey ∂u∂y 
% 
% In cylindrical coordinates or 3d waveguides, a different formulation is
% needed.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if isa(gew,"Cylinder")
    warning('Cylinders do not support this function yet. The results will be wrong.');
end

S = cell(gew.geom.nLay, 1); % allocate for each layer
for l = 1:gew.geom.nLay
    sizeF = [size(dat.w), gew.geom.N(l), gew.geom.Nudof(l), gew.geom.Nudof(l)]; % size of F = grad u
    F = zeros(sizeF); % allocate for displacement gradient F = grad u
    lay = gew.lay(l);
    Dy = 1/lay.h*lay.D1; % differentiation matrix
    iku = 1i*dat.k.*dat.u{l}; % ex.F = ik*u
    uu = permute(dat.u{l}, [1, 2, 5, 3, 4]); % additional dimension for mult. with diff. mat.
    dyu = sum(shiftdim(Dy, -2).*uu, 4); % ey.F = ∂u/∂y, dimension 4 is singleton
    F(:,:,:,1,:) = iku;  % assign components ex.F
    F(:,:,:,2,:) = dyu;  % assign components ey.F
    S{l} = 1/2*(F + permute(F, [1,2,3,5,4])); % symmetric part of gradient
end 

end % end function
