function [S] = strain(gew,dat)
% strain - Strain tensor S.
% 
% The strain S is the symmetric part of the displacement gradient F = grad(u), i.e.,
% 
% S = 1/2(F + F^T)
% 
% In a plate (Cartesian coordinates) the displacement gradient can be computed as
% 
% F = ex ik u +  ey ∂u∂y .
% 
% In a cylinder (Cylindrical coordinates) it is given by
% 
% F = ex ik u +  er ∂u∂r ephi 1/r (i*n*I + A)*u , 
% where I is the 3x3 identity matrix and A = ephi er - er ephi . This means that
% the displacement gradient has an angular component even if the displacement u
% does not.
% 
% In 3d waveguides, a different formulation is again needed.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) strain(gewObj, datObj); % function to apply
    S = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

S = cell(gew.geom.nLay, 1); % allocate for each layer
F = displGrad(gew,dat);     % each type of Waveguide "gew" knows how to compute their displacement gradient
for l = 1:gew.geom.nLay
    S{l} = 1/2*(F{l} + permute(F{l}, [1,2,3,5,4])); % symmetric part of gradient
end 

end % end function
