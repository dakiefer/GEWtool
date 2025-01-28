function [u] = fieldComponents(dat)
% fieldComponents - restructure eigenvectors Psi into each of the field components. 
% Example: if Psi = [3*N x 1] is a vector containing ux, uy and uz components,
% then u = [N x 3] array where u(:,i) = ui. In practice we have eigenvectors at
% Nk wavenumbers with Nw frequencies. Then the above Psi is of size [Nk x Nw x 3*N], 
% while u is of size [Nk x Nw x N x 3]. The actual u returned by this function
% is a cell array with components for each layer. The ux displacements across
% the layer lind = 3 of the mode at kind = 5 and wind = 10 is
% >> u{lind}(kind,wind,:,1)
% 
% Usage: 
% > u = fieldComponents(dat);
% 
% Arguments: 
% - dat:    GEWdat object that stores the computed modes. 
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    u = arrayfun(@fieldComponents,dat,'UniformOutput',false); % apply to every object in the arrays "dat"
    return;
end

Psi = zeros(dat.Nk, dat.Nw, dat.gew.geom.Ndof); % allocate with zeros
Psi(:,:,dat.gew.geom.gdofFree) = dat.Psi;       % expand Psi with Dirichlet-BC nodes

u = cell(dat.gew.geom.nLay, 1); % allocate 
for l = 1:dat.gew.geom.nLay    
    ulay = Psi(:,:,dat.gew.geom.gdofOfLay{l});
    sz =  [dat.Nk, dat.Nw, dat.gew.geom.N(l), dat.gew.geom.Nudof(l)]; 
    u{l} = reshape(ulay, sz);
end

end
