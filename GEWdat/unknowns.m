function [u] = unknowns(gew,dat)
% unknowns - restructure eigenvectors Psi into each of the unknown components. 
% Example: if Psi = [3*N x 1] is a vector containing ux, uy and uz components,
% then u = [N x 3] array where u(:,i) = ui. In practice we have eigenvectors at
% Nk wavenumbers with Nw frequencies. Then the above Psi is of size [Nk x Nw x 3*N], 
% while u is of size [Nk x Nw x N x 3]. The actual u returned by this function
% is a cell array with components for each layer. The ux displacements across
% the layer lind = 3 of the mode at kind = 5 and wind = 10 is
% >> u{lind}(kind,wind,:,1)
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

u = cell(gew.geom.nLay, 1); % initialize
for l = 1:gew.geom.nLay
    ulay = dat.Psi(:,:,gew.geom.gdofOfLay{l});
    sz =  [size(dat.k), gew.geom.N(l), gew.geom.Nudof(l)]; 
    u{l} = reshape(ulay, sz);
end

end
