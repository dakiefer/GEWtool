function psi = eigenVecs(gew, dat)
% EIGENVEC - Convert displacement array to eigenvector. 
% 
% This function is used for postprocessing, e.g., in groupVel(). The eigenvector
% psi is a re-arrangement of the displacement array dat.u. Furthermore, the
% degrees of freedom fixed by Dirichlet boundary conditions are removed.
% 
% Arguments: 
% - gew:     (Waveguide) Model of the waves
% - dat:     (struct) Data structure with field "u", as returned by computeK() or 
%            computeW()
% 
% Return value: 
% - psi:     (nK x nW x Ndof array) psi(i,j,:) is the eigenvector at the i-th wavenumber
%            and the j-th frequency.
% 
% % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
    compute = @(gewObj,datObj) eigenVecs(gewObj, datObj); % function to apply
    psi = arrayfun(compute,gew,dat,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
    return;
end

sizeU = size(dat.u{1});
psi = zeros([sizeU(1:2), gew.geom.Ndof]); % allocate
gdofsAccum = [];
for l = 1:gew.geom.nLay % convert structured u into unstructured u
    gdofLay = gew.geom.gdofOfLay{l}; % where to put into the global u vector
    [gdofNew, ldofNew] = setdiff(gdofLay, gdofsAccum); % remove coincident nodes
    gdofsAccum = [gdofsAccum, gdofNew]; % remember already treated dofs
    ulay = reshape(dat.u{l}, sizeU(1), sizeU(2), []);
    psi(:,:,gdofNew) = ulay(:,:,ldofNew);
end
psi(:,:,gew.geom.gdofDBC) = []; % remove homogeneous DBC nodes

end
