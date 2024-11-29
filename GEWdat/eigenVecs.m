function psi = eigenVecs(gew, u)
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

% if ~isscalar(gew) % compute recursively for every waveguide problem in the vector "gew"
%     compute = @(gewObj,ui) eigenVecs(gewObj, ui); % function to apply
%     psi = arrayfun(compute,gew,u,'UniformOutput',false); % apply to every object in the arrays "gew" and "dat"
%     return;
% end

if isscalar(gew)
    u = {u}; % convert into cell of size 1x1 corresponding to the one "gew"
end
if length(gew) ~= length(u)
    error('GEWtool:eigenVecs:wrongSize', 'The length of ''gew'' and ''u'' needs to be equal.'); 
end

psi = cell(size(gew)); % allocate empty cell
for i = 1:length(gew)
    ui = u{i}; 
    sizeU = size(ui{1});
    psii = zeros([sizeU(1:2), gew(i).geom.Ndof]); % allocate
    gdofsAccum = [];
    for l = 1:gew(i).geom.nLay % convert structured u into unstructured u
        gdofLay = gew(i).geom.gdofOfLay{l}; % where to put into the global u vector
        [gdofNew, ldofNew] = setdiff(gdofLay, gdofsAccum); % remove coincident nodes
        gdofsAccum = [gdofsAccum, gdofNew]; % remember already treated dofs
        ulay = reshape(ui{l}, sizeU(1), sizeU(2), []);
        psii(:,:,gdofNew) = ulay(:,:,ldofNew);
    end
    psii(:,:,gew(i).geom.gdofDBC) = []; % remove homogeneous DBC nodes
    psi{i} = psii; 
end

if isscalar(gew) % unpack if scalar
    psi = psi{1}; 
end

end
