function Pmn = crossPowerFlux(dat)
% crossPowerFlux - Calculate the cross-power-flux between all of the modes as defined 
% by Auld. 
% 
% B. A. Auld, Acoustic Fields and Waves in Solids 2, 2nd ed., vol. 2, 2 vols. 
% Malabar: Krieger Publishing Company, 1990.

if ~isscalar(dat) % compute recursively for every waveguide problem in the vector "dat"
    Pmn = arrayfun(@crossPowerFlux,dat,'UniformOutput',false); % apply to every object in the array "dat"
    return;
end

if any(dat.w ~= dat.w(1,:), 'all') % orthogonality for constant frequency! 
    error('GEWTOOL:crossPowerFlux:nonconstfreq', 'Cross power flux must be computed for wavenumbers at constant frequency.');
end

s = size(dat.k); nK = s(1); nF = s(2);
v = velocity(dat);
T = stress(dat);
Imn = cell(1, dat.gew.geom.nLay);
for l = 1:length(dat.gew.lay)
    vn = permute(v{l}, [5 1 2 3 4]); % order as: [singleton x nK x nF x zi x v]
    TLay = permute(T{l}, [6 1 2 3 5 4]); % order as: [singleton x nK x nF x zi x T']
    txn = TLay(:,:,:,:,:,1); % traction ex.T = T'.ex of size [singleton x nK x nF x zi x tx]
    
    vm = permute(vn, [2 1 3 4 5]); % transpose: [nK x singleton x nF x zi x v]
    txm = permute(txn, [2 1 3 4 5]); % transpose: [nK x singleton x nF x zi x T']
    
    Imn{l} = sum(-conj(vn).*txm - vm.*conj(txn), 5); % cross power flux density
end

Pmn = 1/4*GEWintegrate(dat.gew, Imn, [], 4); % size: [nK, nK, nF]

end