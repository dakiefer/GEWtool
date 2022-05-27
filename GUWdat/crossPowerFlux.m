function Pmn = crossPowerFlux(guw, dat)
% crossPowerFlux Calculate the cross-power-flux of the waves as defined by Auld. 
% 
% B. A. Auld, Acoustic Fields and Waves in Solids 2, 2nd ed., vol. 2, 2 vols. 
% Malabar: Krieger Publishing Company, 1990.

if any(dat.w ~= dat.w(:,1), 'all') % orthogonality for constant frequency! 
    error('GUWtool: cross power flux must be computed for wavenumbers at constant frequency.');
end

s = size(dat.k); nF = s(1); nK = s(2);
Pmn = zeros(nK, nK, nF); % allocate: for each frequency nF compute Pnm
v = velocity(dat);
T = stress(guw, dat);
for l = 1:length(guw.lay)
    vn = permute(v{l}, [5 2 1 3 4]); % order as: [singleton x nK x nF x yi x v]
    TLay = permute(T{l}, [6 2 1 3 5 4]); % order as: [singleton x nK x nF x yi x T']
    txn = TLay(:,:,:,:,:,1); % traction ex.T = T'.ex of size [singleton x nK x nF x yi x tx]
    
    vm = permute(vn, [2 1 3 4 5]); % transpose: [nK x singleton x nF x yi x v]
    txm = permute(txn, [2 1 3 4 5]); % transpose: [nK x singleton x nF x yi x T']
    
    I = sum(-conj(vn).*txm - vm.*conj(txn), 5); % power flux density
    PmnLay = 1/4*chebintegrate(I, guw.geom.yItf(l,:), 4); % total power flux in layer l
    Pmn = Pmn + PmnLay; % cumulate 
end


end