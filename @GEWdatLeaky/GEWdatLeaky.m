classdef GEWdatLeaky < GEWdat
% GEWdatLeaky - Store and postprocess quasi-guided waves solutions.
% 
% Quasi-guided waves are governed by the following EVP that is nonlinear in the
% eigenvalue ik: 
% 
% [ (ik)^2*L2 + ik*L1 + L0 + w^2*M + ∑j ibetaj*Rj ]*q = 0
% 
% where q:       eigenvector containing the displacements in the plate and the "nA"
%                bulk wave amplitudes "Aj".
%       ik:      eigenvalue -> horizontal wavenumber (times i)  
%       ibetaj:  vertical wavenumber (times i) of the jth bulk wave
%       w:       angular frequency (parameter) 
%       Li,M,Rj: n x n-matrices stored in "obj.gew.opNonlin"
% 
% The above problem is nonlinear in ik because the vertical wavenumbers satisfy
% dispersion relations of the form 
% 
% ik^2 + ibeta^2 = (iw/cj)^2 with bulk wave velocities "cj". 
% 
% The nonlinear eigenvalue problem is transformed to a polynomial one by
% introducing the new eigenvectors:
% 
% Psi = [ ibetaj q ]
%       [    q     ] 
% 
% The above is applied recursively for each bulk wave j, i.e., with 2 bulk waves
% the eigenvector has four blocks: 
% 
% Psi = [ ibeta1 ibeta2 q ]  <- block 1
%       [    ibeta2 q     ]  <- block 2
%       [    ibeta1 q     ]  <- block 3
%       [        q        ]  <- block 4
% 
% For radiation into fluid media, the above transformation leads to a quadratic
% eigenvalue problem while radiation into solids leads to a cubic eigenvalue
% problem. In the latter case, GEWtool applies a partial companion linearization
% to reduce the problem to a quadratic one. In both cases, GEWtool finaly solves
% an eigenvalue problem of the form
% 
% [ (ik)^2*L2 + ik*L1 + L0 + w^2*M ]*Psi = 0
% 
% with Psi: eigenvector stored in "obj.Psi"
%      ik: eigenvalue (as before) stored in "obj.k"
%      Li,M: matrices stored in "obj.gew.op"
%
% 2025 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

properties 
    nA    % (scalar integer) number of bulk waves in the exterior half spaces
    nBeta % (scalar integer) number of *non-eliminated* vertical wave numbers
end

properties (Dependent)
    A      % [Nk x Nw x nA] bulk wave amplitudes in the exterior half spaces
    beta   % [Nk x Nw x nA] vertical wavenumbers in the exterior half spaces
    q      % [Nk x Nw x n]  eigenvectors of the nonlinear eigenvalue problem
    pTop   % [Nk x Nw] outward power flux density through the top surface
    pBot   % [Nk x Nw] outward power flux density through the bottom surface
end

methods
    function obj = GEWdatLeaky(gew,k,w,Psi)
        % GEWdatLeaky - Create GEWdatLeaky object. 
        % Arguments: 
        % - gew:   object of class 'Waveguide'
        % - k:     [Nk x Nw] or [Nk x 1] wavenumbers in rad/m
        % - w:     [Nk x Nw] or [1 x Nw] angular frequencies in rad/s
        % - Psi:   [Nk x Nw x state-space-degreesOfFreedom] eigenvectors 
        if ~isa(gew,'PlateLeaky')
            error('GEWtool:GEWdatLeaky', 'The first argument needs to be an object of type PlateLeaky are a subclass thereof.'); 
        end
        obj = obj@GEWdat(gew,k,w,Psi);
        obj.nA = size(gew.opNonlin.L0,2) - numel(gew.geom.gdofFree); % number of bulk waves
        obj.nBeta = GEWdatLeaky.getNumberOfNonEliminatedBeta(gew.halfSpaces);
    end
    function q = get.q(obj)
        blockInd_q = 2^obj.nBeta;
        q = obj.getBlock(blockInd_q);
    end
    function A = get.A(obj)
        indA = [obj.gew.halfSpaces.dofA];
        A = obj.q(:,:,indA);
    end
    function beta = get.beta(obj)
        beta = nan([size(obj.k) obj.nBeta]); % allocate
        qHq = sum(conj(obj.q).*obj.q,3);
        for i = 1:obj.nBeta
            blockInd_ibq = getBlockIndexOfBetai(obj,i); % which one of the 2^nA blocks to chose
            ibq = obj.getBlock(blockInd_ibq);
            betaih = -1i*sum(conj(obj.q).*ibq,3)./qHq;
            beta(:,:,i) = betaih/obj.gew.np.h0; % vertical wavenumbers for the ith bulk wave
        end
    end
    function pTop = get.pTop(obj)
        pTop = obj.powerFluxThroughSurf("top"); 
    end
    function pBot = get.pBot(obj)
        pBot = obj.powerFluxThroughSurf("bottom"); 
    end

    function beta = getBetasAt(obj,surf)
        warning('temporary implementation')
        if surf == "bottom"
            ind = 1; 
        else 
            ind = 2; 
        end
        beta = obj.beta(:,:,ind);
    end

    function pOut = powerFluxThroughSurf(obj,surf)
        loading = obj.getLoadingAt(surf); 
        if isempty(loading)
            pOut = 0; 
            return; 
        end
        Asurf = getBulkWaveAmplitudeAt(obj,surf);
        uSurf = obj.q(:,:,loading.dofU); 
        vSurf = -1i*obj.w.*uSurf; 
        if isa(loading.mat,'MaterialFluid')
            rhof = loading.mat.rho; 
            tSurf = -rhof*obj.w.^2.*Asurf; % assume top surface for now
            pOut = -1/2*real( conj(vSurf).*tSurf ); 
        elseif isa(loading.mat,'MaterialIsotropic')
            coupl = obj.gew.couplingMatricesSolid(loading,obj.gew.np);
            udof = obj.gew.udof;
            k = obj.k; g = obj.beta(:,:,1); e = obj.beta(:,:,2); w = obj.w;
            % en = obj.gew.dofOutofplane(udof);
            Tk2 = shiftdim(coupl.Tk2(udof,udof),-2);
            Tkg = shiftdim(coupl.Tkg(udof,udof),-2);
            Tke = shiftdim(coupl.Tke(udof,udof),-2);
            Tw2 = shiftdim(coupl.Tw2(udof,udof),-2);
            dotA = @(T) sum( T.*permute(Asurf,[1 2 4 3]), 4 ); % contraction of T with Asurf: = T.Asurf
            tauz = k.^2.*dotA(Tk2) + k.*g.*dotA(Tkg) + k.*e.*dotA(Tke) - w.^2.*dotA(Tw2); 
            pOut = -1/2*real(sum( conj(vSurf).*tauz,3 ));
        else
            error('GEWdatLeaky:powerFluxThroughSurf','Unknown loading type at %s. It should be of type MaterialFluid or MaterialIsotropic.',surf); 
        end
        if surf == "bottom"
            pOut = -pOut; % correct sign
        end
    end

    function Asurf = getBulkWaveAmplitudeAt(obj,surf)
        loading = obj.getLoadingAt(surf); 
        Asurf = obj.q(:,:,loading.dofA); %.*exp(-1i*beta*zItf);
    end

    function loading = getLoadingAt(obj,surf)
        isAtSurf = {obj.gew.halfSpaces.at}; % collect into a cell array
        ind = cellfun(@(x)x==surf, isAtSurf); % compare each of the entries
        loading = obj.gew.halfSpaces(ind); 
    end

    function block = getBlock(obj,j)
        n = size(obj.gew.opNonlin.L0,2); % size of each block 
        if j > 2^obj.nBeta
            error('GEWdatLeaky:getBlock', 'Index out of range. There exist only %d blocks.', obj.nBeta + 1); 
        end
        ind_block_j = (1:n) + (j-1)*n; 
        block = obj.Psi(:,:,ind_block_j);
    end

    function indBetaBlock = getBlockIndexOfBetai(obj,i)
        j = obj.nBeta:-1:i; 
        indBetaBlock = sum((2.^j)/2); 
    end
end

methods (Static)
    function nBeta = getNumberOfNonEliminatedBeta(loading)
        nBeta = 0; 
        for i = 1:length(loading)
            nBeta_i = length(loading(i).dofA); % number of beta == number of bulk wave amplitudes
            if ~loading(i).eliminated
                nBeta = nBeta + nBeta_i;
            end
        end
    end
end

end