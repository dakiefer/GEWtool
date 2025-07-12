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
    nA % (scalar integer) number of bulk waves in the exterior half spaces
end

properties (Dependent)
    A     % [Nk x Nw x nA] bulk wave amplitudes in the exterior half spaces
    beta  % [Nk x Nw x nA] vertical wavenumbers in the exterior half spaces
    q     % [Nk x Nw x n]  eigenvectors of the nonlinear eigenvalue problem
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
    end
    function q = get.q(obj)
        blockInd_q = 2^obj.nA;
        q = obj.getBlock(blockInd_q);
    end
    function A = get.A(obj)
        n = size(obj.gew.opNonlin.L0,2);
        indA = (n-obj.nA+1):n;
        A = obj.q(:,:,indA);
    end
    function beta = get.beta(obj)
        beta = nan([size(obj.k) obj.nA]); % allocate
        qHq = sum(conj(obj.q).*obj.q,3);
        for i = 1:obj.nA
            blockInd_ibq = getBlockIndexOfBetai(obj,i); % which one of the 2^nA blocks to chose
            ibq = obj.getBlock(blockInd_ibq);
            betaih = -1i*sum(conj(obj.q).*ibq,3)./qHq;
            beta(:,:,i) = betaih/obj.gew.np.h0; % horizontal wavenumbers for the ith bulk wave
        end
    end
    function block = getBlock(obj,j)
        n = size(obj.gew.opNonlin.L0,2); % size of each block 
        if j > 2^obj.nA
            error('GEWdatLeaky:getBlock', 'Index out of range. There exist only %d blocks.', obj.nA + 1); 
        end
        if n*(2^obj.nA) ~= size(obj.Psi,3)
            error('GEWdatLeaky:getBlock', 'Size inconsistency. I expected %d blocks.', n*(obj.nA+1));
        end
        ind_block_j = (1:n) + (j-1)*n; 
        block = obj.Psi(:,:,ind_block_j);
    end
    function indBetaBlock = getBlockIndexOfBetai(obj,i)
        j = obj.nA:-1:i; 
        indBetaBlock = sum((2.^j)/2); 
    end
end

end