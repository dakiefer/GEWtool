function [kzgv, wzgv, k0s] = ZGV_Sylv_MFRDScan(L2, L1, L0, M, opts)
% ZGV_Sylv_MFRDScan - Compute ZGV points via iterative shift and search using fast eigenvalue solver.
%
% Returns ZGV points (k, w) with the smallest k for 
% [ (i*k)^2*L2 + i*k*L1 + L0 + mu*M ]*u = 0, where mu = w^2.
% The Method of Fixed Relative Distance (MFRD) proposed by Jarlebring, Kvall and
% Michiels [1] is used to this end. The method computes candidate wavenumbers close 
% to a shift k0 and refines them with a Newton-type iteration. By iteratively adapting
% the shift k0, a wavenumber range [kmin, kend] is searched. Details can be found in [2].
%
% It is important that matrix M is nonsingular and that the system is
% transformed in such way that k and mu are real.
% Input:
% - L2,L1,L0,M: matrices of the model 
% - options in opts:
%      - Neigs: (6) Number of eigenvalues to compute using eigs for a given shift.
%        Reducing Neigs can considerably speed up the computation, especially if
%        you have large matrices. If it is too small, you will miss solutions.
%      - kEnd: (40) Maximum normalized wavenumber to scan to. Reduce kEnd to speed up
%        the computation. ZGV points with kZGV >~ kEnd will not be found.
%      - Dk: (kEnd/50): Increment for the target k0 where to search for ZGV
%        points. Increase Dk to speed up the computation. Smaller Dk make it
%        more likely that all ZGV points will be found. 
%      - kStart: (2*Dk): Initial shift for normalized wavenumber k. Searching at
%        low wavenumbers is super slow, avoid it.
%      - DeltaPert: (1e-3): Regularization parameter. At ZGV points, we have two
%        very close solutions k and (1+DeltaPert)*k. These are computed as
%        initial guess where the ZGV point is finally computed.
%      - wmax: (inf, inactive) maximum angular frequency that defines kEnd
%              via the wave speed. This option is an extension for computeZGVScan().
%      - MaxIter: (100) Break search when reaching this number of iterations.
%      - show: (false) whether to display information during calculation 
%
% [1] E. Jarlebring, S. Kvaal, and W. Michiels, "Computing all Pairs (λ,μ) Such that 
% λ is a Double Eigenvalue of A+μB," SIAM J. Matrix Anal. Appl., vol. 32, no. 3, 
% pp. 902–927, Jul. 2011, doi: 10.1137/100783157.
%
% [2] D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 03, 2022. doi: 10.48550/arXiv.2211.01995.
%
% [3] K. Meerbergen, B. Plestenjak, "A Sylvester-Arnoldi type method for the 
% generalized eigenvalue problem with two-by-two operator determinants, 
% Numer. Linear Algebra Appl. 22 (2015) 1131-1146
%
% 2023 - Bor Plestenjak, adapted by Daniel A. Kiefer

narginchk(4, 5);
if nargin < 5, opts = []; end
if isfield(opts,'Neigs'),        Neigs = opts.Neigs;              else, Neigs = 6;          end
if isfield(opts,'kEnd'),         kEnd = opts.kEnd;                else, kEnd = 40;          end
if isfield(opts,'Dk'),           Dk = opts.Dk;                    else, Dk = kEnd/50;       end % minimum increment in kh after each iteration
if isfield(opts,'kStart'),       kStart = opts.kStart;            else, kStart = 2*Dk;      end % initial shift 
if isfield(opts,'DeltaPert'),    DeltaPert = opts.DeltaPert;      else, DeltaPert = 1e-3;   end % regularization parameter: leads to error if too small because sylvester has no unique solution
if isfield(opts,'wmax'),         wmax = opts.wmax;                else, wmax = inf;         end % maximum normalized frequency of interest
if isfield(opts,'MaxIter'),      MaxIter = opts.MaxIter;          else, MaxIter = 100;      end
if isfield(opts,'show'),         show = opts.show;                else, show = false;       end % display iteration results
if ~isfield(opts, 'hermitian'), opts.hermitian = [];        end
if isempty(opts.hermitian)
    if ishermitian(L2) && ishermitian(1i*L1) && ishermitian(L0) && ishermitian(M)
        opts.hermitian = true;
    else 
        opts.hermitian = false; 
    end
end

m = size(L0,1);
n = m^2;

found = 0;   % number of ZGV points found
k0 = kStart; % wavenumber shift close to which we compute candidates
kmin = 1e-4*kStart; % search candidates above this value
iter = 0;    % number of iterations
kzgv = []; wzgv = []; % ZGV wavenumbers and frequencies
k0s = nan(1,MaxIter);    % to collect used target values
optsNewton.beta_corr = true; optsNewton.show = false; optsNewton.maxsteps = 10; % options for ZGVNewtonBeta()
warningStatus = warning('query', 'MATLAB:eigs:NotAllEigsConverged');
warning('off', 'MATLAB:eigs:NotAllEigsConverged')
% k0list = kStart:Dk:kEnd;
while k0<kEnd && iter<MaxIter
    iter = iter + 1;
    % compute candidates:
    [z,Lbd] = ZGV_eigs(L2,L1,L0,M,DeltaPert,Neigs,1i*k0); % solves MEP for ZGV wavenumbers
    ks = k0 - 1i./diag(Lbd);
    ind = find(abs(imag(ks)./real(ks))<1e-3 & real(ks) > kmin); % candidates are solutions that are real
    % for each candidate in ks(ind), compute mu = w^2 and apply Newton iteration: 
    for j = ind.'
        mu = multDelta2(z(:,j))/multDelta0(z(:,j)); % Rayleigh quotient for mu
        if abs(imag(mu))<1e-2*abs(real(mu)) % refine solution using Newton-type method
            w = real(sqrt(mu)); % angular frequency (nondimensional)
            % [kR,wR,~,isConverged] = ZGVNewtonBeta(full(L2), full(L1), full(L0), full(M),...
            %     real(ks(j)), w, [], optsNewton);
            if opts.hermitian
                [kR,wR,~,isConverged] = ZGVNewtonBeta(full(L2), full(L1), full(L0), full(M),...
                real(ks(j)), w, [], optsNewton);
            else 
                [kR,wR,~,~,isConverged] = ZGVNewtonComplex(full(L2), full(L1), full(L0), full(M),...
                real(ks(j)), w, [], [], optsNewton);
            end
            % add converged solutions to list of ZGV points:
            if isConverged && kR>kmin && wR<wmax
                if isempty(kzgv)
                    kzgv = kR; wzgv = wR; found = 1;
                else % only add if not yet in list
                    difSquare = (kzgv - kR).^2 + (wzgv - wR).^2; % zgvs - [kR, wR];
                    minNormDif = sqrt(min(difSquare));
                    if minNormDif > 1e-8*sqrt(kR^2 + wR^2) % relative to norm [kR, wR]
                        kzgv = [kzgv; kR]; wzgv = [wzgv; wR]; found = found + 1;
                    end
                end
            end
        end
    end
    if show
        fprintf('%d. finished search at k0*h=%g: found %d ZGV points in total.\n', iter, k0, found);
    end
    
    % updated shift k0 and minimum kmin for next search: 
    if ~isempty(ks(ind)) % avoid error: max of empty array is an empty array
        kmin = max(real(ks(ind)))*0.95; % update kmin. Note that there can be serverl ZGV with similar k.
    end
    k0s(iter) = k0;
    k0 = max(k0 + Dk, kmin + Dk/2);
end
k0s = k0s(~isnan(k0s)); % crop nan entries

if (iter >= MaxIter)
    warning('ZGV_Sylv_MFRDScan:maxIterReached', ...
        'Maximum iteration number reached. If you are missing ZGV points, increase opts.MaxIter which is currently %d.', MaxIter);
end
if strcmp(warningStatus.state, 'on'), warning('on', 'MATLAB:eigs:NotAllEigsConverged'); end

[wzgv, ind] = sort(wzgv); % sort in frequency before returning
kzgv = kzgv(ind);

function x = multDelta0(y)
    Y1 = reshape(y(1:n),m,m);
    Y2 = reshape(y(n+1:2*n),m,m);
    X1 = M*(Y1*transpose(L1) + Y2*transpose(L2)) - ((1+DeltaPert)*L1*Y1 + (1+DeltaPert)^2*L2*Y2)*transpose(M);
    X2 = M*Y1*transpose(L2) - (1+DeltaPert)^2*L2*Y1*transpose(M);
    x = y'*[reshape(X1,m^2,1); reshape(X2,m^2,1)];
end

function x = multDelta2(y)
    Y1 = reshape(y(1:n),m,m);
    Y2 = reshape(y(n+1:2*n),m,m);
    X1 = ((1+DeltaPert)*L1*Y1+(1+DeltaPert)^2*L2*Y2)*transpose(L0) - L0*(Y1*transpose(L1)+Y2*transpose(L2));
    X2 = (1+DeltaPert)^2*L2*(Y1*transpose(L0)+Y2*transpose(L1)) - (L0*Y1 + (1+DeltaPert)*L1*Y2)*transpose(L2);
    x = y'*[reshape(X1,m^2,1); reshape(X2,m^2,1)];
end

end % ZGV_MFRDScan
