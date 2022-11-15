function [kzgv, wzgv] = ZGV_MFRDScan(L2, L1, L0, M, opts)
% ZGV_MFRDScan - Compute ZGV points via iterative shift and search.
%
% Returns ZGV points (k, w) with the smallest k for 
% [ (i*k)^2*L2 + i*k*L1 + L0 + mu*M ]*u = 0, where mu = w^2.
% The Method of Fixed Relative Distance (MFRD) proposed by Jarlebring, Kvall and
% Michiels [1] is used for this end. The method computes candidate wavenumbers close 
% to a shift k0 and refines them with a Newton-type iteration. By iteratively adapting
% the shift k0, a wavenumber range [kmin, kmax] is searched. Details can be found in [2].
%
% It is important that matrix M is nonsingular amd that system is
% transformed in such way that k and mu are real.
% Input:
%    - L2,L1,L0,M: matrices of the model (for larger than 50 x 50, use sparse matrices)
%    - options in opts:
%         - Neigs: (8) how many eigenvalues we compute using eigs for one shift
%         - MaxPoints: (20) maximum number ofZGV points we want
%         - MaxIter: (10) maximum number of iterations with different shift
%         - ShiftFactor: (1.1), factor for shift increment for next iteration
%         - DeltaPert: (1e-6), factor for perturbation is 1+DeltaPert
%         - kStart: (1): initial shift for normalized wavenumber k
%         - kMax: maximum normalized wavenumber to scan within 
%         - show: (false) display information during calculation 
%
% [1] E. Jarlebring, S. Kvaal, and W. Michiels, "Computing all Pairs (λ,μ) Such that 
% λ is a Double Eigenvalue of A+μB," SIAM J. Matrix Anal. Appl., vol. 32, no. 3, 
% pp. 902–927, Jul. 2011, doi: 10.1137/100783157.
%
% [2] D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 03, 2022. doi: 10.48550/arXiv.2211.01995.
%
% 2022 - Bor Plestenjak, adapted by Daniel A. Kiefer

narginchk(4, 5);
if nargin < 5, opts = []; end

if isfield(opts,'Neigs'),        Neigs = opts.Neigs;              else, Neigs = 8;          end
if isfield(opts,'MaxPoints'),    MaxPoints = opts.MaxPoints;      else, MaxPoints = 20;     end
if isfield(opts,'MaxIter'),      MaxIter = opts.MaxIter;          else, MaxIter = 10;       end
if isfield(opts,'ShiftFactor'),  ShiftFactor = opts.ShiftFactor;  else, ShiftFactor = 1.1;  end
if isfield(opts,'DeltaPert'),    DeltaPert = opts.DeltaPert;      else, DeltaPert = 1e-6;   end
if isfield(opts,'kStart'),       kStart = opts.kStart;            else, kStart = 1;    end
if isfield(opts,'kMax'),         kMax = opts.kMax;                else, kMax = inf;    end
if ~isfield(opts,'show'),        opts.show = false;    end

if opts.show, disp(opts), end

% matrices for the three-parameter eigenvalue problem:
L2p = L2*(1+DeltaPert)^2;
L1p = L1*(1+DeltaPert);
C2 = [1 0; 0 0];
C1 = [0 -1;-1 0];
C0 = [0 0; 0 1];
% assemble corresponding Delta-matrices (operator determinants):
Delta0 = deltaMatrix(L2, L1, M, L2p, L1p, M, C2, C1, 0*C0);
Delta1 = -1*deltaMatrix(L2, L0, M, L2p, L0, M, C2, C0, 0*C0);
DeltaM = -1*deltaMatrix(L2, L1, L0, L2p, L1p, L0, C2, C1, C0);

found = 0;   % number of ZGV points found
k0 = kStart; % wavenumber shift close to which we compute candidates
kmin = 1e-8; % search candidates above this value
iter = 0;    % number of iterations
kzgv = []; wzgv = []; % ZGV wavenumbers and frequencies
optsNewton.beta_corr = true; optsNewton.show = false; optsNewton.maxsteps = 10; % options for ZGVNewtonBeta()
warningStatus = warning('query', 'MATLAB:eigs:NotAllEigsConverged');
warning('off', 'MATLAB:eigs:NotAllEigsConverged')
while k0<kMax && found<MaxPoints && iter<MaxIter
    iter = iter + 1;
    % compute candidates:
    [z,lbd] = eigs(Delta1,Delta0,Neigs,1i*k0); % compute candidates
    ks = -1i*diag(lbd); % wavenumber candidates (nondimensions)
    ind = find(abs(imag(ks))<1e-4 & real(ks) > kmin); % candidates are solutions that are real
    % for each candidate in ks(ind), compute mu = w^2 and apply Newton iteration: 
    for j = ind.'
        mu = (z(:,j)'*DeltaM*z(:,j))/(z(:,j)'*Delta0*z(:,j)); % Rayleight quotient estimation for mu
        if abs(imag(mu))<1e-2*abs(real(mu)) % refine solution using Newton-type method
            w = real(sqrt(mu)); % angular frequency (nondimensional)
            [kR,wR,~,isConverged] = ZGVNewtonBeta(full(L2), full(L1), full(L0), full(M),...
                real(ks(j)), w, [], optsNewton);
            % add converged solutions to list of ZGV points:
            if isConverged && (kR>kmin) 
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
    if opts.show
        fprintf('%d. finished search at k0*h=%g: found %d ZGV points in total.\n', iter, k0, found);
    end
    % updated shift k0 and minimum kmin for next search: 
    if k0 < 1e-1, k0 = 0.8; end   % avoid getting stuck for initial shift close to 0
    if ~isempty(ks(ind)) % avoid error: max of empty array is an empty array
        kmin = max(real(ks(ind))); % update kmin
        k0 = max(k0,kmin)*ShiftFactor; % update shift
    else
        if opts.show, fprintf('    no candidates.\n'); end
        k0 = k0*ShiftFactor; % update shift
    end
end
if strcmp(warningStatus.state, 'on'), warning('on', 'MATLAB:eigs:NotAllEigsConverged'); end

end % ZGV_MFRDScan

% helper function:
function Delta = deltaMatrix(A2, A1, A0, B2, B1, B0, C2, C1, C0)
    Delta = kron(A2, kron(B1, C0)) + kron(A1, kron(B0, C2)) + kron(A0, kron(B2, C1))...
        - kron(A0, kron(B1, C2)) - kron(A2, kron(B0, C1)) - kron(A1, kron(B2, C0));
    % when C0 = 0-matrix, this might be faster:
    % Delta0 = kron(L1, kron(M,C2)) + kron(M, kron(L2p,C1))...
    %     - kron(M, kron(L1p, C2)) - kron(L2, kron(M, C1));
    % Delta1 = -1*( kron(L0, kron(M, C2)) + kron(M, kron(L2p, C0))...
    %     - kron(M, kron(L0, C2)) - kron(L2, kron(M, C0)) );
end
