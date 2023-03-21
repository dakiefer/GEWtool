function [kzgv, wzgv] = ZGV_Sylv_MFRDScan(L2, L1, L0, M, opts)
% ZGV_Sylv_MFRDScan - Compute ZGV points via iterative shift and search using fast eigenvalue solver.
%
% Returns ZGV points (k, w) with the smallest k for 
% [ (i*k)^2*L2 + i*k*L1 + L0 + mu*M ]*u = 0, where mu = w^2.
% The Method of Fixed Relative Distance (MFRD) proposed by Jarlebring, Kvall and
% Michiels [1] is used to this end. The method computes candidate wavenumbers close 
% to a shift k0 and refines them with a Newton-type iteration. By iteratively adapting
% the shift k0, a wavenumber range [kmin, kmax] is searched. Details can be found in [2].
%
% It is important that matrix M is nonsingular and that the system is
% transformed in such way that k and mu are real.
% Input:
%    - L2,L1,L0,M: matrices of the model 
%    - options in opts:
%         - Neigs: (8) number of eigenvalues to compute using eigs for one shift
%         - MaxPoints: (20) maximum number of requested ZGV points 
%         - MaxIter: (10) maximum number of iterations with different shift
%         - ShiftFactor: (1.1), factor for shift increment for next iteration
%         - DeltaPert: (1e-6), factor for perturbation is 1+DeltaPert
%         - kStart: (1): initial shift for normalized wavenumber k
%         - kMax: maximum normalized wavenumber to scan 
%         - show: (false) whether to display information during calculation 
%         - wmax: (inf, inactive) maximum angular frequency that defines kMax
%                 via the wave speed. This option is an extension for computeZGVScan().
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

if isfield(opts,'Neigs'),        Neigs = opts.Neigs;              else, Neigs = 8;          end
if isfield(opts,'MaxPoints'),    MaxPoints = opts.MaxPoints;      else, MaxPoints = 50;     end % number of ZGV points to find
if isfield(opts,'MaxIter'),      MaxIter = opts.MaxIter;          else, MaxIter = 30;       end
if isfield(opts,'ShiftFactor'),  ShiftFactor = opts.ShiftFactor;  else, ShiftFactor = 1.1;  end % relative increment for k0
if isfield(opts,'DeltaPert'),    DeltaPert = opts.DeltaPert;      else, DeltaPert = 1e-6;   end % regularization parameter 
if isfield(opts,'kStart'),       kStart = opts.kStart;            else, kStart = 1;         end % initial shift kStart, search will be done for k > kStart 
if isfield(opts,'kMax'),         kMax = opts.kMax;                else, kMax = inf;         end
if ~isfield(opts,'show'),        opts.show = false;                                         end % display iteration results

if opts.show, disp(opts), end

m = size(L0,1);
n = m^2;

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
    [z,Lbd] = ZGV_eigs(L2,L1,L0,M,DeltaPert,Neigs,1i*k0); % solves MEP for ZGV wavenumbers
    ks = k0 - 1i./diag(Lbd);
    ind = find(abs(imag(ks))<1e-4 & real(ks) > kmin); % candidates are solutions that are real
    % for each candidate in ks(ind), compute mu = w^2 and apply Newton iteration: 
    for j = ind.'
        mu = multDelta2(z(:,j))/multDelta0(z(:,j)); % Rayleigh quotient for mu
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

function x = multDelta0(y)
    Y1 = reshape(y(1:n),m,m);
    Y2 = reshape(y(n+1:2*n),m,m);
    X1 = M*(Y1*transpose(L1) + Y2*transpose(L2)) - ((1+DeltaPert)*L1*Y1 + (1+DeltaPert)^2*L2*Y2)*M;
    X2 = M*Y1*transpose(L2) - (1+DeltaPert)^2*L2*Y1*M;
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
