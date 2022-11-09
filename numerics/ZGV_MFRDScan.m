function [k, w] = ZGV_MFRDScan(L2, L1, L0, M, opts)

% Returns ZGV points (kh, omega) with the smallest kh for a twoparameter pencil 
% (K0 + kh*K1 + kh^2*K2 + mu*C)x = 0, where mu = 2*pi*omega^2/f0
% It is important that matrix M is nonsingular amd that system is
% transformed in such way that kh and mu are real.
% Input:
%    - L0,L1,L2,M: matrices of the model (for larger than 50 x 50, use sparse matrices)
%    - f0: model parameter to compute omega from mu
%    - options in opts:
%         - Neigs: (8) how many eigenvalues we compute using eigs for one shift
%         - MaxPoints: (20) how many ZGV points we want
%         - MaxIter: (10) how many iterations with different shift
%         - ShiftFactor: (1.05), factor for shift increment for next iteration
%         - DeltaPert: (1e-6), factor for perturbation is 1+DeltaPert
%         - ZeroShift: (1): initial shift for kh
%
% We start scanning with kh = ZeroShift, compute candidates for the ZGV
% points using eigs on the corresponding Delta matrices from quadratic 
% two-parameter eigenvalue problem linearized as a three-parameter 
% eigenvalue problem. Each candidate is refined with Gauss-Newton method
%
% We use solutions of the original and the perturbed equation similar to
% the method by Jarlebring, Kvaal, Michiels
%
% 2022 - Bor Plestenjak, adapted by Daniel A. Kiefer

narginchk(4, 5);
if nargin < 5, opts = []; end

if isfield(opts,'Neigs'),           Neigs = opts.Neigs;                  else, Neigs = 8;          end
if isfield(opts,'MaxPoints'),       MaxPoints = opts.MaxPoints;          else, MaxPoints = 20;     end
if isfield(opts,'MaxIter'),         MaxIter = opts.MaxIter;              else, MaxIter = 10;       end
if isfield(opts,'ShiftFactor'),     ShiftFactor = opts.ShiftFactor;      else, ShiftFactor = 1.1;  end
if isfield(opts,'DeltaPert'),       DeltaPert = opts.DeltaPert;          else, DeltaPert = 1e-6;   end
if isfield(opts,'kStart'),          kStart = opts.kStart;                else, kStart = 0.5;    end
if ~isfield(opts,'show'),           opts.show = false;    end

% matrices for the linearization of the perturbed pair into a 3P eigenvalue problem
A1 = [0 0; 0 1];
B1 = [0 -1;-1 0];
C1 = [0 0; 0 0];
D1 = [1 0; 0 0];

A2 = L0;
B2 = L1;
C2 = M;
D2 = L2;

A3 = L0;
B3 = L1*(1+DeltaPert);
C3 = M;
D3 = L2*(1+DeltaPert)^2;

[Delta0,Delta1,Delta2,~] = threepar_delta(-A1,B1,C1,D1,-A2,B2,C2,D2,-A3,B3,C3,D3);

found = 0;   % number of ZGV points found
k0 = kStart; % wavenumber shift close to which we compute candidates
iter = 0;    % number of iterations
zgvs = [];  % ZGV points as [k, w]-pairs
opts.beta_corr = true; opts.show = false; opts.maxsteps = 10; % options for ZGVNewtonBeta()
warningStatus = warning('query', 'MATLAB:eigs:NotAllEigsConverged');
warning('off', 'MATLAB:eigs:NotAllEigsConverged')
while found<MaxPoints && iter<MaxIter
    iter = iter + 1;
    [z,lbd] = eigs(Delta1,Delta0,Neigs,1i*k0); % compute candidates
    ks = -1i*diag(lbd); % wavenumber candidates
    ind = find(abs(imag(ks))<1e-4 & real(ks) > kStart); % candidates are solutions that are real
    for j = ind.'
        mu = (z(:,j)'*Delta2*z(:,j))/(z(:,j)'*Delta0*z(:,j)); % Rayleight quotient estimation for mu
        if abs(imag(mu))<1e-3 % refine solution using Gauss-Newton method
            w = real(sqrt(mu));
            [kR,wR,~,isConverged] = ZGVNewtonBeta(full(L2), full(L1), full(L0), full(M), real(ks(j)), w, [], opts);
            if isConverged && (kR>kStart) % add to saved list of ZGV points
                if isempty(zgvs)
                    zgvs = [kR, wR]; found = 1;
                else % only add if not yet in list
                    dif = zgvs - [kR, wR];
                    minNormDif = sqrt(min(dif(:,1).^2 + dif(:,2).^2));
                    if minNormDif > 1e-8*sqrt(kR^2 + wR^2)
                        zgvs = [zgvs; [kR wR]];
                        found = found + 1;
                    end
                end
            end
        end
    end
    if k0 < 1e-1, k0 = 0.8; end   % avoid getting stuck for initial shift close to 0
    if ~isempty(ks(ind))    % max of empty array is an empty array
        k0 = max(k0,max(real(ks)))*ShiftFactor; 
    else
        k0 = k0*ShiftFactor;
    end
end
if strcmp(warningStatus.state, 'on'), warning('on', 'MATLAB:eigs:NotAllEigsConverged'); end

% return values:
if ~isempty(zgvs)
    k = real(zgvs(:,1)); w = real(zgvs(:,2));
else
    k = []; w = [];
end
