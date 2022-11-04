function [k, w] = ZGV_MFRDScan(L0,L1,L2,M, opts)

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

% We start scanning with kh = ZeroShift, compute candidates for the ZGV
% points using eigs on the corresponding Delta matrices from quadratic 
% two-parameter eigenvalue problem linearized as a three-parameter 
% eigenvalue problem. Each candidate is refined with Gauss-Newton method

% We use solutions of the original and the perturbed equation similar to
% the method by Jarlebring, Kvaal, Michiels

% Bor Plestenjak 15.5.2022

narginchk(4, 5);
if nargin < 5, opts = []; end

if isfield(opts,'Neigs'),           Neigs = opts.Neigs;                  else, Neigs = 8;          end
if isfield(opts,'MaxPoints'),       MaxPoints = opts.MaxPoints;          else, MaxPoints = 20;     end
if isfield(opts,'MaxIter'),         MaxIter = opts.MaxIter;              else, MaxIter = 10;       end
if isfield(opts,'ShiftFactor'),     ShiftFactor = opts.ShiftFactor;      else, ShiftFactor = 1.1; end
if isfield(opts,'DeltaPert'),       DeltaPert = opts.DeltaPert;          else, DeltaPert = 1e-6;   end
if isfield(opts,'ZeroShift'),       ZeroShift = opts.ZeroShift;          else, ZeroShift = 0.5;      end
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

found = 0;
shift = ZeroShift;
iter = 0;
pairs = [];
opts.beta_corr = true; opts.show = false; opts.maxsteps = 10;
warningStatus = warning('query', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix')
while found<MaxPoints && iter<MaxIter
    iter = iter + 1;
    [X,D] = eigs(Delta1,Delta0,Neigs,1i*shift);
    ks = -1i*diag(D);
    ind = find(abs(imag(ks))<1e-4 & real(ks) > 1e-6); % candidates are solutions that are real and positive
    ind = intersect(ind, find(real(ks)>ZeroShift)); % candidates have larger wavenumber than what we have already found
    for j = ind'
        k = real(ks(j));
        mu = (X(:,j)'*Delta2*X(:,j))/(X(:,j)'*Delta0*X(:,j)); % Rayleight quotient estimation for mu
        if abs(imag(mu))<1e-3 
            % we refine and verify the solution using Gauss-Newton method
            w = real(sqrt(mu));
            [kR,wR,~,flag] = ZGVNewtonBeta(full(L2), full(L1), full(L1), full(M), k, w, [], opts);
            muR = wR^2;
            if (flag==true) && (kR>ZeroShift)
                dif = 1e10*norm([kR muR]);
                if ~isempty(pairs)
                    tmp = pairs - [kR muR];
                    dif = sqrt(min(abs(tmp(:,1)).^2+abs(tmp(:,1)).^2));
                end
                % we add point if it is not already in the set    
                if dif>1e-8*norm([kR muR])
                    pairs = [pairs; [kR wR]];
                    found = found + 1;
%                     fprintf('iter %3d | shift %10.3e,  ZGV point %3d  lambda (%16.9e,%16.9e) | mu (%16.9e,%16.9e) \n',iter,shift,found,real(kR),imag(kR),real(wR),imag(wR));
                end
            end
        end
    end
    if shift < 1e-1, shift = 0.8; end   % avoid getting stuck for initial shift close to 0
    if ~isempty(ks(ind))    % max of empty array is an empty array
        shift = max(shift,max(real(ks)))*ShiftFactor; 
    else
        shift = shift*ShiftFactor;
    end
end
if strcmp(warningStatus.state, 'on'), warning('on', 'MATLAB:nearlySingularMatrix'); end

k = real(pairs(:,1));
w = real(pairs(:,2));
