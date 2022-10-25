function [k, w] = ZGV_TwoParamEVP(L0,L1,L2,M,opts)

% [lambda, mu] = ZGV_QEP(K0,K1,K2,M) finds points (lambda,mu) for the
% eigenvalue problem (K0 + lambda*K1 + lambda^2*K2 + mu*M)x = 0, 
% where the curve mu(lambda) is flat, i.e., derivative mu'(lambda)=0
%
% The problem is transformed into a singular two-parameter eigenvalue 
% problem and solved using the function singtwopar. 

% Keep the size of the problem small, i.e., matrices should not be larger
% than 20 x 20.

% Bor Plestenjak 2022

if nargin<5, opts=[]; end

class_t = superiorfloat(L0,L1,L2,M);

n = size(L0,2);
Z = zeros(n, class_t);

A1 = [0 0;0 1];
B1 = [0 -1;-1 0];
C1 = [0 0; 0 0];
D1 = [1 0; 0 0];

A3 = L0;
B3 = L1;
C3 = M;
D3 = L2;

A2 = [L0 Z; L1 L0];
B2 = [L1 Z; 2*L2 L1];
C2 = [M Z; Z M];
D2 = [L2 Z; Z L2];

[Delta0,Delta1,Delta2,Delta3] = threepar_delta(-A1,B1,C1,D1,-A2,B2,C2,D2,-A3,B3,C3,D3);

lambda = numeric_t([],class_t); 
mu = numeric_t([],class_t); 

if isfield(opts,'sc_steps'), sc_steps = opts.sc_steps; else, sc_steps=2;      end
if isfield(opts,'membtol'),  membtol = opts.membtol;   else, membtol = 1e-6;  end
if ~isfield(opts,'showrank'), opts.showrank = 1;    end
if ~isfield(opts,'rrqr'),      opts.rrqr=1;       end
if ~isfield(opts,'method'),  opts.method = 'project';  end

% optional initial staircase reduction of singular pencils
if sc_steps>0
    Delta = {Delta0,Delta1,Delta2,Delta3};
    for i=1:sc_steps
        Delta = staircase_step_cr_np(Delta,opts);
    end
    Delta0 = Delta{1}; Delta1 = Delta{2}; % we only need Delta0 and Delta1
end

% We solve a singular two-parameter eigenvalue problem by first applying
% singgep to the pencil (Delta1,Delta0) and then by inserting computed
% lambda's in the GEP and checking the ZGV criteria for the computed mu's. 
lambdaCand = singgep(Delta1,Delta0,opts);
for i = 1:length(lambdaCand)
    L = L0 + lambdaCand(i)*L1 + lambdaCand(i)^2*L2;
    [X,muCand,Y] = eig(L,-M, 'vector'); % return eigenvalues muCand as vector
    for j = 1:length(muCand)
        Ld = L1 + 2*lambdaCand(i)*L2;
        if ~isinf(muCand(j)) && abs(Y(:,j)'*Ld*X(:,j))/(norm(Y(:,j))*norm(X(:,j)))<membtol
           lambda = [lambda; lambdaCand(i)];
           mu = [mu; muCand(j)];
        end
    end
end

% extract meaningful solutions - lambda should be imaginary
ind = intersect(find(abs(real(lambda))<1e-3),find(imag(lambda)>1e-1));
ind = intersect(ind, find(isfinite(mu)));
k = real(-1i*lambda(ind));
w = sqrt(real(mu(ind)));

