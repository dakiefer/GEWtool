function [k, w] = ZGVDirect(L2,L1,L0,M,opts)
% ZGVDirect - Compute ZGV points by solving the singular three-parameter EVP.
%
% Returns ZGV points (k, w) of the dispersion curves defined by 
% [ (i*k)^2*L2 + i*k*L1 + L0 + w^2*M ]*u = 0 . 
% A second eigenvalue problem describes the exceptional mode that appears at the
% ZGV points. These two coupled problems form a singular two-parameter eigenvalue
% problem that is quadratic in i*k. To avoid increasing the problem size when
% linearizing, a singular three-parameter eigenvalue problem is finally solved
% as described in [1].
% 
% Options in opts: 
%    - sc_steps: (2) number of staircase reduction steps.
%    - membtol: (1e-6) tolerance for residuum of ZGV condition. 
%    - showrank: (true) show rank of Delta matrices during reduction. 
%    - show: (true) show output 
%    - rrqr: (true) See singgep().
%    - method: ('project') method for solving the singular EVP.
%
% Literature:
% [1] D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, "Computing 
% zero-group-velocity points in anisotropic elastic waveguides: globally and locally 
% convergent methods." arXiv, Nov. 2022. doi: 10.48550/arXiv.2211.01995.
%
% 2022 - Bor Plestenjak, adapted by Daniel A. Kiefer

if nargin<5, opts=[]; end
class_t = superiorfloat(L0,L1,L2,M);
n = size(L0,2);
if ~exist('threepar_delta', 'file') % check if function is on Matlab path
    error('GEWtool:installMultiParEig',...
        'This functionality requires you to first install MultiParEig from\nhttps://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig');
end

% define matrices of three-parameter EVP:
Z = zeros(n, class_t);
L2d = [L2 Z; Z L2];
L1d = [L1 Z; 2*L2 L1];
L0d = [L0 Z; L1 L0];
Md  = [M Z; Z M];
C2  = [1 0; 0 0];
C1  = [0 -1;-1 0];
C0  = [0 0;0 1];

% assemble Delta-matrices (operator determinants):
[Delta0,Delta1,Delta2,Delta3] = threepar_delta(-C0,C1,0*C0,C2,-L0d,L1d,Md,L2d,-L0,L1,M,L2);

% define algorithm options to be used in the following:
if isfield(opts,'sc_steps'),  sc_steps = opts.sc_steps;  else, sc_steps=2;      end
if isfield(opts,'membtol'),   membtol = opts.membtol;    else, membtol = 1e-6;  end
if ~isfield(opts,'showrank'), opts.showrank = true;      end
if ~isfield(opts,'show'),     opts.show = true;          end
if ~isfield(opts,'rrqr'),     opts.rrqr=true;               end
if ~isfield(opts,'method'),   opts.method = 'project';   end

% optional initial staircase reduction of singular pencils
if sc_steps>0
    Delta = {Delta0,Delta1,Delta2,Delta3};
    for i=1:sc_steps
        Delta = staircase_step_cr_np(Delta,opts);
    end
    Delta0 = Delta{1}; Delta1 = Delta{2}; % we only need Delta0 and Delta1
end

% We solve a singular three-parameter eigenvalue problem by first applying singgep
% to the pencil (Delta1,Delta0) to obtain lambda = 1i*k. The frequencies w are then 
% obtained by inserting the lambda-values into the guided-wave problem. Solutions are
% only those where indeed cg=0 and k>0.
lambda = numeric_t([],class_t); % lambda = 1i*k
mu = numeric_t([],class_t);     % mu = w^2
lambdaCand = singgep(Delta1,Delta0,opts);
for i = 1:length(lambdaCand)
    L = lambdaCand(i)^2*L2 + lambdaCand(i)*L1 + L0;
    [X,muCand,Y] = eig(-L,M, 'vector'); % return square-frequencies muCand as vector
    for j = 1:length(muCand)
        iWd = 1i*2*lambdaCand(i)*L2 + 1i*L1; % note that 1i*Wd is Hermitean
        if ~isinf(muCand(j)) && abs(Y(:,j)'*iWd*X(:,j))/(norm(Y(:,j))*norm(X(:,j)))<membtol
           lambda = [lambda; lambdaCand(i)];
           mu = [mu; muCand(j)];
        end
    end
end

% extract meaningful solutions - lambda should be imaginary
isImag = imag(lambda)>1e-3 & ( abs(real(lambda)) < abs(imag(lambda))/1e4 );
ind = isImag & isfinite(mu);
k = real(-1i*lambda(ind));
w = sqrt(real(mu(ind)));

