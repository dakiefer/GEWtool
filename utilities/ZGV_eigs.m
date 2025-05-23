function [X,D] = ZGV_eigs(L2,L1,L0,M,delta,Neigs,shift)
% ZGV_eigs - compute eigenvalues close to a target of a MEP.
% 
% 2023 - written by Bor Plestenjak for ZGV_Sylv_MFRDScan().

m = size(L0,1);
n = m^2;

if m > 70
    L2 = sparse(L2); L1 = sparse(L1); L0 = sparse(L0); M = sparse(M);
end

L = L0 + shift*L1 + shift^2*L2;
Ldel = L0 + (1+delta)*shift*L1 + (1+delta)^2*shift^2*L2;
LM = transpose(L)/transpose(M);    % = B in A*X + X*B = C
MLdel = M\Ldel;         % = A in A*X + X*B = C

[Q,R] = schur(full(-MLdel),'complex');
[U,S] = schur(full(LM),'complex');

F1 = transpose(L1+shift*L2)/transpose(M);
F2 = transpose(L2)/transpose(M);
F3 = M\((1+delta)*L1+shift*(1+delta)^2*L2);
F4 = (1+delta)^2*(M\L2);

F1=full(F1);
F2=full(F2);
F3=full(F3);
F4=full(F4);

[X,D] = eigs(@(y) multGamma(y,F1,F2,F3,F4,Q,R,U,S,m,n,shift),2*m^2,Neigs,'largestabs','Tolerance',1e-4);
end



% auxiliary function computes x = (Delta1 - shift*Delta0)\(Delta0*y)
% in an efficient way using transformation into a Sylvester equation
% where Delta1, Delta0 are as in scan_ZGV_modes
function x = multGamma(y,F1,F2,F3,F4,Q,R,U,S,m,n,shift)
Y1 = reshape(y(1:n),m,m);
Y2 = reshape(y(n+1:2*n),m,m);
C = Y1*F1 - F3*Y1 + Y2*F2 - F4*Y2;
Z1 = sylvester(R,S,-Q'*C*U);
X1 = Q*Z1*U';
X2 = Y1 + shift*X1;
x = [reshape(X1,m^2,1); reshape(X2,m^2,1)];
end