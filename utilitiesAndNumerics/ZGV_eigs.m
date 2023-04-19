function [X,D] = ZGV_eigs(L2,L1,L0,M,delta,Neigs,shift)
% ZGV_eigs - compute eigenvalues close to a target of a MEP.
% 
% 2023 - written by Bor Plestenjak for ZGV_Sylv_MFRDScan().

m = size(L0,1);
n = m^2;

L = L0 + shift*L1 + shift^2*L2;
Ldel = L0 + (1+delta)*shift*L1 + (1+delta)^2*shift^2*L2;
LM = transpose(L)/M;    % = B in A*X + X*B = C
MLdel = M\Ldel;         % = A in A*X + X*B = C

[Q,R] = schur(-MLdel,'complex');
[U,S] = schur(LM,'complex');

F1 = transpose(L1+shift*L2)/M;
F2 = transpose(L2)/M;
F3 = M\((1+delta)*L1+shift*(1+delta)^2*L2);
F4 = (1+delta)^2*(M\L2);

[X,D] = eigs(@multGamma,2*m^2,Neigs);

% auxiliary function computes x = (Delta1 - shift*Delta0)\(Delta0*y)
% in an efficient way using transformation into a Sylvester equation
% where Delta1, Delta0 are as in scan_ZGV_modes
function x = multGamma(y)
    Y1 = reshape(y(1:n),m,m);
    Y2 = reshape(y(n+1:2*n),m,m);
    C = Y1*F1 - F3*Y1 + Y2*F2 - F4*Y2;
    Z1 = sylvester(R,S,-Q'*C*U);
    X1 = Q*Z1*U';
    X2 = Y1 + shift*X1;
    x = [reshape(X1,m^2,1); reshape(X2,m^2,1)];
end

end