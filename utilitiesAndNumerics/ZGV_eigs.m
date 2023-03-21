function [X,D] = ZGV_eigs(L2,L1,L0,M,epsilon,Neigs,shift)
% ZGV_eigs - compute eigenvalues close to a target of a MEP.
% 
% 2023 - written by Bor Plestenjak for ZGV_Sylv_MFRDScan().

m = size(L0,1);
n = m^2;

S = L0 + shift*L1 + shift^2*L2;
Seps = L0 + (1+epsilon)*shift*L1 + (1+epsilon)^2*shift^2*L2;
S1 = transpose(S)/M;
S1eps = M\Seps;

[Q1,R1] = schur(-S1eps,'complex');
[Q2,R2] = schur(S1,'complex');

F1 = transpose(L1+shift*L2)/M;
F2 = transpose(L2)/M;
F3 = M\((1+epsilon)*L1+shift*(1+epsilon)^2*L2);
F4 = (1+epsilon)^2*(M\L2);

[X,D] = eigs(@multGamma,2*m^2,Neigs);

% auxiliary function computes x = (Delta1 - shift*Delta0)\(Delta0*y)
% in an efficient way using transformation into a Sylvester equation
% where Delta1, Delta0 are as in scan_ZGV_modes
function x = multGamma(y)
    Y1 = reshape(y(1:n),m,m);
    Y2 = reshape(y(n+1:2*n),m,m);
    C = Y1*F1 - F3*Y1 + Y2*F2 - F4*Y2;
    Z1 = sylvester(R1,R2,-Q1'*C*Q2);
    X1 = Q1*Z1*Q2';
    X2 = Y1 + shift*X1;
    x = [reshape(X1,m^2,1); reshape(X2,m^2,1)];
end

end