function [cg] = groupVel(gew, dat)
%GROUPVEL compute the group velocities of guided wave solutions.

if isa(gew,"CylinderCircumferential")
    warning('Circumferential waves do not support this function yet. The results might be wrong.');
end
M = gew.op.M; % mass matrix 
L2 = gew.op.L2; % stiffness matrix in k^2
L1 = gew.op.L1; % stiffness matrix in k
w = dat.w/gew.np.fh0*gew.np.h0; % normalize frequency like the above operators
k = dat.k*gew.np.h0; % wavenumbers
u = eigenVecs(gew,dat);

% v = conj(u); % left eigenvectors;
% u = permute(u, [1 2 4 3]); % right eigenvectors (contract with second dim of M,L2,L1)
cg = zeros(size(k)); 
for i = 1:size(k,1)
    for j = 1:size(k,2)
        u0 = squeeze(u(i,j,:));
        cg(i, j) = (u0'*(2*k(i,j).*L2 - 1i*L1)*u0) ./ (2*w(i,j)*u0'*M*u0);
    end
end
% Q = sum(v.*sum((2*k.*L2 - 1i*L1).*u, 4), 3); % nominator
% P = sum(v.*sum((2*w.*M).*u, 4), 3); % denominator 
cg = cg*gew.np.fh0;

end
