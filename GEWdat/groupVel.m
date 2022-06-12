function [cg] = groupVel(gew, dat)
%GROUPVEL compute the group velocities of guided wave solutions.

M = shiftdim(gew.op.M,-2); % mass operator such that dimension are: [nF, nK, M]
L2 = shiftdim(gew.op.L2,-2); % stiffness op in k^2
L1 = shiftdim(gew.op.L1,-2); % stiffness op in k
w = dat.w/gew.np.fh0*gew.np.h0; % normalize frequency like the above operators
k = dat.k*gew.np.h0; % wavenumbers
u = zeros([size(dat.k), gew.geom.Ndof]); % allocate
gdofsAccum = [];
for l = 1:gew.geom.nLay % convert structured u into unstructured u
    gdofLay = gew.geom.gdofOfLay{l}; % where to put into the global u vector
    [gdofNew, ldofNew] = setdiff(gdofLay, gdofsAccum); % remove coincident nodes
    gdofsAccum = [gdofsAccum, gdofNew]; % remember already treated dofs
    ulay = reshape(dat.u{l}, size(k,1), size(k,2), []);
    u(:,:,gdofNew) = ulay(:,:,ldofNew);
end
v = conj(u); % left eigenvectors;
u = permute(u, [1 2 4 3]); % right eigenvectors (contract with second dim of M,L2,L1)
Q = sum(v.*sum((2*k.*L2 - 1i*L1).*u, 4), 3); % nominator
P = sum(v.*sum((2*w.*M).*u, 4), 3); % denominator 
cg = Q./P*gew.np.fh0;

end
