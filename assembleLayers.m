function [Mglob, L0glob, L1glob, L2glob] = assembleLayers(geom, lays, n, c0, h0, rho0)
udof = 1:3;

L2glob = []; L1glob = []; L0glob = []; Mglob = [];
for lay=lays
    [L0lay, L1lay, L2lay] = lay.stiffnessOp(udof, n);
    Mlay = lay.massOp(udof);
    c0l = lay.mat.tensor(1,2,1,2); rho0l = lay.mat.rho; hl = lay.h;
    L2glob = blkdiag(L2glob, L2lay*c0l/c0);
    L1glob = blkdiag(L1glob, L1lay*c0l/c0*h0/hl);
    L0glob = blkdiag(L0glob, L0lay*c0l/c0*(h0/hl)^2);
    Mglob =  blkdiag(Mglob, Mlay*rho0l/rho0);
end

Nitf = geom.nNodes - 2; % inner interfaces (no outer boundaries)
for ii = 1:Nitf
    lay1 = lays(ii); lay2 = lays(ii+1);
    ldofBC1 = geom.ldofBC{ii}(2,:); % top interface dofs
    ldofBC2 = geom.ldofBC{ii+1}(1,:); % bottom interface dofs
    gdofBC = [geom.gdofOfElem{ii}(ldofBC1), geom.gdofOfElem{ii+1}(ldofBC2)];
    gdofs = [geom.gdofOfElem{ii}(:); geom.gdofOfElem{ii+1}(:)];
    ldofBC1 = 2:2:2*length(udof); ldofBC2 = 1:2:2*length(udof); % redefine to fit the reduced BC matrices 
    c0l1 = lay1.mat.tensor(1,2,1,2); c0l2 = lay2.mat.tensor(1,2,1,2); 
    [B0l1, B1l1] = lay1.tractionOp(udof, n); 
    B0l1 = B0l1(ldofBC1, :)*c0l1/c0*h0/lay1.h; B1l1 = B1l1(ldofBC1, :)*c0l1/c0;
    [B0l2, B1l2] = lay2.tractionOp(udof, n);
    B0l2 = B0l2(ldofBC2, :)*c0l2/c0*h0/lay2.h; B1l2 = B1l2(ldofBC2, :)*c0l2/c0;
    [Ul1] = lay1.displacementOp(udof);
    [Ul2] = lay2.displacementOp(udof);
    CC0 = [B0l1,              -B0l2;  % traction
            Ul1(ldofBC1, :),  -Ul2(ldofBC2, :)]; % displ
    CC1 = [B1l1,              -B1l2;  % traction
         0*B1l1,              0*B1l2]; % displ
    L0glob(gdofBC, gdofs) = CC0; L1glob(gdofBC, gdofs) = CC1; 
    L2glob(gdofBC, gdofs) = 0; Mglob(gdofBC, gdofs) = 0;
end

end
