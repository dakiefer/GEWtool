function [M, L0, L1, L2] = freeBCs(M, L0, L1, L2, geom, lays, n, c0, h0)
    udof = 1:3;

    lay1 = lays(1); lay2 = lays(end);
    ldofBC1 = geom.ldofBC{1}(1,:); % bottom interface dofs
    ldofBC2 = geom.ldofBC{end}(2,:); % top interface dofs
    gdofBC1 = geom.gdofOfElem{1}(ldofBC1);
    gdofBC2 = geom.gdofOfElem{end}(ldofBC2);
    gdofs1 = geom.gdofOfElem{1}(:);
    gdofs2 = geom.gdofOfElem{end}(:);
    ldofBC1 = 1:2:2*length(udof); ldofBC2 = 2:2:2*length(udof); % redefine to fit the reduced BC matrices 
    c0l1 = lay1.mat.tensor(1,2,1,2); c0l2 = lay2.mat.tensor(1,2,1,2); 
    [B0l1, B1l1] = lay1.tractionOp(udof, n); 
    B0l1 = B0l1(ldofBC1, :)*c0l1/c0*h0/lay1.h; B1l1 = B1l1(ldofBC1, :)*c0l1/c0;
    [B0l2, B1l2] = lay2.tractionOp(udof, n);
    B0l2 = B0l2(ldofBC2, :)*c0l2/c0*h0/lay2.h; B1l2 = B1l2(ldofBC2, :)*c0l2/c0;
    % bottom BCs:
    L0(gdofBC1, gdofs1) = B0l1; L1(gdofBC1, gdofs1) = B1l1; 
    L2(gdofBC1, gdofs1) = 0; M(gdofBC1, gdofs1) = 0;
    % top BCs:
    L0(gdofBC2, gdofs2) = B0l2; L1(gdofBC2, gdofs2) = B1l2; 
    L2(gdofBC2, gdofs2) = 0; M(gdofBC2, gdofs2) = 0;
end
