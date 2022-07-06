function [op] = freeBCs(obj, udof, n)
% freeBCs - Implement traction-free boundary conditions. This is a reminecent of
% the spectral collocation implementation.

geom = obj.geom; lays = obj.lay;
op = obj.op;
c0 = obj.np.c0; h0 = obj.np.h0; 

lay1 = lays(1); lay2 = lays(end);
ldofBC1 = geom.ldofBC{1}(1,:); % bottom interface dofs
ldofBC2 = geom.ldofBC{end}(2,:); % top interface dofs
gdofBC1 = geom.gdofOfLay{1}(ldofBC1);
gdofBC2 = geom.gdofOfLay{end}(ldofBC2);
gdofs1 = geom.gdofOfLay{1}(:);
gdofs2 = geom.gdofOfLay{end}(:);
ldofBC1 = 1:2:2*length(udof); ldofBC2 = 2:2:2*length(udof); % redefine to fit the reduced BC matrices 
c0l1 = lay1.mat.c(1,2,1,2); c0l2 = lay2.mat.c(1,2,1,2); 
[B0l1, B1l1] = lay1.tractionOp(udof, n); 
B0l1 = B0l1(ldofBC1, :)*c0l1/c0*h0/lay1.h; B1l1 = B1l1(ldofBC1, :)*c0l1/c0;
[B0l2, B1l2] = lay2.tractionOp(udof, n);
B0l2 = B0l2(ldofBC2, :)*c0l2/c0*h0/lay2.h; B1l2 = B1l2(ldofBC2, :)*c0l2/c0;
% bottom BCs:
op.L0(gdofBC1, gdofs1) = B0l1; 
op.L1(gdofBC1, gdofs1) = B1l1; 
op.L2(gdofBC1, gdofs1) = 0; 
op.M(gdofBC1, gdofs1) = 0;
% top BCs:
op.L0(gdofBC2, gdofs2) = B0l2; 
op.L1(gdofBC2, gdofs2) = B1l2; 
op.L2(gdofBC2, gdofs2) = 0; 
op.M(gdofBC2, gdofs2) = 0;

% % save into op for implementing alternative BCs: (Not being used at the time)
% % Carefule: layers are of different size -> do not simply concatenate
% op.B0 = zeros(size(op.M));
% op.B0(gdofBC1, gdofs1) = B0l1;
% op.B0(gdofBC2, gdofs2) = B0l2;
% op.B1 = zeros(size(op.M));
% op.B0(gdofBC1, gdofs1) = B1l1;
% op.B0(gdofBC2, gdofs2) = B1l2;
% TODO row order needs to match gdofBC! maybe separate into top and bottom BC? 

end % function
