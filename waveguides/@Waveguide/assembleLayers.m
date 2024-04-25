function obj = assembleLayers(obj, udof, n)
% assembleLayers - Assembles the layer matrices into one global system.
% Arguments:
% - udof: vector specifying the displacement degrees of freedom (polarization)
% - n:    circumferential order for cylindrical waveguide
% This function modifies the argument "obj" passed by reference.
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% geometry and material:
geom = obj.geom; lays = obj.lay;
c0 = obj.np.c0; rho0 = obj.np.rho0; h0 = obj.np.h0; 

% initialize:
L2 = zeros(geom.Ndof); L1 = zeros(geom.Ndof); 
L0 = zeros(geom.Ndof); M  = zeros(geom.Ndof);
for l = 1:geom.nLay
    % get operators of the layer l:
    lay = lays(l); % layer l
    hl = lay.h/h0; % normalized layer thickness
    [L0lay, L1lay, L2lay] = lay.stiffnessOp(udof, hl, n);
    Mlay = lay.massOp(udof, hl);
    % assemble into global matrices:
    dof = geom.gdofOfLay{l}; % global degrees of freedom for layer l
    cl = lay.mat.c(1,2,1,2)/c0; rhol = lay.mat.rho/rho0; % normalized material parameters
    L2(dof,dof) = L2(dof,dof) + L2lay*cl;
    L1(dof,dof) = L1(dof,dof) + L1lay*cl;
    L0(dof,dof) = L0(dof,dof) + L0lay*cl;
     M(dof,dof) =  M(dof,dof) +  Mlay*rhol;
end

% reset boundary conditions, if any have been set: 
obj.geom.gdofDBC = [];

% assign result:
op.M = M; op.L0 = L0; op.L1 = L1; op.L2 = L2;
obj.op = op;

end % function
