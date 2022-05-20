function op = assembleLayers(obj, udof, n)

% geometry and material:
geom = obj.geom; lays = obj.lay;
c0 = obj.np.c0; rho0 = obj.np.rho0; h0 = obj.np.h0; 

% initialize:
L2 = zeros(geom.Ndof); L1 = zeros(geom.Ndof); 
L0 = zeros(geom.Ndof); M  = zeros(geom.Ndof);
for l = 1:geom.nLay
    lay = lays(l); % layer l
    dof = geom.gdofOfLay{l}; % global degrees of freedom for layer l
    [L0lay, L1lay, L2lay] = lay.stiffnessOp(udof, n);
    Mlay = lay.massOp(udof);
    c0l = lay.mat.c(1,2,1,2); rho0l = lay.mat.rho; hl = lay.h; % normalization params
    % assemble into global matrices:
    L2(dof,dof) = L2(dof,dof) + L2lay*c0l/c0;
    L1(dof,dof) = L1(dof,dof) + L1lay*c0l/c0*h0/hl;
    L0(dof,dof) = L0(dof,dof) + L0lay*c0l/c0*(h0/hl)^2;
     M(dof,dof) =  M(dof,dof) +  Mlay*rho0l/rho0;
end

% return as structure:
op.M = M; op.L0 = L0; op.L1 = L1; op.L2 = L2;

end % function
