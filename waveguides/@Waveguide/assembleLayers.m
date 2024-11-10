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
np = obj.np;

% initialize:
L2 = zeros(geom.Ndof); L1 = zeros(geom.Ndof); 
L0 = zeros(geom.Ndof); M  = zeros(geom.Ndof);
for l = 1:geom.nLay
    % get operators of the layer l:
    lay = lays{l}; % layer l
    [L0lay, L1lay, L2lay] = lay.stiffnessOp(udof, np, lay.h, n);
    Mlay = lay.massOp(udof, np, lay.h);
    % assemble into global matrices:
    dof = geom.gdofOfLay{l}; % global degrees of freedom for layer l
    L2(dof,dof) = L2(dof,dof) + L2lay;
    L1(dof,dof) = L1(dof,dof) + L1lay;
    L0(dof,dof) = L0(dof,dof) + L0lay;
     M(dof,dof) =  M(dof,dof) +  Mlay;
end

% reset boundary conditions, if any have been set: 
obj.geom.gdofDBC = [];

% use electrically open BCs per default (the V-term in the paper):
if isa(lays{1}, 'LayerPlatePiezo') % at the bottom
    if ~isdiag(lays{1}.PP), error('TODO: only implemented for GLL Lagrange elements.'); end
    dof = obj.geom.gdofBC{1}; % boundary dofs of the first layer
    nBotPhi = dof(length(udof)+1, 1); % dof of electric potential at bottom of this layer
    L1(nBotPhi,nBotPhi) = L1(nBotPhi,nBotPhi) - 1i*lays{1}.mat.eps0/np.eps0;
end
if isa(lays{end}, 'LayerPlatePiezo') % at the top
    if ~isdiag(lays{l}.PP), error('TODO: only implemented for GLL Lagrange elements.'); end
    l = length(lays);
    dof = obj.geom.gdofBC{l}; % boundary dofs of the last layer
    nTopPhi = dof(length(udof)+1, 2); % dof of electric potential at top of this layer
    L1(nTopPhi,nTopPhi) = L1(nTopPhi,nTopPhi) - 1i*lays{l}.mat.eps0/np.eps0;
end

% assign result:
op.M = M; op.L0 = L0; op.L1 = L1; op.L2 = L2;
obj.op = op;

end % function
