function [S] = strain(gew,dat)
%STRAIN Compute the strain S.

S = cell(gew.geom.nLay, 1); % allocate for each layer
for i = 1:gew.geom.nLay
    % initialize 
    lay = gew.lay(i);
    D1 = 1/lay.h*lay.D1; % TODO: multilayer
    Nudof = gew.geom.Nudof(i); % Lamb and/or SH
    es = eye(Nudof); % unit directional vectors
    
    % compute du/dy:
    uu = permute(dat.u{i}, [1, 2, 5, 3, 4]); % additional dimension for mult. with diff. mat.
    ud = sum(shiftdim(D1, -2).*uu, 4);       % du/dy
    uu = permute(uu, [1,2,4,3,5]);           % same size as ud
    
    % compute displacement gradient
    es = shiftdim(es, -3);
    ex = es(:,:,:,:,1);   % propagation direction
    F = 1i*dat.k.*ex.*uu; % inialize displacement gradient 
    for j=2:Nudof % every dimension of waveguide cross section
        ej = es(:,:,:,:,j);
        F = F + ej.*ud;
    end
    
    % compute stress
    S{i} = 1/2*(F + permute(F, [1,2,3,5,4])); % symmetric part of gradient
end 

end % end function
