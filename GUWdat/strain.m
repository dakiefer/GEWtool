function [S] = strain(wguide,dat)
%STRAIN Compute the strain S.

% initialize 
D1 = 1/wguide.geom.h*wguide.lay.D1; % TODO: multilayer
Nudof = wguide.geom.Nudof; % Lamb and/or SH
es = eye(Nudof); % unit directional vectors

% compute du/dy:
uu = permute(dat.u, [1, 2, 5, 3, 4]); % assume this to be done already? 
ud = sum(shiftdim(D1, -2).*uu, 4);
uu = permute(uu, [1,2,4,3,5]);

% compute displacement gradient
es = shiftdim(es, -3);
ex = es(:,:,:,:,1);
F = 1i*dat.k.*ex.*uu; % inialize displacement gradient 
for j=2:Nudof
    ej = es(:,:,:,:,j);
    F = F + ej.*ud;
end

% compute stress
S = 1/2*(F + permute(F, [1,2,3,5,4]));

end
