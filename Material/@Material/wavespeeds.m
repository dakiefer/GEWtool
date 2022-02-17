function [cs, eu] = wavespeeds(obj, ek)
    % WAVESPEEDS Compute the wave speeds by solving the Kelvin-Christoffel
    % equation (eigenvalue problem for cl, ct1, ct2).
    if nargin < 2, ek = [1; 0; 0]; end
    validateattributes(ek, {'numeric'}, {'vector', 'numel', 3});
    ek = ek(:)/norm(ek(:)); % normalized column vector
    ek = ek(:); ekS = shiftdim(ek, -3); % prepare for contraction
    D = 1/obj.rho*sum(sum(ek.*obj.c.*ekS, 4), 1); % Kristoffel tensor: 1/rho ek.c.ek
    D = squeeze(D); % 3x3 matrix
    [eu, cs2] = eig(D, 'vector'); 
    % cosTheta = abs(sum(ek.*eu, 1)); % angle with propagation direction
    % [~, ind] = sort(cosTheta, 'descend'); % sort: quasi-long., quasi-transv1, quasi-transv2
    % ^ the above might lead to "confiution" of modes (see slowness curves)
    [~, ind] = sort(cs2, 'descend');
    cs = sqrt(cs2(ind)); % wave speeds [cl, ct1, ct2]
    eu = eu(:,ind); % polarization vectors
end
