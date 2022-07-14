function B = rotateEuler(A, a, b, g)
% ROTATEEULER - Rotate the nth-order tensor "A" sequentially around x-y-z.
% Extrinsic, passive rotation around x-y-z (in that order). Thereby, 
% x-y-z is the original fixed coordinate system.
% For more details see: https://en.wikipedia.org/wiki/Euler_angles
%
% Arguments:
% - A:       Tensor of nth-order (n-dimensional array).
% - a,b,g:   Euler/Cardan angles of rotation (around x-y-z, respectively).
% 
% Literature: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022 - Daniel Kiefer
% Institut Langevin, Paris, France

% some error checking:
if isrow(A)
    error('GEWTOOL:transformBasis:rowVector', 'A should be a column vector.');
end
sA = size(A);
if iscolumn(A), sA = sA(1); end % weird matlab size
if ~all(sA == 3)
    error('GEWTOOL:rotateEuler:wrongSize', 'All dimensions of A should be of size 3.');
end

% generic rotation about z-axis (passive, i.e., the axes are rotated):
Z0 = @(phi) [cos(phi), sin(phi), 0; 
            -sin(phi), cos(phi), 0; 
            0,         0,        1];

Z = Z0(g);                    % rotation matrix around z by angle g
X = circshift(Z0(a), [1, 1]); % rotation matrix around x by angle a
Y = circshift(Z0(b), [2, 2]); % rotation matrix around y by angle b

Q = X*Y*Z; % extrinsic, passive rotation around x-y-z (in that order)
B = transformBasis(A, Q);

end % function 
