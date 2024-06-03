function R = eulerAnglesToRotationMatrix(varargin)
% eulerAnglesToRotationMatrix - Convert set of angles-axes to a rotation matrix.
%
% Arguments:
% - angle:   (scalar numeric) angle in radian to rotate
% - axis:    (one of 'x', 'y', 'z') axis to rotate about
% 
% You can specify as many sets of angle-axis pairs as you desire. The rotation
% will be performed in the provided order in the fixed initial coordinate system
% (extrinsic).
% 
% usage: 
% R = eulerAnglesToRotationMatrix(angle, axis); 
% R = eulerAnglesToRotationMatrix(angle1, axis1, angle2, axis2); % rotate around
%          "axis1", then around "axis2".
% 
% Literature: 
% [1] D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
% [2] W. L. Bond, “The mathematics of the physical properties of crystals,” The
% Bell System Technical Journal, vol. 22, no. 1, pp. 1–72, Jan. 1943, doi:
% 10.1002/j.1538-7305.1943.tb01304.x.
% [3] https://en.wikipedia.org/wiki/Euler_angles
% 
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% some error checking:
if length(varargin) == 3 
    if isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
        varargin(2,:) = {'x','y','z'};
    end
end
if mod(numel(varargin),2)
    error('GEWTOOL:rotateEuler:wrongArg',...
        'Provide always a combination of angle and axis of rotation, e.g.: pi, ''x''.');
end

angles = cell2mat(varargin(1:2:end));
axes = varargin(2:2:end);
Nrots = length(angles);

% generic rotation about z-axis (passive, i.e., the axes are rotated):
Z0 = @(phi) [cos(phi), sin(phi), 0; 
            -sin(phi), cos(phi), 0; 
            0,         0,        1];

R = eye(3);  % initialize rotation matrix
for i = 1:Nrots 
    switch axes{i}
        case 'x'
            ax = 1; 
        case 'y'
            ax = 2; 
        case 'z' 
            ax = 0; % 3 == 0 (don't shift matrix)
    end
    Ri = circshift(Z0(angles(i)), [ax, ax]); % rotation matrix around ax by angle i
    R = R*Ri;   % extrinsic, passive rotation axes{1} -> axes{2}, etc (in that order). 
                % Reversing matrix multiplication changes to intrinsic rotation.
end

end % function 
