function B = rotateEuler(A, varargin)
% ROTATEEULER - Rotate the nth-order tensor "A".
%
% Arguments:
% - A:       (3x3x...x3 numeric) nth-order tensor (all dimensions 3)
% - angle:   (scalar numeric) angle in radian to rotate
% - axis:    (one of 'x', 'y', 'z') axis to rotate about
% You can specify as many sets of angle-axis pairs as you desire. The rotation
% will be performed in the provided order in the fixed initial coordinate system
% (extrinsic).
% 
% Alternative Arguments:
% - A:       (3x3x...x3 numeric) nth-order tensor (all dimensions 3)
% - angleX:  (scalar numeric) angle in rad to turn around axis x
% - angleY:  (scalar numeric) angle in rad to turn around axis y'
% - angleZ:  (scalar numeric) angle in rad to turn around axis z''
% This is an extrinsic, passive rotation around x-y-z (in that order). Thereby, 
% x-y-z is the original fixed coordinate system. Note that these are improper
% Euler angles (also called Cardan angles or Tait–Bryan angles).
% For more details see: https://en.wikipedia.org/wiki/Euler_angles
% 
% usage: 
% R = rotateEuler(A, angle, axis); 
% R = rotateEuler(A, angle1, axis1, angle2, axis2); % rotate around
%          "axis1", then around "axis2".
% R = rotateEuler(A, angleX, angleY, angleZ); % rotate x, then y, then z.
% 
% Literature: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    % some error checking:
    if isrow(A)
        error('GEWTOOL:transformBasis:rowVector', 'A should be a column vector.');
    end
    sA = size(A);
    if iscolumn(A), sA = sA(1); end % weird matlab size
    if ~all(sA == 3)
        error('GEWTOOL:rotateEuler:wrongSize', 'All dimensions of A should be of size 3.');
    end
    
    R = eulerAnglesToRotationMatrix(varargin{:});
    B = transformBasis(A, R);

end % function 
