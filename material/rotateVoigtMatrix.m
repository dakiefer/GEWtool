function Crot = rotateVoigtMatrix(C, varargin)
% rotateVoigtMatrix - Rotate the 6x6 Voigt matrix C around the specified Euler angles.
%
% Arguments:
% - C:       (6x6 numeric) Voigt notated stiffness matrix
% - angle:   (scalar numeric) angle in radian to rotate
% - axis:    (one of 'x', 'y', 'z') axis to rotate about
% You can specify as many sets of angle-axis pairs as you desire. The rotation
% will be performed in the provided order in the fixed initial coordinate system
% (extrinsic).
% 
% See also: 
% [1] https://en.wikipedia.org/wiki/Euler_angles
% [2] B. A. Auld, Acoustic Fields and Waves in Solids, 2nd ed., vol. 1. Malabar,
% Fla: Krieger Publishing Company, 1990.
% [3] W. L. Bond, “The mathematics of the physical properties of crystals,”
% The Bell System Technical Journal, vol. 22, no. 1, pp. 1–72, Jan. 1943, doi:
% 10.1002/j.1538-7305.1943.tb01304.x.
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%        Gatien Clement, Institut Langevin, Paris, France

    % some error checking:
    validateattributes(C,{'numeric'},{'size',[6 6]});
    
    R = eulerAnglesToRotationMatrix(varargin{:});
    M = toVoigtRotation(R);
    Crot = M*C*M.'; 

end % function 


function M = toVoigtRotation(a)
% toVoigtRotation - Convert 3x3 rotation matrix to 6x6 rotation matrix.
% The 3x3 rotation matrix "a" is converted to the 6x6 rotation matrix "M" such
% that the Voigt notated stiffness matrix "C" can be rotated to the new axes by
% Cnew = M*C*M'; 
% 
% Literature: 
% [1] B. A. Auld, Acoustic Fields and Waves in Solids, 2nd ed., vol. 1. Malabar,
% Fla: Krieger Publishing Company, 1990.
% [2] W. L. Bond, “The mathematics of the physical properties of crystals,”
% The Bell System Technical Journal, vol. 22, no. 1, pp. 1–72, Jan. 1943, doi:
% 10.1002/j.1538-7305.1943.tb01304.x.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    M = zeros(6);
    p0 = [1,2,3]; p1 = [3,1,2]; p2 = [2,3,1];
    M(1:3,1:3) = a.^2; % (1,1)-block
    M(4:6,1:3) = a(p2,p0).*a(p1,p0); % (2,1)-block
    M(1:3,4:6) = 2*a(p0,p2).*a(p0,p1); % (1,2)-block
    M(4:6,4:6) = a(p1,p1).*a(p2,p2) + a(p1,p2).*a(p2,p1); % (2,2)-block
    
end
