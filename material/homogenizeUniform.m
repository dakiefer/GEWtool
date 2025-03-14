function [voigt,reuss] = homogenizeUniform(mat, N)
% HOMOGENIZEUNIFORM - Avarage stiffness with a uniform orientation distribution.
% Mainly used to compute the universal anisotropy index AU.
%
% Arguments:
% - mat:     (class Material) Material to be homogenized 
% - N:       (scalar integer, default: 1200) Number of orientations for averaging.
% 
% See also
% AUanisotropyIndex
% 
% Literature
% S. Hirsekorn, “Elastic properties of polycrystals: a review,” Textures and
% Microstructures, vol. 12, no. 1–3, pp. 1–14, 1990.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%        Gatien Clement, Institut Langevin, Paris, France
%        Claire Prada, Institut Langevin, Paris, France

if nargin < 2
    N = 1200; 
end

Cnorm = norm(mat.C); 
C0 = mat.C/Cnorm; 
pointsXYZ = sphereLattice(4, N); % nearly equidistant distribution of points on 4d-sphere in Cartesian coordinates
[theta,phi,psi] = quaternionToEulerAngles(pointsXYZ,1,3,1); % Euler angles for x-z-x-rotation

% uniform Voigt and Reuss homogenization: 
Cvoigt = zeros(6,6); Sreuss = zeros(6,6);
for i = 1:N
    Ci = rotateVoigtMatrix(C0, theta(i),'x', phi(i),'z', psi(i),'x'); % faster than tensor rot, use same as above!
    Cvoigt = Cvoigt + Ci;
    Sreuss = Sreuss + inv(Ci);
end
Cvoigt = Cvoigt/N; Sreuss = Sreuss/N;
Creuss = inv(Sreuss); % invert close to unity (divide by N before)!

if norm(Cvoigt - Cvoigt.')/2 > 1e4*eps || norm(Creuss - Creuss.')/2 > 1e4*eps
    error('GEWTOOL:Material:homogenize', 'The homogenized tensors are not symmetric.');
end

voigt =  Material('Voigt', (Cvoigt + Cvoigt.')/2*Cnorm, mat.rho);
reuss =  Material('Reuss', (Creuss + Creuss.')/2*Cnorm, mat.rho); 

end % function 
