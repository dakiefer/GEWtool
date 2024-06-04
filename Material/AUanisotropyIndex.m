function AU = AUanisotropyIndex(mat, N)
% AUanisotropyIndex - Universal anisotropy index.
%
% Arguments:
% - mat:     (class Material) Material to be homogenized 
% - N:       (integer, default: 1200) Number of directions on which to average
% 
% Literature: 
% S. I. Ranganathan and M. Ostoja-Starzewski, “Universal Elastic Anisotropy
% Index,” Phys. Rev. Lett., vol. 101(5), Aug. 2008, doi: 10.1103/PhysRevLett.101.055504.
% 
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%        Gatien Clement, Institut Langevin, Paris, France
%        Claire Prada, Institut Langevin, Paris, France

if nargin < 2
    N = 1200; 
end

[voigt,reuss] = homogenizeUniform(mat, N); 
AU = sum(voigt.C.*inv(reuss.C),'all') - 6; % AU = Cv : Sr - 6

end % function 
