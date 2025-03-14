function B = transformBasis(A, Q)
% TRANSFORMBASIS - Change of basis "Q" applied to the nth-order tensor "A".
% 
% Arguments:
% - A:   Tensor of nth-order (n-dimensional array).
% - Q:   Orthogonal matrix discribing the transformation, e.g., a rotation.
% 
% Literature: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

% some error checking:
sQ = size(Q); spaceDim = sQ(1); % space dimensionality, e.g., 3d or 2d
if ~ismatrix(Q) || isvector(Q) || sQ(1) ~= sQ(2)
    error('GEWTOOL:transformBasis:invalidTransformationMatrix', 'Q should be a square matrix.');
end
if mod(numel(A), spaceDim) ~= 0
    error('GEWTOOL:transformBasis:wrongSize', 'All tensor dimension of A must be of the same length and compatible with Q.')
end
if isrow(A)
    error('GEWTOOL:transformBasis:rowVector', 'A should be a column vector.');
end

sA = size(A);
if iscolumn(A), sA = sA(1); end % correct weird matlab size
n = length(sA); % order of tensor

% transform:
for dim = 1:n % transform each of the basis vectors of the tensor A one by one
    A = permute(A, [1:dim-1, n+1, dim:n]); % add one dimension before "dim"
    QQ = shiftdim(Q, -(dim-1)); % contraction of Q with the appropriate dimension of A
    A = squeeze(sum(QQ.*A, dim+1)); % permute has shifted "dim" one to the right
end

B = A;

end % function 
