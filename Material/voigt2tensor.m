function c = voigt2tensor(C)
% VOIGT2TENSOR returns the 4th order stiffness tensor c corresponding to the Voigt
% notated matrix C.
%
% 2022 - Daniel Kiefer
% Institut Langevin, Paris, France

validateattributes(C,{'numeric'},{'size',[6 6]});
if ~issymmetric(C), error('Voigt matrix C must be symmetric.'); end
c = zeros(3, 3, 3, 3);
IJ = [ {1}, {1}; {2}, {2}; {3}, {3}; {2}, {3}; {3}, {1}; {1}, {2}]; % voigt -> tensor indices
for i = 1:6
    for j = 1:6
        c(IJ{i,:}, IJ{j,:}) = C(i, j);
        c(IJ{i,end:-1:1}, IJ{j,:}) = C(i, j);
        c(IJ{i,:}, IJ{j,end:-1:1}) = C(i, j);
        c(IJ{i,end:-1:1}, IJ{j,end:-1:1}) = C(i, j);
    end
end

end
