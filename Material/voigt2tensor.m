function c = voigt2tensor(C)
% VOIGT2TENSOR returns the 4th order stiffness tensor c corresponding to the Voigt
% notated matrix C.
% 
% see also: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

validateattributes(C,{'numeric'},{'size',[6 6]});
if ~issymmetric(C), error('Voigt matrix C must be symmetric.'); end
c = zeros(3, 3, 3, 3);
ij = [ {1},{1}; {2},{2}; {3},{3}; {2},{3}; {3},{1}; {1},{2}]; % voigt -> tensor indices
ji = [ {1},{1}; {2},{2}; {3},{3}; {3},{2}; {1},{3}; {2},{1}]; % voigt -> tensor indices
for I = 1:6
    for J = I:6 % start from I (symmetry in C) -> explicitly assign major symmetry below
        % minor symmetries: 
        c(ij{I,:}, ij{J,:}) = C(I, J); % cijkl
        c(ij{I,:}, ji{J,:}) = C(I, J); % cijlk
        c(ji{I,:}, ij{J,:}) = C(I, J); % cjikl
        c(ji{I,:}, ji{J,:}) = C(I, J); % cjilk
        % major symmetry:
        c(ij{J,:}, ij{I,:}) = C(I, J); % cklij
        c(ij{J,:}, ji{I,:}) = C(I, J); % cklji
        c(ji{J,:}, ij{I,:}) = C(I, J); % clkij
        c(ji{J,:}, ji{I,:}) = C(I, J); % clkji
    end
end

end % function
