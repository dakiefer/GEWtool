function C = tensor2voigt(c)
% TENSOR2VOIGT - returns the 6x6 Voigt notated stiffness corresponding to the 
% 4th order stiffness tensor c.
% 
% see also: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% % (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022 - Daniel A. Kiefer
% Institut Langevin, Paris, France
% 

validateattributes(c,{'numeric'},{'size',[3 3 3 3]});
C = zeros(6, 6);
ij = [{1},{1}; {2},{2}; {3},{3}; {2},{3}; {3},{1}; {1},{2}]; % voigt -> tensor indices
for I = 1:6
    for J = I:6 % start from I (symmetry in C) 
        C(I, J) = c(ij{I,:}, ij{J, :});
        C(J, I) = c(ij{I,:}, ij{J, :});
    end
end

end % function 
