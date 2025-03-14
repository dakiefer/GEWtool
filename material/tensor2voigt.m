function M = tensor2voigt(t)
% TENSOR2VOIGT - returns the 6x6 or 3x6 Voigt notated matrix corresponding to the 
% 4th or 3rd order tensor t.
% 
% See also: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    ss = size(t); % size of matrix 
    if ndims(t) == 4 & all(ss == [3 3 3 3])
        M = tensor2voigt4(t);
    elseif ndims(t) == 3 & all(ss == [3 3 3])
        M = tensor2voigt3(t);
    else
        error('GEWTOOL:tensor2voigt:wrongarg','Only 4th and 3rd order tensors are supported. Provide an array of size 3x3x3x3 or 3x3x3.');
    end

end % function 


function C = tensor2voigt4(c)
    % tensor2voigt4 - convert 4th-order tensors
    %
    % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

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


function E = tensor2voigt3(e)
    % tensor2voigt4 - convert 4th-order tensors
    %
    % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    E = zeros(3, 6);
    jk = [ {1},{1}; {2},{2}; {3},{3}; {2},{3}; {3},{1}; {1},{2}]; % voigt -> tensor indices
    % assign E_i{J} from e_i{jk}, where {jk} in {11,22,33,23,32,13,31,12,21} -> J in {1..6} :
    for i = 1:3
        for J = 1:6
            E(i,J) = e(i,jk{J,:});
        end
    end
end % function
