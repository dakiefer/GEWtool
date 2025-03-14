function t = voigt2tensor(M)
% VOIGT2TENSOR returns the tensor t corresponding to the Voigt notated matrix M.
% This function can handle 4th and 3rd order tensors, which correspond to 6x6 or 3x6
% Voigt-notated matrices, respectively.
% 
% See also: D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 
% (Elastic waves in solids), vol. 1, 2 vols. London: ISTE éditions, 2021.
%
% 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    ss = size(M); % size of matrix 
    if all(ss == [6 6])
        t = voigt2tensor4(M);
    elseif all(ss == [3 6])
        t = voigt2tensor3(M);
    else
        error('GEWTOOL:voigt2tensor:wrongarg','Voigt matrix must be of size 6x6 or 3x6');
    end

end % function


% % convert 4th-order tensors
function c = voigt2tensor4(C)
    % voigt2tensor4 - convert 4th-order tensors
    %
    % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

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


function e = voigt2tensor3(E)
    % voigt2tensor3 - convert 3rd-order tensors
    %
    % 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

    validateattributes(E,{'numeric'},{'size',[3 6]});
    e = zeros(3, 3, 3);
    jk = [ {1},{1}; {2},{2}; {3},{3}; {2},{3}; {3},{1}; {1},{2}]; % voigt -> tensor indices
    kj = [ {1},{1}; {2},{2}; {3},{3}; {3},{2}; {1},{3}; {2},{1}]; % voigt -> tensor indices
    % assign e_i{jk} from E_i{J}, where J in {1..6} -> {jk} in {11,22,33,23,32,13,31,12,21} :
    for i = 1:3
        for J = 1:6
            e(i, jk{J,:}) = E(i, J);
            e(i, kj{J,:}) = E(i, J);
        end
    end
end % function