function T = stress(obj, F, np, udof)
% stress - Compute stress for the displacement gradient F = grad(u).
% The last two dimensions of F are interpreted as the dimensions of the 2nd rank 
% displacement gradient F. Only the degrees of freedom "udof" are considered.
% 
% Literature: 
% K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive
% Testing of Materials: Theoretical Foundations (translated from German),
% 1st ed. Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
%
% Arguments:
% - obj:   Material object 
% - F:     displacement gradient F = grad(u): [Nudof x Nudof] with possibly preceeding dimensions
% - np:    (struct) normalization parameters, i.e., SI value to interpret as unit values
% - udof:  [Nudof x 1] indicating which components of stiffness tensor to use
%
% Return values:
% - T:     (same size as F) stress tensors
% 
% 2026 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, CNRS, France

if nargin < 4, udof = [1;2;3]; end

dims = 1:ndims(F)+2; % new dimensions
dims = [dims(1:end-4) dims(end-1:end) flip(dims(end-3:end-2))]; % e.g. [1 2 3 6 7 5 4], where 6 and 7 are new dimensions
%        ↑ leading     ↑ new          ↑ transpose of F-dimensions

c = shiftdim( obj.c(udof,udof,udof,udof), -(length(dims)-4))/np.c0; % norm. stiffness: shift to correct dimension
F = permute(F, dims); % arrange so that dimensions of F and c overlap that are to be contracted
last   = length(dims); 
before = length(dims)-1;
T = sum(sum( c.*F , last), before); % doulbe contraction of last two dimensions

end
