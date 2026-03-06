function [kk_sorted, idx] = reorderByProximity(kk, wm, wp)
% reorderByProximity - reorder modes such that the eigenvalues change smoothly.
% 
% For each column of kk, the entries are re-ordered such that kk(n,i) and
% kk(n,i+1) are the closest possible in magnitude and phase. Usually, kk(n,i) is
% the n-th mode's complex wavenumber at frequency wi.
% 
% The function is based on the Hungarian (assignment) algorithm implemented in
% Matlab's matchpairs().
% 
% Arguments: 
% - kk:    matrix to be sorted/traced
% - wm:    weight for magnitude-matching (default: 1)
% - wp:    weight for phase-matching (default: 1)
%
% Return values: 
% - kk_sorted:  matrix kk but with each column permuted (possibly differently)
% - idx:        indices such that for i = 1:size(kk,2), kk_sorted(:,i) = kk(idx(:,i),i); end
% 
% This function was entirely written by an AI-based code generation system, 
% prompted and adapted by D. A. Kiefer.
% 
% 2026 - Daniel A. Kiefer, Institut Langevin, CNRS, ESPCI Paris, France

if nargin < 2, wm = 1; end % default weight for magnitude-matching
if nargin < 3, wp = 1; end % default weight for phase-matching

[nRows, nCols] = size(kk);

kk_sorted = kk;
idx = zeros(nRows, nCols);

idx(:,1) = 1:nRows;

BIG = 1e6;   % penalty cost

for i = 2:nCols
    prev = kk_sorted(:,i-1);
    curr = kk(:,i);

    valid_prev = isfinite(prev);
    valid_curr = isfinite(curr);

    % Broadcast validity mask
    valid = valid_prev & valid_curr.';

    % --- magnitude term
    mag_prev = abs(prev);
    mag_curr = abs(curr).';

    scale = max(mean(mag_prev(valid_prev)), eps);

    Dmag = abs(mag_curr - mag_prev) / scale;

    % --- phase term (safe ratio)
    ratio = curr.' ./ prev;
    Dphi = abs(angle(ratio)) / pi;

    % --- combined cost
    C = wm * Dmag + wp * Dphi;

    % Penalize invalid comparisons
    C(~valid) = BIG;

    % assignment
    pairs = matchpairs(C, BIG*10); % BIG*10 -> do not allow unmatched -> leads to array size error

    perm = zeros(nRows,1);
    perm(pairs(:,1)) = pairs(:,2);

    kk_sorted(:,i) = curr(perm);
    idx(:,i) = perm;
end

end
