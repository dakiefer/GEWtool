function plotSlownessCurve(obj, erot, varargin)
    % plotSlownessCurve - plot the slowness curve around axis erot. 
    % Iteratively calls "wavespeeds()" for different propagation directions and
    % plots the result into a polar plot.
    %
    % Usage:
    % plotSlownessCurve(mat) % rotate around z-axis
    % plotSlownessCurve(mat, erot) % rotate around erot [3x1]
    % plotSlownessCurve(mat, [], 'k.') % additionally pass plotting arguments
    % 
    % See also: wavespeeds, Material.
    %
    % 2022-2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France
    
    if nargin >= 2, erot = erot(:)/norm(erot); end % normalize
    if nargin < 2 | isempty(erot)
        erot = [0;0;1]; % default is rotation in x-y-plane
    end

    alpha = linspace(0, 2*pi, 300); % default angles
    [s, e0] = obj.slownessCurve(alpha, erot);
    polarplot(alpha, s, varargin{:});
    legend({'sl', 'st1', 'st2'})
    title(sprintf('%s: slowness curves around rotation axis [%g, %g, %g]\nangle 0: [%g, %g, %g]',...
        obj.name, erot(1), erot(2), erot(3), e0(1), e0(2), e0(3)))
end
