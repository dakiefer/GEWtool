function plotSlownessCurve(varargin)
    % plotSlownessCurve - plot the slowness curve around axis erot. 
    % Iteratively calls "wavespeeds()" for different propagation directions and
    % plots the result into a polar plot.
    %
    % Usage:
    % plotSlownessCurve(mat) % rotate around z-axis
    % plotSlownessCurve(mat, erot) % rotate around erot [3x1]
    % 
    % See also: wavespeeds, Material.
    %
    % 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
    
    if nargin==1
        obj = varargin{1}; erot = [0;0;1];
    elseif nargin>=2
        obj = varargin{1};
        erot = varargin{2}; erot = erot(:)/norm(erot);
    end
    if nargin < 2,  end
    eks = null(erot.'); % svd to compute orthogonal vectors to erot
    e0 = eks(:,1); % use one of the vectors

    N = 300; % number of points
    p = linspace(0, 2*pi, N);
    vs = zeros(3, N);
    for ii = 1:N
        p0 = p(ii);
        % Rodrigues' rotation formula (last term is zero due to e0.erot = 0): 
        enew = cos(p0)*e0 + sin(p0)*cross(erot, e0); %  + (1 - cos(p0))*sum(e0.*erot)*erot;
        vs(:,ii) = obj.wavespeeds(enew);
    end
    s = 1./vs.'; % slowness

    polarplot(p, s)
    legend({'sl', 'st1', 'st2'})
    title(sprintf('%s: slowness curves around rotation axis [%g, %g, %g]\nangle 0: [%g, %g, %g]',...
        obj.name, erot(1), erot(2), erot(3), e0(1), e0(2), e0(3)))
end
