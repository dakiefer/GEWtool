function [s, e0] = slownessCurve(obj, alpha, erot)

if nargin == 3
    erot = erot(:)/norm(erot); % normalize
elseif nargin < 3
    erot = [0;0;1]; % default is rotation in x-y-plane
elseif nargin < 2
    error('GEWTOOL:slownessCurve:wrongNumberOfArgs', 'Specify the rotation angles as second argument (vector).');
end
if ~isvector(alpha)
    error('GEWTOOL:slownessCurve:wrongArg', 'The angles alpha should be specified as a vector.');
end

eks = null(erot.'); % svd to compute orthogonal vectors to erot
e0 = eks(:,1); % use one of the vectors

N = length(alpha);
cs = zeros(3, N);
for ii = 1:N
    alpha0 = alpha(ii);
    % Rodrigues' rotation formula (last term is zero due to e0*erot = 0): 
    enew = cos(alpha0)*e0 + sin(alpha0)*cross(erot, e0); %  + (1 - cos(p0))*sum(e0.*erot)*erot;
    cs(:,ii) = obj.wavespeeds(enew);
end
s = 1./cs.'; % slowness

end