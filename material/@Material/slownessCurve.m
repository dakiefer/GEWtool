function [s, ea] = slownessCurve(obj, alpha, erot)

if nargin == 3
    erot = erot(:)/norm(erot); % normalize
elseif nargin < 3
    erot = [0;0;1]; % default is rotation in x-y-plane
elseif nargin < 2
    error('GEWTOOL:slownessCurve:wrongNumberOfArgs', 'Specify the rotation angles as second argument (vector).');
end
if ~isvector(alpha)
    error('GEWTOOL:slownessCurve:wrongArg', 'The set of angles alpha should be specified as a vector.');
end

ebasis = null(erot.'); % svd to compute orthogonal vectors to erot
ea = ebasis(:,1); % first basis vector spanning the plane normal to erot
eb = ebasis(:,2); % second basis vector spanning the plane normal to erot

N = length(alpha);
cs = zeros(3, N);
for i = 1:N
    % Rodrigues' rotation formula (last term is zero due to ea*erot = 0): 
    eki = cos(alpha(i))*ea + sin(alpha(i))*eb; %  + (1 - cos(p0))*sum(ea.*erot)*erot;
    cs(:,i) = obj.wavespeeds(eki);
end
s = 1./cs.'; % slowness

end