function [s, ea, eb] = slownessCurve(obj, alpha, erot)

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

% % compute basis vectors ea and eb that span the plane orthogonal to erot:
% % (the below code works, but ea and eb are arbitrarily oriented in this plane)
% ebasis = null(erot.'); % svd to compute orthogonal vectors to erot
% ea = ebasis(:,1); % first basis vector spanning the plane normal to erot
% eb = ebasis(:,2); % second basis vector spanning the plane normal to erot

% Compute basis vectors ea and eb that span the plane orthogonal to erot:
% Choose eb such that it is orthogonal to erot and ex; and ea orthogonal to 
% erot and ea. This will lead to ea closely aligned with ex (or ey if rotation axis almost coincides with ex):
ex = [1;0;0];
ey = [0;1;0];
if abs(erot'*ex) < 0.8 % erot is not collinear to ex (let's turn ea towards ex): 
    nb = ex; % shall be normal to direction eb
else % erot is (almost) aligned with ex (let's turn ea towards ey): 
    nb = ey; % shall be normal to direction eb
end
eb = cross(erot,nb); eb = eb/norm(eb);
ea = cross(eb,erot); % is colinear to nb

N = length(alpha);
cs = zeros(3, N);
for i = 1:N
    % Rodrigues' rotation formula (last term is zero due to ea*erot = 0): 
    eki = cos(alpha(i))*ea + sin(alpha(i))*eb; %  + (1 - cos(p0))*sum(ea.*erot)*erot;
    cs(:,i) = obj.wavespeeds(eki);
end
s = 1./cs.'; % slowness

end