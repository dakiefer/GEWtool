% Run using: runtests()
% Test the decoupling of Lamb (longitudinal) and SH (torsional)
% polarizations in a plate (cylinder).
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France

r = [(21.1-3.56)*1e-3, 21.1e-3];    % inner and outer radii
h = r(end) - r(1);                  % thickness
N = 12;                             % number of discretization points
k = linspace(1e-1, 10, 50)/h;       % wavenumbers to solve for
nModes = 10;                        % number of modes to consider
mat = Material('aluminum');

%% test decoupling in plate
plate = Plate(mat,h,N);
gew = plate.fullyCoupled;
datFull = computeW(gew,k,nModes);
gew = plate.Lamb;
datLamb = computeW(gew,k,nModes);
gew = plate.sh;
datSH = computeW(gew,k,nModes);

wwLambSH = sort([datLamb.w, datSH.w],2); 
wwLambSH = wwLambSH(:,1:nModes); % crop 
wwCoupled = datFull.w;
errRel = (wwCoupled - wwLambSH)./mean(wwCoupled(:));

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    figure; hold on; xlim([0,3000]); ylim([0, 6e3]/h);
    title('plate: dispersion curves')
    plot(datFull.k,datFull.w/2/pi,'Color',[1 1 1]*0.5);
    plot(datLamb.k,datLamb.w/2/pi,'.b');
    plot(datSH.k, datSH.w/2/pi, '.r');

    figure; 
    title('relative error in frequency')
    plot(abs(errRel), 'k.');
end

assert(max(abs(errRel(:))) < 1e-6,'The maximum deviation between the fully coupled solution and Lamb+SH solutions is too big.');

%% test decoupling in cylinder
% Note: longitudinal and torsional modes only decouple for n = 0.
n = 0; % only zeroth circumferential order decouple!
cyl = Cylinder(mat,r,N);
gew = cyl.fullyCoupled(n);
datFull = computeW(gew,k,nModes);
gew = cyl.Lamb(n);
datLamb = computeW(gew,k,nModes);
gew = cyl.sh(n);
datSH = computeW(gew,k,nModes);

wwLambSH = sort([datLamb.w, datSH.w],2); 
wwLambSH = wwLambSH(:,1:nModes); % crop 
wwCoupled = datFull.w;
errRel = (wwCoupled - wwLambSH)./mean(wwCoupled(:));

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    figure; hold on; xlim([0,3000]); ylim([0, 6e3]/h);
    title('cylinder: dispersion curves')
    plot(datFull.k,datFull.w/2/pi,'Color',[1 1 1]*0.5);
    plot(datLamb.k,datLamb.w/2/pi,'.b');
    plot(datSH.k, datSH.w/2/pi, '.r');

    figure; 
    title('relative error in frequency')
    plot(abs(errRel), 'k.');
end

assert(max(abs(errRel(:))) < 1e-6,'The maximum deviation between the fully coupled solution and Lamb+SH solutions is too big.');
