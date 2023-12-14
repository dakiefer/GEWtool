% Run using: runtests()
% Test computing a target region of the dispersion curves.
%
% see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, France

mat = Material('aluminum');
h = 1e-3; 
N = 25; 
plate = Plate(mat,h,N);
gew = plate.Lamb;
k = linspace(1e-2, 20, 100)*1e3; % use sufficient points to compare k- and w-solutions
datFull = computeW(gew,k,20);

%% target in frequency
wtarget = 2*pi*10e6;
opts.subspace = true; opts.target = wtarget;
dat = computeW(gew,k,3,opts);

% look at distance between the two solutions:
pointsFull = [datFull.k(:)/1e3, datFull.w(:)/2/pi/1e6];
pointsRegion = [dat.k(:)/1e3, dat.w(:)/2/pi/1e6];
distances = inf(size(pointsRegion,1),1);
for i = 1:size(pointsRegion,1)
    diff = pointsFull - pointsRegion(i,:);
    distsq = diff(:,1).^2 + diff(:,2).^2;
    distances(i) = sqrt(min(distsq));
end
distMax = max(distances);
assert(distMax < 1e-4);

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    figure; hold on; 
    plot(datFull.k/1e3, datFull.w/2/pi/1e6, '-', 'Color',[1 1 1]*0.5);
    plot(dat.k/1e3, dat.w/2/pi/1e6, '.b', 'MarkerSize',8);
    yline(wtarget/2/pi/1e6,'r')
    % display max deviation: 
    distMax
end

%% target in wavenumber
w = 2*pi*linspace(1e-2,20,300)*1e6;
ktarget = 10e3;
opts.subspace = true; opts.target = ktarget;
dat = computeK(gew,w,2,opts);

% look at distance between the two solutions:
pointsFull = [datFull.k(:)/1e3, datFull.w(:)/2/pi/1e6];
pointsRegion = [dat.k(:)/1e3, dat.w(:)/2/pi/1e6];
indPos = real(pointsRegion(:,1)) > 1e-2;
pointsRegion = pointsRegion(indPos,:);
indReal = abs(imag(pointsRegion(:,1))) < 0.01*abs(real(pointsRegion(:,1)));
pointsRegion = pointsRegion(indReal,:);
distances = inf(size(pointsRegion,1),1);
for i = 1:size(pointsRegion,1)
    diff = pointsFull - pointsRegion(i,:);
    distsq = diff(:,1).^2 + diff(:,2).^2;
    distances(i) = sqrt(min(distsq));
end
distMax = max(abs(distances));
assert(distMax < (k(2)-k(1))/1e3);

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    figure; hold on; 
    plot(datFull.k/1e3, datFull.w/2/pi/1e6, '-', 'Color',[1 1 1]*0.5);
    plot(real(dat.k)/1e3, dat.w/2/pi/1e6, '.b', 'MarkerSize',8);
    xline(ktarget/1e3,'r')
    % display max deviation: 
    distMax
end

%% target in wavenumber (linearized)
w = 2*pi*linspace(1e-2,20,300)*1e6;
ktarget = 10e3;
opts.subspace = true; opts.target = ktarget;
gewLin = copy(gew);
gewLin.linearizeInK;
dat = computeK(gewLin,w,2,opts);

% look at distance between the two solutions:
pointsFull = [datFull.k(:)/1e3, datFull.w(:)/2/pi/1e6];
pointsRegion = [dat.k(:)/1e3, dat.w(:)/2/pi/1e6];
indPos = real(pointsRegion(:,1)) > 1e-2;
pointsRegion = pointsRegion(indPos,:);
indReal = abs(imag(pointsRegion(:,1))) < 0.01*abs(real(pointsRegion(:,1)));
pointsRegion = pointsRegion(indReal,:);
distances = inf(size(pointsRegion,1),1);
for i = 1:size(pointsRegion,1)
    diff = pointsFull - pointsRegion(i,:);
    distsq = diff(:,1).^2 + diff(:,2).^2;
    distances(i) = sqrt(min(distsq));
end
distMax = max(abs(distances));
assert(distMax < (k(2)-k(1))/1e3);

% plot if the variable "show" has been set to true
if exist('show', 'var') && show
    figure; hold on; 
    plot(datFull.k/1e3, datFull.w/2/pi/1e6, '-', 'Color',[1 1 1]*0.5);
    plot(real(dat.k)/1e3, dat.w/2/pi/1e6, '.b', 'MarkerSize',8);
    xline(ktarget/1e3,'r')
    % display max deviation: 
    distMax
end
