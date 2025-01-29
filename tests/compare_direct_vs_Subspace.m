% We use a FEM discretization with three-noded elements:
h = 1e-3; 
mat = Material('steel'); 
nLay = 5:5:55; 
N = 3; 
nModes = 10;
k = linspace(1e-2, 25, 100)/h; 
plate = Plate(repmat(mat,1,nLay(1)),h/nLay(1),N); 
gew = plate.Lamb;
% N = 10:10:80; %% will be looped over (define after generating the above plate + gew)

% % common options for both cases: 
standardEVP = true;  % default: true
eigenvecs =   true;  % default: true

%% plot to see it 
dat = computeW(gew, k, nModes); 
figure; plot(dat);

%% time the solvers for standardEVP = true (the default) 
% setup options
clear optQRM; 
optQRM.subspace = false; % use QR method instead of subspace method
optQRM.sparse = false; % QR Method does not support sparse matrices
optQRM.standardEVP = standardEVP; 
optQRM.eigenvecs = eigenvecs; 
clear optSub; 
optSub.subspace = true;  
optSub.sparse = true; 
optSub.standardEVP = standardEVP; 
optSub.eigenvecs = eigenvecs; 

% print the parsed options once at the beginning (if show = true): 
disp('The options for the direct method are:'); optQRM.show = true;
datQRM = computeW(gew, k, nModes, optQRM); 
disp('The options for the subspace method are:'); optSub.show = true;
datSub = computeW(gew, k, nModes, optSub);  

% compute timing:
optQRM.show = false; optSub.show = false; % disable display for the loop
tQRM = zeros(size(nLay));
tSub = zeros(size(nLay));
n = zeros(size(nLay));
for i = 1:length(nLay)
    mats = repmat(mat,1,nLay(i));
    plate = Plate(mats,h/nLay(i),N); 
    gew = plate.Lamb; 
    n(i) = size(gew.op.M,1);
    computeQRM = @() computeW(gew, k, nModes, optQRM);
    tQRM(i) = timeit(computeQRM);
    computeSub = @() computeW(gew, k, nModes, optSub);
    tSub(i) = timeit(computeSub);
end

%% plot timing 
figure;
bar(n, [tQRM(:), tSub(:)]);
legend({'direct', 'subspace'}, 'Location','northwest')
xlabel('matrix size $n$'); 
ylabel('computing time in s');
title('Standard EVP')
