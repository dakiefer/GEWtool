
mat = MaterialIsotropic('steel');   % load from database (or create your own)
h = 1e-3;                        % thickness in m
k = linspace(1e-2, 12, 20)/h;   % wavenumbers to solve for
nModes = 4; 

N1 = 7;
plate1 = Plate([mat mat mat], h/3, N1); 
gew1 = plate1.Lamb; 

N2 = 8;
plate2 = Plate([mat mat], h/2, N2);
gew2 = plate2.Lamb;

gew = [gew1 gew2]; 
dat = computeW(gew, k, nModes);
w0 = 2*pi*gew1.np.fh0/gew1.np.h0; 

i = 1; % gew(i) will be used for scalar tests 

if exist('show','var') && show 
    figure(1); clf; hold on
    plot(dat(1).k/1e3, dat(1).w/2/pi/1e6,'*','SeriesIndex',1,'DisplayName','bi-layer');
    plot(dat(2).k/1e3, dat(2).w/2/pi/1e6,'r.','SeriesIndex',2,'DisplayName','tri-layer');
    ylim([0, dat(1).w(end,1)/2/pi/1e6*1.1]); 
    xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
    legend(legendUnq, 'Location', 'southeast')
end

% define list of functions and properties: 
% each row is a test consisting of 
% %         function Handle,     integral qty?,  number of field components, args }
funcList = {
            @energyDensityElastic,      false,      1,      {}; 
            @energyDensityKinetic,      false,      1,      {}; 
            @energyElastic,             true,       1,      {}; 
            @energyKinetic,             true,       1,      {}; 
            @energyTotal,               true,       1,      {}; 
            @energyVel,                 true,       1,      {}; 
            @energyVelAxial,            true,       1,      {}; 
            @energyVelMag,              true,       1,      {}; 
            @energyVelTransverse,       true,       1,      {}; 
            @energyVelVec,              true,       2,      {}; 
            @excitabilityLUS,           true,       1,      {'top'};
            @groupVel,                  true,       1,      {};
            @groupVelAxial,             true,       1,      {}; 
            @powerFlux,                 true,       2,      {};
            @powerFluxAxial,            true,       1,      {};
            @powerFluxMag,              true,       1,      {};
            @powerFluxTransverse,       true,       1,      {};
            @poyntingVec,               false,      2,      {};
            @poyntingVecAxial,          false,      1,      {};
            @poyntingVecTransverse,     false,      1,      {};
            @strain,                    false,    [2 2],    {};
            @stress,                    false,    [2 2],    {};
            @fieldComponents,           false,      2,      {}; 
            @velocity,                  false,      2,      {};
            };  

% extra wurscht für: 
warning('Implement tests for crossPowerFlux(), normalizeComplex(), normalizeReal()!')

%% two problems - same result
err = abs(dat(1).w - dat(2).w)/w0; 
if exist('show','var') && show  
    errMax = max(err(:))
end
assert(max(err(:)) < 1e-6);

%% test list of functions 
pass = false(1,size(funcList,1));
for j = 1:size(funcList,1) 
    pass(j) = testFunction(dat, funcList{j,:});
end
if any(~pass)
    fHandles = funcList(~pass,1);
    str = 'The following functions failed:\n'; 
    for j = 1:length(fHandles)
        str = [str, char(fHandles{j}), '()\n']; 
    end
    fprintf(2,str)
    error('GEWTOOL/GEWdat', str);
else
    disp('All functions computed correctly. Testing some extra functions in the following.')
end

%% eigenVecs()
% Check behavior of scalar input:
result = eigenVecs(gew(i), dat(i).u);
% Check that the size is correct: 
% Two layers with N1 nodes = 2*N1-1 distinct nodes. Lamb waves with 2*(2*N1-1)
% unknowns.
sz = [length(k), nModes, (gew(i).geom.nLay*gew(i).geom.N(i)-(gew(i).geom.nLay-1))*2]; % expected size
validateattributes(result, {'numeric'}, {'3d', 'size', sz})  
% Check multiple Waveguide input: 
result = eigenVecs(gew, {dat.u}); % pass u as cell array for each gew(i) 
validateattributes(result, {'cell'}, {'vector', 'numel', length(gew)})  

%% extractModes()
indk = 3:5; 
indw = 2; 
result = extractModes(dat(i), indk, indw); 
assert(result.Nk == length(indk))
assert(result.Nw == length(indw))
validateattributes(result.Psi, {'numeric'}, {'size', [length(indk), length(indw), size(dat(i).Psi,3)]})
result = extractModes(dat, indk, indw); 
for j = 1:length(gew)
    assert(result(j).Nk == length(indk))
    assert(result(j).Nw == length(indw))
    validateattributes(result(j).Psi, {'numeric'}, {'size', [length(indk), length(indw), size(dat(j).Psi,3)]})
end


%% 
% clear dd; 
% clc
% dd = excitabilityLUS(dat,'top')
% dd{1}

function pass = testFunction(dat, fHandle, integr, nComponents, args)
    fprintf('testing function %s()\n', char(fHandle));
    % test computation on a single GEWdat object:
    i = 1; % test on this one
    pass = true; 
    try
        result = fHandle(dat(i), args{:}); 
        if integr
            validateIntegral(result, dat(i), nComponents); 
        else
            validateContinuous(result, dat(i), nComponents); 
        end
        % test computation on a vector of GEWdat objects:
        result = fHandle(dat, args{:}); 
        if integr
            validateIntegral(result, dat, nComponents); 
        else 
            validateContinuous(result, dat, nComponents); 
        end
    catch
        pass = false; 
    end
end

function validateIntegral(result, dat, ndims)
    i = 1; % only test the i-th waveguide problem dat(i)
    if isscalar(dat)
        validateattributes(result, {'numeric'}, {'size', [dat.Nk, dat.Nw, ndims]})
    else 
        validateattributes(result, {'cell'}, {'vector', 'numel', length(dat)})
        validateattributes(result{i}, {'numeric'}, {'size', [dat(i).Nk, dat(i).Nw, ndims]})
    end
end

function validateContinuous(result, dat, ndims)
    i = 1; % only test the i-th waveguide problem dat(i)
    if isscalar(dat)
        validateattributes(result, {'cell'}, {'vector', 'numel', dat.gew.geom.nLay})
        validateattributes(result{1}, {'numeric'}, {'size', [dat.Nk, dat.Nw, dat.gew.geom.N(1), ndims]})
    else 
        validateattributes(result, {'cell'}, {'vector', 'numel', length(dat)})
        validateattributes(result{i}, {'cell'}, {'vector', 'numel', dat(i).gew.geom.nLay})
        validateattributes(result{i}{1}, {'numeric'}, {'size', [dat(i).Nk, dat(i).Nw, dat(i).gew.geom.N(1), ndims]})
    end
end