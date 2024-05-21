% Test that computeW() finds the same dispersion curves with different options. 
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

N = 12;
mat = MaterialIsotropic('steel');
h = 1e-3;
k = linspace(1e-1, 12, 20)/h;
plate = Plate(mat, h, N);
gew = plate.LambS;
nModes = 4;

datDefault = computeW(gew, k, nModes);

if exist('show', 'var') && show
    figure(1); clf; hold on
    hS = plot(datDefault.k/1e3, datDefault.w/2/pi/1e6, 'k'); % plot symmetric
    ylim([0, 6]); 
    xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
    title(sprintf('Lamb waves in %gmm %s', h/1e-3, mat.name))
end

%% test options
for a = 0:2^4-1
    clear opts; 
    opts.standardEVP = bitget(a,4); 
    opts.eigenvecs = bitget(a,3);
    opts.subspace = bitget(a,2); 
    opts.sparse = bitget(a,1); 
    
    if opts.sparse && ~opts.subspace
        continue; % sparse matrices only in combination with subspace method
    end
    dat = computeW(gew, k, nModes, opts);
    Dw = abs(dat.w - datDefault.w)*h/mat.ct; % difference in dimensionless freq
    assert(all(Dw < 1e-6, 'all'))
    if opts.eigenvecs
        % MATLAB  does not use one consistent normalization of eigenvectors:
        sU = size(dat.u{1}); 
        udat = reshape(dat.u{1}, sU(1), sU(2), sU(3)*sU(4)); % combine ux-uy
        udef = reshape(datDefault.u{1}, sU(1), sU(2), sU(3)*sU(4));
        udat = udat./vecnorm(udat,2,3); % normalize to unit 2-norm
        udef = udef./vecnorm(udef,2,3); % normalize to unit 2-norm
        Du = (abs(udat) - abs(udef)).^2; % compare only magnitude (eigenvecs defined up to complex scalar))
        dev = vecnorm(Du,2,3); % deviation: 2-norm of the differences
        assert(all(dev < 1e-6, 'all'))
    end
    if exist('show', 'var') && show
        disp('passed with options'); disp(opts);
    end
end
