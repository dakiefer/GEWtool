addpath('../experimental')
mat = Material('aluminum');
h = 1; 
w0 = 2*pi*5000; % for reference and mode selection
k0 = 5.59596;  % for reference and mode selection (also: 3.6180, 3.8674, 5.59596, 7.8648)
N = 5:5:40;
k = linspace(1e-3, 15, 200);

solv = @(guw) computeW(guw,k); % Note: this function shall not call tic/toc!!
impls = {@(N) implCurrent(mat,h,N), @(N) Lamb_matrices_SEM(mat,h,N)};

%% current implementation 
figure(1), hold off
chron = zeros(length(N)+1,length(impls));
implsName = {};
for n = 1:length(impls)
    impl = impls{n};
    func = functions(impl); 
    implsName{n} = func.function;
    fprintf('Testing %s:\n', implsName{n})
    tic; 
    for ii = 1:length(N)
        guw = impl(N(ii));
        dat = solv(guw); 
        chron(ii+1,n) = toc;
        fprintf('dof = %d. elapsed time: %g.\n', size(guw.op.M,1),  chron(ii+1,n)-chron(ii,n));
    end
    chron(ii+1,n) = toc;
    fprintf('Total elapsed time: %g.\n\n', chron(ii+1,n));
    plot(dat.k(:), dat.w(:)/2/pi, '.'); ylim([0, 6e3]/h); hold on;
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz'), drawnow
end

figure(2)
bar(N, diff(chron))
legend(implsName, 'Location','northwest')
xlabel('Discretization order N')
ylabel('Runtime in s')


%% auxilary function
function guw = implCurrent(mat,h,N)
    plate = Plate(mat, h, N);
    guw = plate.Lamb();
end
