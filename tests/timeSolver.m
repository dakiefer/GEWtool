%% computational time 
% Times computational cost required to solve for dispersion curves with different
% discretization orders and different solvers.
% 
% 2022 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

mat = Material('steel');
h = 1; 
N = 5:5:50;
k = linspace(1e-1, 15, 200);
w = k*mat.ct;

% list of solvers to compare:
solver_list = {@(gew) computeW(gew,k)};

%% current implementation 
figure(1), hold off
chron = zeros(length(N),length(solver_list));
n = zeros(length(N), length(solver_list)); % allocate for matrix size (depends on problem solved)
solverName = {};
for i = 1:length(solver_list)
    solver = solver_list{i};
    func = functions(solver); 
    solverName{i} = func.function; % get the name of the solver function
    fprintf('Testing %s:\n', solverName{i})
    for j = 1:length(N)
        plate = Plate(mat, h, N(j));
        gew = plate.Lamb();
        tic;
        dat = solver(gew); 
        chron(j,i) = toc;
        n(j,i) = size(gew.op.M,1);
%         chron(ii,n) = timeit(@()solver(gew)); % not really necessary as we do lots of for loops in the solver function
        fprintf('dof = %d. solution time: %g.\n', size(gew.op.M,1),  chron(j,i));
    end
    plot(dat.k(:), dat.w(:)/2/pi, '.'); ylim([0, 6e3]/h); hold on;
    xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz'), drawnow
end

figure(1)
legend(solverName, 'Location','northwest')

figure(2), hold on,
bar(n, chron)
legend(solverName, 'Location','northwest')
xlabel('matrix size n')
ylabel('solution time in s')
