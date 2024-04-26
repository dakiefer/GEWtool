% Compare the complex spectrum of axial waves in a solid rod to the solutions
% by Zemanek [1].
%
% [1] J. Zemanek Jr., “An Experimental and Theoretical Investigation of Elastic
% Wave Propagation in a Cylinder,” The Journal of the Acoustical Society of
% America, vol. 51, no. 1B, pp. 265–283, Jan. 1972, doi: 10.1121/1.1912838.
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France
%        Clemens Grünsteidl, RECENDT, Austria

rho = 2700; %should not matter
E = 6.4721e+10; nu = 0.3317; % same as Zemanek ref.
[lbd, mu] = MaterialIsotropic.Enu2lame(E,nu);
alu = MaterialIsotropic('alu', lbd, mu, rho);
r = 1.00e-3;                          % radius r
N = 22;                               % number of nodes (dictates accuracy)
w = alu.ct/r*linspace(1e-3, 15, 200); % frequency list
cyl = Cylinder(alu, [0, r], N);       % create waveguide description 

% compute longitudinal waves (ux and ur displacements only)
long = cyl.longitudinal;              % assembles matrices 
dat = computeK(long, w, 16);          % compute 16 modes (there are inf many complex)

%% Compare to Zemanek
pic=imread('data/rod_Zemanek.png');   % load reference
fig = figure(1); clf; hold on
image([0,10],[15,0.],pic); axis xy
xlabel('$\gamma r$'); ylabel('$\omega r / c_t $')
title(sprintf('reference: Zemanek doi 10.1121/1.1912838', r/1e-3))

% plot GEWtool result on top
phGEW = plot(real(dat.k)*r, dat.w*r/alu.ct, 'r.', 'MarkerSize', 8);
ylim([0, 15]); xlim([0 10]);
legend([phGEW(1)], {'GEWtool'}, 'Location','south east')

% request user to evaluate test:
assert( userTestConfirmation(fig) )
