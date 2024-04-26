% Compare circumferential waves in a tube to Qu et al. [1].
%
% [1] J. Qu, Y. Berthelot, and Z. Li, “Dispersion of Guided Circumferential
% Waves in a Circular Annulus,” in Review of Progress in Quantitative
% Nondestructive Evaluation: Volume 15A, D. O. Thompson and D. E. Chimenti,
% Eds., Boston, MA: Springer US, 1996, pp. 169–176. doi:
% 10.1007/978-1-4613-0383-1_21.
%
% Run using: runtests(). see also: 
% https://fr.mathworks.com/help/matlab/matlab_prog/write-script-based-unit-tests.html
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

rho = 2700; % should not matter
E = 6.4721e+10; 
nu = 0.25; % nu as in Qu et al. [1]
[lbd, mu] = MaterialIsotropic.Enu2lame(E,nu);
mat = MaterialIsotropic('prototype', lbd, mu, rho);
eta = 0.1; a = 1e-3; b = a/eta; h = b-a; % a: inner radius, b: outer radius
N = 10;                                  % number of nodes (dictates accuracy)
k = linspace(1e-3, 10, 70)/h;            % wavenumber list to solve for
cyl = CylinderCircumferential(mat, [a b], N); % waveguide 

% compute circumferential waves (ur-uphi-polarized)
gew = cyl.Lamb;              % assembles matrices 
dat = computeW(gew, k, 8);   % compute 

%% Compare to Zemanek
pic=imread('data/Qu_circumferential.png');   % load reference
fig = figure(1); clf; hold on
image([0,10],[10,0],pic); axis xy
xlabel('$k h$'); ylabel('$\omega h / c_t$')
title('reference: Qu doi 10.1007/978-1-4613-0383-1_21')

% plot GEWtool result on top
phGEW = plot(dat.k*h, dat.w*h/mat.ct, 'r.', 'MarkerSize', 8);
ylim([0, 10]); xlim([0 10]);
legend([phGEW(1)], {'GEWtool'}, 'Location','south east')

% request user to evaluate test:
assert( userTestConfirmation(fig) )
