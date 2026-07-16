%% Compute leaky Lamb waves in a plate radiating into a fluid medium
% This example shows how to compute leaky Lamb waves. 
%
%                                cf, rhof  ↑ fluid medium (infinite)
%  _↑z____________________________________ top surface (continuitiy conditions)
%   →x     --> k  guided wave      c,rho   plate (stiffness c, mass density rho)
%  _______________________________________ bottom surface (traction free)
%                                          ↓ vacuum
%
% 2026 - Daniel A. Kiefer, Institut Langevin, CNRS, France

mat = Material('brass'); 
ext = MaterialFluid('water');
h = 1e-3;  % thickness in m
N = 12;    % number of nodes
w = 2*pi*[linspace(0.01, 3, 500)].'*1e6; % angular frequency in rad/s
plate = Plate({mat ext}, [h inf], N);    % leaky waveguide model
gew = plate.Lamb; tic;                   % choose Lamb-polarization
dat = computeK(gew, w, 4*N); toc;        % solve

%% plot 
figure(1); clf; 
subplot(2,1,1); 
plot(dat);
subplot(2,1,2); 
plot(dat.w/2/pi/1e6, imag(dat.k)/1e3, 'k.'); ylim([0, 0.2]);
xlabel('frequency $\omega/2\pi$ in MHz'); ylabel('attenuation $\Im k$ in rad/mm')
