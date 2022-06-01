function [u] = rayleighLambU(mat, wh, kh, y, sym)
%RAYLEIGHLAMBU Computes the displacement corresponding to a solution of the
%Rayleigh-Lamb equation.
% - mat: material structure containing mat.c and mat.rho
% - wh: angular frequency-thickness product
% - kh: wavenumber-thickness product
% - y: thickness coordinate(s) normalized to [-0.5, 0.5];
% - sym: select between symmetric (S) and anti-symmetric (A) modes
% 
% [1] D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 1: propagation 
% (Elastic waves in solids 1: propagation), vol. 1, 2 vols. London: ISTE éditions, 2021.

cl = sqrt((mat.c(1,1,2,2) + 2*mat.c(2,3,2,3))/mat.rho); % better accuracy than mat.cl
ct = sqrt(mat.c(2,3,2,3)/mat.rho);
khl = wh/cl; kht = wh/ct;
khly = sqrt(khl.^2 - kh.^2);
khty = sqrt(kht.^2 - kh.^2);
y = y(:); 

switch sym
    case {'S', 'sym'}
        a = 0; 
    case {'A', 'anti'}
        a = pi/2;
end

ux = khty*(cos(khty*y + a) - 2*kh^2/(kh.^2-khty^2)*cos(khty/2 + a)/cos(khly/2 + a)*cos(khly*y + a));
uy = -1i*kh*(sin(khty*y + a) + 2*khly*khty/(kh^2 - khty^2)*cos(khty/2 + a)/cos(khly/2 + a)*sin(khly*y + a));

u = [ux, uy];

end
