function [res] = rayleighLambA(mat, wh, kh)
%RAYLEIGHLAMBS Compute the residuum to the Rayleigh-Lamb frequency equation for 
% antisymmetric modes.
% [1] D. Royer and T. Valier-Brasier, Ondes élastiques dans les solides 1: propagation 
% (Elastic waves in solids 1: propagation), vol. 1, 2 vols. London: ISTE éditions, 2021.


cl = mat.cl; ct = mat.ct;
khl = wh/cl; kht = wh/ct;
khly = sqrt(khl.^2 - kh.^2);
khty = sqrt(kht.^2 - kh.^2);

res = tan(khty/2)./tan(khly/2) + (khty.^2 - kh.^2).^2./(4*kh.^2.*khty.*khly);

end
