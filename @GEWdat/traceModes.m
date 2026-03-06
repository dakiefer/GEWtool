function datTraced = traceModes(dat)
% traceModes - reorder modes such that the eigenvalues change smoothly.
% 
% This will enable you to plot the complex wavenumber spectrum as curves for
% each mode, e.g., 
% >> plot3(real(datTraced.k.'), imag(datTraced.k.'), datTraced.w.', '-'); 
% 
% If you find unexpected results, increase nModes and N. 
% 
% 2026 - Daniel A. Kiefer, Institut Langevin, CNRS, ESPCI Paris, France

[kk, ind] = reorderByProximity(dat.k); % matches the modes such that the wavenumbers change as little as possible with frequency
ww = dat.w;
for n = 1:size(ww,2), ww(:,n) = ww(ind(:,n),n); end % same ordering as wavenumbers
if ~isempty(dat.Psi)
    Psi = dat.Psi;
    for n = 1:size(ww,2), Psi(:,n,:) = Psi(ind(:,n),n,:); end % same ordering 
end
if isa(dat,"GEWdatLeaky")
    datTraced = GEWdatLeaky(dat.gew,kk,ww,Psi); % note: dat.beta is extracted from Psi -> we are all good
else
    datTraced = GEWdat(dat.gew,kk,ww,Psi);
end

end