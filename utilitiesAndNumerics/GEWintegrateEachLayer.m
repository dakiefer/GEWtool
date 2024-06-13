function [I] = GEWintegrateEachLayer(gew, f, varargin)
% GEWintegrateEachLayer - Integrate f over each layer of the waveguide separately.
% 
% Calls GEWintegrate() iteratively for each layer of the Waveguide object "gew"
% and concatenates the result in the cell array "I" that is returned. I{1} is
% the integrated field "f{1}" of layer 1. 
% 
% See also: GEWintegrate, lglnodes, lgwt
%
% 2024 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

I = cell(1,gew.geom.nLay);
for l = 1:gew.geom.nLay
    I{l} = GEWintegrate(gew, f, l, varargin{:}); 
end

end
