function p = momentum(obj, v, np, ~)
% momentum - Compute the linear momentum for the velocity v.
% The last dimension of v is interpreted as the dimension of the 1st rank
% velocity tensor. Only the degrees of freedom "udof" are considered.
% 
% Literature: 
% K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive
% Testing of Materials: Theoretical Foundations (translated from German),
% 1st ed. Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
%
% Arguments:
% - obj:   Material object 
% - v:     velocity v = -iw*u [Nudof x 1] with possibly preceeding dimensions
% - np:    (struct) normalization parameters, i.e., SI value to interpret as unit values
% - udof:  [Nudof x 1] indicates which components of displacements are involved
%
% Return values:
% - p:     (same size as v) momentum 
% 
% 2026 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris | PSL, CNRS, France

rho = obj.rho/np.rho0; 
p = rho*v; % field v is already normalized

end
