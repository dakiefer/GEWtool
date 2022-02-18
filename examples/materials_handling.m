%% Represent and load anisotropic material data
% 
% 2022 Daniel Kiefer
% Institut Langevin, Paris, France
% 

clear all
%% The "Material" class can be used to represent and manage material data. The
% simplest way to create an material object is by loading data from a ".json"
% file by passing the file name: 
mat = Material('zircaloy'); % load database/zircaloy.json

% Look at the properties:
mat         % list all properties
C = mat.C   % save the Voigt stiffness matrix into "C"

% The elastic data is always stored in the property "C", the Voigt notated
% stiffness matrix. The 4th order stiffness tensor "c" is created from there on
% request. 
% To create a new material simply create a new json file with similar syntax as
% the existing ones. Note that the data can be provided either in isotropic (Lamé 
% parameters) or anisotropic form (Voigt matrix C). For anisotropic materials
% you should take into account that ex is the wave propagation direction, ey the
% cross-sectional coordinate (of the plate) and ez the out-of plane direction. 
% If you need to switch the ex and ez base vectors, you can use the function 
% permute13():
mat31 = mat.permute13();
C31 = mat31.C

%% The wave speeds (cl, ct, ct2) are calculated for the general anisotropic case
% by solving the Kelvin-Christoffel equation [1] for phase propagation in direction 
% ek = ex = [1;0;0]. The function wavespeeds() can be used to compute
% propagation velocities in an arbitrary direction:
vs = mat.wavespeeds([1;1;1]) % the direction vector will be normalized

% The slowness curves around the axis "ez" can be plotted: 
plot(mat) % calls mat.plotSlownessCurve()

% To plot the slowness curves around the axis "erot" = [1,1,1] use:
plot(mat, [1,1,1])

%% 
% Bibliography:
% [1] K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive 
% Testing of Materials: Theoretical Foundations (translated from german), 1st ed. 
% Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
%
