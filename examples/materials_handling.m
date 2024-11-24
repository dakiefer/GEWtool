%% Represent and load anisotropic material data
%  
% Bibliography:
% [1] K.-J. Langenberg, R. Marklein, and K. Mayer, Ultrasonic Nondestructive 
% Testing of Materials: Theoretical Foundations (translated from german), 1st ed. 
% Boca Raton: CRC Press, 2012. doi: 10.1201/b11724.
%
% 2022-2023 - Daniel A. Kiefer, Institut Langevin, ESPCI Paris, France

%% The "Material" class can be used to represent and manage material data. 
% see also: 
% > help Material
% > help MaterialIsotropic
% The simplest way to create an material object is by loading data from a ".json"
% file by passing the file name: 
mat = Material('silicon'); % load database/silicon.json

% Look at the properties:
mat             % list all properties
CVoigt = mat.C  % save the Voigt stiffness matrix into "CVoigt"

% The elastic data is always stored in the 4th order stiffness tensor "c". This
% allows to represent nonlinear (linearized) cases that break symmetries in the tensor. 
% To create a new material simply create a new json file with similar syntax as
% the existing ones. Note that the data in this file can be provided either in isotropic (Lamé 
% parameters) or anisotropic form (Voigt matrix C). For anisotropic materials
% you should take into account that ex is the wave propagation direction, ey the
% cross-sectional coordinate (of the plate) and ez the out-of plane direction (SH polarization). 
% If you need to switch the ex and ez base vectors, you can use the function 
% permute(). Its default transformation is ex-ey-ez -> ez-ex-ey, 
% i.e., 3->1, 1->2, 2->3, 6->4, 4->5, 5->6, given by perm = [3,1,2,6,4,5]:
mat31 = mat.permute();

% The material reference system can be rotated (see 'help Material.rotateEuler'): 
mat15 = mat.rotateEuler(15/180*pi, 'z'); % rotates around ez (the plate normal)

%% The wave speeds (cl, ct, ct2) are calculated for the general anisotropic case
% by solving the Kelvin-Christoffel equation [1] for phase propagation in direction 
% ek = ex = [1;0;0]. The function wavespeeds() can be used to compute
% propagation velocities in an arbitrary direction:
wavespeeds = mat.wavespeeds([1;1;1]) % the direction vector will be normalized

% The slowness curves around the axis "ez" can be plotted: 
plot(mat) % calls mat.plotSlownessCurve()

% To plot the slowness curves around the axis "erot" = [1,1,1] use:
plot(mat, [1,1,1])
