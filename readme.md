# GEWtool [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10114243.svg)](https://doi.org/10.5281/zenodo.10114243)


**Compute guided elastic wave (GEW) dispersion in MATLAB.**

`GEWtool` provides extremely fast and reliable computation of guided elastodynamic waves (GEWs) in plates and cylinders. It is simple to use yet provides full access to the computational results as well as the underlying code. You are welcome to contribute to this open-source project.

**Features**:

- Multi-layered plates, tubes and rods
- Finds all solutions, super fast 
- General anisotropy, dissipation
- Compute real frequencies, complex wavenumbers, or ZGV points
- Choose polarization (Lamb/SH/coupled) and symmetry (S/A)

[![GitHub](resourcesAndDeps/img/logo_github.svg)](https://github.com/dakiefer/GEWtool) Code repository: [https://github.com/dakiefer/GEWtool](https://github.com/dakiefer/GEWtool)

## Example: Lamb waves

```matlab
mat = Material('steel');         % load from database (or create your own)
h = 1e-3;                        % thickness in m
N = 12;                          % number of nodes (dictates accuracy)
k = linspace(1e-2, 12, 100)/h;   % wavenumbers to solve for
plate = Plate(mat, h, N);        % create waveguide description 
gews = plate.LambSA; tic;        % choose S+A Lamb waves (assembles matrices)
dat = computeW(gews, k, 4); toc; % solve and save 4 modes (argument optional)
plot(dat(1).k, dat(1).w/2/pi, 'b'); hold on;        % plot symmetric
plot(dat(2).k, dat(2).w/2/pi, 'r'); ylim([0, 6e6]); % plot anti-symmetric
xlabel('wavenumber k in rad/m'), ylabel('frequency f in Hz')
```

output:
`> Elapsed time is 0.010129 seconds.` 

![Lamb waves in steel](resourcesAndDeps/img/dispersion_lamb_steel.png)

Proceed by inspecting the laser-ultrasonic excitability of the waves computed above (product of tangential and normal displacements ux·uy):

```matlab
k = linspace(1e-2, 12, 200)/h;          % use more wavenumbers
gew = plate.Lamb;                       % choose all Lamb waves
dat = computeW(gew, k, 7);              % compute
exc = excitabilityLUS(gew, dat, 'top'); % ux*uy at top surface (value of 1 at 100x median)
exc = 20*log10(exc);                    % in decibel
scatter(dat.k(:)/1e3, dat.w(:)/2/pi/1e6, 15, exc(:), 'filled'), ylim([0, 6]);
colormap(flipud(colormap)); cb = colorbar; caxis([-50, 0]);
xlabel('wavenumber k in rad/mm'), ylabel('frequency f in MHz')
title('laser-ultrasonic excitability in dB')
```

![LUS excitability](resourcesAndDeps/img/lus_excitability.png)

## Installation 

Add `GEWtool` and its subfolders to the Matlab path and save it for future sessions. To achieve this:

1. change to the `GEWtool` folder (e.g., by navigating or using `cd`)
2. execute `install`

Enjoy!

## Getting started

To get started, explore the `examples` directory. 

You can also display help for all functions and classes, e.g., by typing `help Plate`. The most important ones are

- Material representation: `Material`, `MaterialIsotropic`
- Waveguides:  `Plate`, `Cylinder`
- Solvers: `computeW`, `computeK`, `computeZGV`

## Known limitations 

- Cylinders: The postprocessing tools provided in the folder `GEWdat` are only designed for plates. *Do not use them for cylinders*. Dispersion curves of axial waves in cylinders compute correctly, nonetheless. 
- Cylinders: only axial waves are supported for now. 
- Leaky waves: no support for now.

Contact me if you have questions:  [daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr)

## Mathematical and physical background

GEWtool implements the *Spectral Element Method (SEM)* (higher-order Finite Elements) to solve the *waveguide problem*, i.e., the boundary value problem that describes wave propagation in the structure. Such an approach is commonly qualified as 'semi-analytical'. Contrary to classical root-finding of the characteristic equation, this method does not miss solutions. Moreover, unlike Finite Elements, the Spectral Elements lead to small but dense matrices. 

Solusions are computed with machine precision provided you have set the discretization order `N` sufficiently high. The higher you go in frequency-thickness, the higher `N`  should be. As a rule of thumb: half of the obtained modes will be accurate. The figure below shows the convergence with respect to the Rayleigh-Lamb root of the S1 mode at 5.6 rad/mm in an aluminum plate (solution close to 5 MHz mm). A [Spectral Collocation](https://github.com/dakiefer/GEW_dispersion_script) implementation is shown in comparison. 15 digits accuracy is attained with N = 16 in this case.

![relative error w.r.t. Rayleigh-Lamb root](resourcesAndDeps/img/convergence.png)

For general information on the formulation of the elastic waveguide problem refer to 
> D. A. Kiefer, _Elastodynamic quasi-guided waves for transit-time ultrasonic flow metering_, ser. FAU Forschungen, Reihe B, Medizin, Naturwissenschaft, Technik, vol. 42. Erlangen: FAU University Press, 2022, doi: [10.25593/978-3-96147-550-6](http://doi.org/10.25593/978-3-96147-550-6). [![PDF](resourcesAndDeps/img/icon_file-pdf.svg)](https://dakiefer.net/publication/2022_dissertation_elastodynamic-quasi-guided-waves/2022_dissertation_Elastodynamic%20quasi-guided%20waves.pdf)

For the computation of zero-group-velocity (ZGV) points refer to
> D. A. Kiefer, B. Plestenjak, H. Gravenkamp, and C. Prada, “Computing zero-group-velocity points in anisotropic elastic waveguides: Globally and locally convergent methods,” The Journal of the Acoustical Society of America, vol. 153, no. 2, pp. 1386–1398, Feb. 2023, doi: [10.1121/10.0017252](http://doi.org/10.1121/10.0017252). [![PDF](resourcesAndDeps/img/icon_file-pdf.svg)](https://dakiefer.net/publication/2023_JASA_Computing_ZGV/2023_JASA_Computing_ZGV.pdf).

## Dependencies

GEWtool depends on the functions `collocD` , `lglnodes`, `lgwt` and `barylag` created by Greg von Winckel. They are bundled together with their license files in the `resourcesAndDeps` directory. You may also find them on

> Greg von Winckel, MATLAB Central File Exchange, [https://fr.mathworks.com/matlabcentral/profile/authors/869721](https://fr.mathworks.com/matlabcentral/profile/authors/869721).

The function `computeZGVDirect` depends on the `MultiParEig toolbox` by Bor Plestenjak and Andrej Muhič: 

> Bor Plestenjak (2022). MultiParEig ([https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig](https://www.mathworks.com/matlabcentral/fileexchange/47844-multipareig)), MATLAB Central File Exchange.

## Citing GEWtool

If this code is useful to you, please cite it as (always indicating the DOI):

> D. A. Kiefer (2023). GEWtool. doi: [10.5281/zenodo.10114243](https://doi.org/10.5281/zenodo.10114243) ([https://github.com/dakiefer/GEWtool](https://github.com/dakiefer/GEWtool))

## Contributors

Daniel A. Kiefer, Institut Langevin, ESPCI Paris, Université PSL  
Author and developer

Bor Plestenjak, Faculty of Mathematics and Physics, University of Ljubljana, Slovenia  
Numerical methods to compute ZGV points

Acknowledgments - Many inspiring discussions influenced GEWtool. D. A. Kiefer is thankful to:   
Hauke Gravenkamp, Claire Prada and Michael Ponschab

## Author

2022–2024 – Daniel A. Kiefer, Institut Langevin, ESPCI Paris, Université PSL.

I have several years of experience in waveguide modeling and numerical implementations thereof. In January 2022 I decided to create a new modular and versatile code from scratch. The result is GEWtool. My hope is that it be a valuable research tool and at the same time a helpful educational resource for those interested in numerical methods and elastic waves.

Contact: [daniel.kiefer@espci.fr](mailto:daniel.kiefer@espci.fr) &nbsp; • &nbsp; [dakiefer.net](https://dakiefer.net) &nbsp; • &nbsp; Follow me on [ResearchGate](https://www.researchgate.net/profile/Daniel-Kiefer-5)!

[![Logo Institut Langevin](resourcesAndDeps/img/logo_institut_langevin.svg)](https://www.institut-langevin.espci.fr)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[![Logo ESPCI](resourcesAndDeps/img/logo_espci.svg)](https://www.espci.psl.eu/en/)