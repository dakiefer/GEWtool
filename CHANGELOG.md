# Changelog 

Documents the changes to GEWtool.

## 1.2.0 (2024-05-16)

- FEATURE Waves in rods (solid cylinders) can now be computed
- FEATURE Integration on Gauss-Legendre points (used for rods)
- FEATURE Waveguide classes allow to retrieve waves of any polarization (specified as a vector)
- examples: added example for waves in a rod
- examples: updated Lamb waves example to be more minimalistic
- tests: compare waves in a rod to solutions by disperse and also to a paper by Zemanek
- BUGFIX that prevented to specify a scalar number of nodes for multilayered cylinders
- BUGFIX corrected underlying equations in the reference implementations/cylindrical

## 1.1.0 (2023-12-23)

- FEATURE Provide targets for the angular frequency or the wavenumbers. You can pass the target in `opts.target` when calling `computeW()` or `computeK()`.
- FEATURE Introduce aliases for common terminology used in cylindrical waveguides: `longitudinal` -> `Lamb(0)`; `torsional` -> `sh(0)`, and `flexural(n)` -> `fullyCoupled(n)`.
- BUGFIX Correct decoupling of longitudinal and torsional waves (this is possible of waves of circumferential order zero).
- BUGFIX Allow to iteratively build matrices of symmetric/anti-symmetric waves. 
- EXAMPLEFILE dispersionSurface_silicon_ZGV.m is more robust when computing for other materials. 

## 1.0 (2023-11-11)

- Initial release.