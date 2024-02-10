# Changelog 

Documents the changes to GEWtool.

## 1.1.0 (2023-12-23)

- FEATURE Provide targets for the angular frequency or the wavenumbers. You can pass the target in `opts.target` when calling `computeW()` or `computeK()`.
- FEATURE Introduce aliases for common terminology used in cylindrical waveguides: `longitudinal` -> `Lamb(0)`; `torsional` -> `sh(0)`, and `flexural(n)` -> `fullyCoupled(n)`.
- BUGFIX Correct decoupling of longitudinal and torsional waves (this is possible of waves of circumferential order zero).
- BUGFIX Allow to iteratively build matrices of symmetric/anti-symmetric waves. 
- EXAMPLEFILE dispersionSurface_silicon_ZGV.m is more robust when computing for other materials. 

## 1.0 (2023-11-11)

- Initial release.