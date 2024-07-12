# Changelog 

Documents the changes to GEWtool.

## 1.5 (2024-07-12)

- FEATURE **faster computation of zgv points** with `computeZGVScan()`. It is now based on a Sylvester-Arnoldi method and exploits the structure of the operator determinants of the multi-parameter eigenvalue problem. The k-domain scanning strategy has also been improved.
- FEATURE `GEWdat` functions now support simultaneous postprocessing of **array of waveguide problems** (as returned by `LambSA()`). 
- BUGFIX the excitability of modes is proportional to the particle velocity, not the displacement as previously assumed.
- groupVel(): ignore imaginary part (no meaning). Warn when computing for complex waves.

## 1.4.1 (2024-06-14)

- FEATURE `computeZGV()` and `computeZGVScan()` now allow to pass an array of Waveguide objects. 
- better **backward compatibility** to Matlab releases before R2022a
- BUGFIX in computeZGV() that lead to wrong determination of initial guess. 

## 1.4 (2024-06-14)

- FEATURE `GEWdat` functions now fully support postprocessing in cylindrical coordinates.
- FEATURE `GEWintegrate()`: now allows to specify the layer(s) on which to integrate.
- FEATURE Consistent behavior of `energyVel`, `energyVelVec`, `energyVelAxial`, `energyVelTransverse`, `groupVel`, `groupVelAxial`.
- FEATURE Waveguide class: now remembers circumferential order n of chosen waves in a cylinder (`gew.n`).
- FEATURE `GEWdat`: more efficient power flux and energy computation (exploit equipartition of energy, avoid computing unnecessary components in power flux)
- Cylinder.decouplesLambvsSH() no longer needs the circumferential wavenumber as argument if the waves have already been selected.
- GEWintegrateEachLayer() added: compute integral on each layer separately.
- Plate/displGrad(): compute the displacement gradient in Cartesian coordinates. 
- Cylinder/displGrad(): compute the displacement gradient in cylindrical coordinates. 
- GEWdat/energyVelVec() added: energy velocity vector.
- GEWdat/energyVelMag() added: magnitude of the energy velocity vectors.
- GEWdat/energyVelAxial() added: energy velocity component aligned with the wave vector.
- GEWdat/energyVelTransverse() added: energy velocity component orthogonal to the wave vector.
- GEWdat/energyVel() is now alias to energyVelMag().
- GEWdat/groupVelAxial() added: group velocity component aligned with the wave vector.
- GEWdat/groupVel() now returns the group velocity magnitude when Lamb and sh polarizations decouple. 
- GEWdat/poyntingVecTransverse() added: compute only the ez-component.
- GEWdat/powerFluxMag() added: computes the magnitude of the power flux.
- GEWdat/powerFlux() now returns a vector instead of only the axial component.
- GEWintegrate() now correctly integrates in cylindrical coordinates when an object of class Cylinder is passed.
- GEWdat/strain() and GEWdat/stress() are now based on the Waveguide.displGrad().
- MaterialIsotropic class: now tests range of validity before assigning parameter values.
- eigenVecs() re-assembles the displacement arrays dat.u back into eigenvectors.
- example scripts now respect the default line color order.
- BUGFIX in GEWintegrate() that lead to wrong size of array when computing the power flux vectors.
- BUGFIX in strain() that lead to wrong results for triclinic plates.

## 1.3 (2024-06-07)

- FEATURE faster computations with computeW() and computeK(): new option "eigenvecs" (whether eigenvectors should be computed) and "standardEVP" (whether to convert the generalized eigenvalue problem into a standard eigenvalue problem described by only one matrix).
- FEATURE postprocessing: GEWinterpolate() implements barycentric Lagrange interpolation to obtain the field at any coordinates yi. 
- FEATURE postprocessing: GEWextrapolateCircumference() extrapolates the field in angular direction of a cylindrical waveguide. 
- FEATURE postprocessing: extractModes() allows to reduce the data in the "dat" structure to the desired mode(s).
- FEATURE polarization() allows to choose any arbitrary polarization by selecting the displacement components yourself. 
- FEATURE Material now computes the universal anisotropy index (AU). 
- FEATURE Material uniform Voigt and Reuss homogenization (for anisotropy index).
- FEATURE Material plotSlownessCurve() now supports the standard plot arguments.
- FEATURE Material compute energy velocity of bulk waves. Plot ray curve.
- FEATURE Material rotateEuler() now supports arbitrary rotation sequences.
- FEATURE Material isConservative() test if it is nondissipative.
- FEATURE Material convert quaternions to Euler angles. 
- FEATURE Material rotateVoigtMatrix() can rotate the 6x6 Voigt notated matrix directly. This is faster than rotating the tensor and is useful for fast homogenization.
- Hermitian matrices: now makes sure that the assembled matrices of nondissipative waveguides are all hermitian even in finite precision context. This can make the computations faster, since the LU method can be used instead of QZ method in eig().
- Waveguide now remembers the polarization components of the assembled wave operators.
- Cylinder: more useful error mesages.
- Added example that animates the modal field in a rod.
- Material tests invariance on reflection within 12 digits accuracy.
- Material now warns you when the stiffness is not positive definite or not symmetric.
- BUGFIX Material: transformBasis() now keeps only 12 digits accuracy and avoids breaking the symmetry of the stiffness tensor.
- BUGFIX Waveguide and Plate: correctly keep information on whether the geometry has been symmetrized. Throw an error when it is not possible to symmetrize.

## 1.2 (2024-05-16)

- FEATURE Waves in rods (solid cylinders) can now be computed
- FEATURE Integration on Gauss-Legendre points (used for rods)
- FEATURE Waveguide classes allow to retrieve waves of any polarization (specified as a vector)
- examples: added example for waves in a rod
- examples: updated Lamb waves example to be more minimalistic
- tests: compare waves in a rod to solutions by disperse and also to a paper by Zemanek
- BUGFIX that prevented to specify a scalar number of nodes for multilayered cylinders
- BUGFIX corrected underlying equations in the reference implementations/cylindrical

## 1.1 (2023-12-23)

- FEATURE Provide targets for the angular frequency or the wavenumbers. You can pass the target in `opts.target` when calling `computeW()` or `computeK()`.
- FEATURE Introduce aliases for common terminology used in cylindrical waveguides: `longitudinal` -> `Lamb(0)`; `torsional` -> `sh(0)`, and `flexural(n)` -> `fullyCoupled(n)`.
- BUGFIX Correct decoupling of longitudinal and torsional waves (this is possible of waves of circumferential order zero).
- BUGFIX Allow to iteratively build matrices of symmetric/anti-symmetric waves. 
- EXAMPLEFILE dispersionSurface_silicon_ZGV.m is more robust when computing for other materials. 

## 1.0 (2023-11-11)

- Initial release.