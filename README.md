finite differences 3d code
==========================
- triffy.f90: main code and subroutines
- triffy.dec: declarations of variables
- parameters.f: size of grid (do not alter the -2+npm+6 as this is needed for PML boundaries)
- fpar: parameters of simulations that are read on the run
- dx/dt should be proportional to maximum wave velocity in medium to allow for stability of the code. A value dx/dt > 5 vmax is advised, however code may still run without instability in some cases with sqrt(3) > dx/dt 
- dx should be sufficiently small to sample the minimum wavelength correctly, avoiding numerical dispersion.
Usually the surface or Rayleigh wave is the most numerically dispersive, being the slowest and most localised wave most times.  
- medium velocity for P and S wave is sqrt((mu + 2 lambda)/rho))
and sqrt(mu/rho), respectively. mu is shear modulus, lambda is lame modulus, and rho is mass density. 
- inhomogeneous waves (e.g. Rayleigh) can be slower than S wave, depending on the specific ration of mu and lambda.
- additional routines defined in other f90 files are used for PML absorbing boundaries only, and need not be modified.
- free surface is implemented in first j (y) layers
- absorbingg boundaries are implemented on all other sides of cuboid
- an additional [optional] absortpion coefficient is added on edges and first j layers, this decays exponentially with distance from edge, routine "corner_abs".
The aim is to attenuate surface waves that are not sufficiently reduced by PML.
corner_abs is controlled by 'rlen' and 'dissipation' coefficients.
Inadequate value of dissipation can produce float point exceptionsi and consequenct program crash. 
Reduce dissipation in fpar proportionally to dt, or put dissipation to 0 to eliminate the corner absorption, or bypass call to corner_abs.
