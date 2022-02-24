# Hot Jupiter Evolution code
Numerical code that solves secular equations for the tidal spin-orbit evolution of star+planet system

We consider the creep tide theory to model the tidal interactions. Stellar tides are considered by employing the Constant Rotation Rate approximation of the creep tide theory (see https://arxiv.org/abs/1910.12990 , Section 3). Thus, the stellar rotation value can assume any values (typical initial rotation rate period values are 1.6 days, 2.3 days or 8.0 days, corresponding to fast, moderately fast and slow rotators, respectively). For the planetary tides, we consider the pseudo-synchronous approximate solution presented in https://arxiv.org/abs/1707.09229 , Section 8.3. Thus, the planetary rotation rate is fixed to the pseudo-synchronous solution from the beginning of the simulation. 

The variables which are integrated in the code are:

- Semi-major axis (Stellar and Planetary Tides);
- Eccentricity (Stellar and Planetary Tides);
- Stellar rotation (Stellar Tides and Magnetic Wind Braking of the stellar rotation rate, the latter interaction following https://adsabs.harvard.edu/full/1997A%26A...326.1023B).

A file containing a more detailed description of the code (and some instructions on how to use it) will be added, as well as new scripts to automatically create plots of the results from numerical simulations.
