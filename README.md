# Physics based wave tank simluations

## Content
This repository contains the files utilized for carrying out Lagrangian n-particle based computations to simulate numerical wave tank experiments in order to investigate rapid wave formation in oceans. The method used is broadly Smoothed Particle Hydrodynamics (SPH) which allows for General Purpose Graphics Processing Unit (GPGPU) computing, allowing for massive parallelization and reducing computational times. It is very well suited for hydrodynamic simulations, as it is meshless and can easily model complex, dynamic and violent flows. It does not require additional effort for tracking free surface particles, which can be a problem for other meshed methods. 
The codes contain the simulation framework for two SPH methods. 
- *Weakly Compressible Smoothed Particle Hydrodynamics (WCSPH)*:
  1. For this, different numerical dissipation schemes and kernel functions were investigated for suitability in carrying out long duration numerical wave tank experiments. Following this, wave tank experiments have been simulated to investigate the effect of modulations and the effect of parameters associated with them on resultant wave amplitdues for a regular wave train moving unidirectionally. 3-D experiments have been also carried out to simulate the effect of frequencies and wave train directions in a crossing sea scenario on the effect of resultant wave amplitdues.
  2. A Peregrine type wave has also been simulated to investigate wave energy localization in our constructed numerical wave tank. 
  
- *Corrected Conservative Smoothed Particle Hydrodynamics (CCSPH)*:
  1. A modification the kernel gradient calculation was implemented to better take into account the effect of the neighboring particles in kernel estimation, specifically near the free surface or domain boundaries, to facilitate better conservation of energy during wave tank simulations. Focusing of double irregular wave groups have been carried out to validate this method against the originally developed WCSPH. Following this, CCSPH has been utilized for simulations of waves using real wave spectra conditions obtained from wave monitoring instruments.

 ##
