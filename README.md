# MCRT_Muon_Simulation
Monte Carlo simulations are a powerful computational tool for simulating the transport of particles through matter. They are inherently probabilistic in nature which makes them useful for cascade physics. This method adapts existing code for propagating photons through various stellar objects to model cosmic-ray muons traveling through materials of different Z-number. The code simplifies multiple coulomb scattering into purely elastic collisions using Rutherford scattering. The screened model calculated the standard deviation of the spread of exit angles to be 10.1 milliradians for iron and 20.7 milliradians for lead. Borozdin et al. simulated the spread in iron to be 11 milliradians in iron and 20 milliradians in lead. The results agree with previous simulations using GEANT4 and can be easily modified for different materials and inclusions. A point of closest approach algorithm is used to reconstruct and identify the placement of high Z materials embedded in a volume. However, more advanced reconstructive algorithms should be used for higher accuracy. For the future of imaging nuclear waste, this code provides a simple and accurate tool.

The FORTRAN files are:

	density.f
	gridset.f
        iarray.f
        mcpolar2.f
        ScreenedPot.f
        sourceph.f
        tauint2.f

Include files are:

	grid.txt  
	photon.txt  

Input parameters are in:

	input.params

The file that compiles the code and creates the executable file 'mcgrid' is:

	Makefile

Adapted from Dr. Kenneth Woods
