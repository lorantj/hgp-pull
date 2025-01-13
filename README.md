# hgp-pull

The Harmonic Guiding Potential pulling tool for non-equilibrium Steered Molecular Dynamics pulls of specified selection(s) along specified path(s) (reaction coordinate(s)) in NAMD [1, 2]. This tool computes on-the-fly in parallel the work(s) along (multiple) independent pulling RC(s) required to determine (multiple) Potential of Mean Force(s) (PMF(s)) along the specified (multiple) independent RC(s).

Requirements: NAMD 2.12 or newer [1, 2] - https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD

# Running a sample of 4 parralel pulls in the GlpF channel

One example of running pulls of water molecules along the GlpF channel is placed in the common directory together with an example of the input configuration file. The output log files will contain additional lines to the typical NAMD log file with the following two lines:
HGP11:
HGP12:
with the specified frequency.

# Analysis - extract works and compute PMF

Firstly, extract the HGP12 line with the work column.
Secondly, convert the timestep to your distance or angle RC based on the constant velocity pull used in the SMD simulation.
Third, calculate the mean and standard deviations for the work in the forward and reverse pulls from a significant number of pulls [3, 4].
Fourth, obtain the PMF along the specified RC(s) according to the FR method [3, 4].

# References:

1. J. C. Phillips, R. Braun, W. Wang, J. Gumbart, E. Tajkhorshid, E. Villa, C. Chipot, R. D. Skeel, L. Kale, and K. Schulten.
Journal of Computational Chemistry 26, 1781 (2005). doi:10.1002/jcc.20289

2. James C. Phillips, David J. Hardy, Julio D. C. Maia, John E. Stone, Joao V. Ribeiro, Rafael C. Bernardi, Ronak Buch, Giacomo Fiorin, Jerome Henin, Wei Jiang, Ryan McGreevy, Marcelo C. R. Melo, Brian K. Radak, Robert D. Skeel, Abhishek Singharoy, Yi Wang, Benoit Roux, Aleksei Aksimentiev, Zaida Luthey-Schulten, Laxmikant V. Kale, Klaus Schulten, Christophe Chipot, and Emad Tajkhorshid.
Scalable molecular dynamics on CPU and GPU architectures with NAMD.
Journal of Chemical Physics, 153, 044130 (2020). doi:10.1063/5.0014475

3. I Kosztin, B Barz, L Janosi.
Calculating potentials of mean force and diffusion coefficients from nonequilibrium processes without Jarzynski’s equality
The Journal of Chemical Physics, 124(6), 064106 (2006). doi:10.1063/1.2166379

4. MW Forney, L Janosi, I Kosztin.
Calculating free-energy profiles in biomolecular systems from fast nonequilibrium processes
Physical Review E—Statistical, Nonlinear, and Soft Matter Physics 78(5), 051913 (2008). doi:10.1103/PhysRevE.78.051913

