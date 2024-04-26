# Neutron stars and the cosmological constant problem

The NSCC is a code to construct solutions for static, spherically symmetric neutron stars described by a new prescription for equation of states (EOSs), also computed in the code. 
The EOS is modelled by three regions: 
- at low-density, for $\rho<=2\rho_0$, where $\rho_0\simeq2.7 \times 10^{14} \text{g}/\text{cm}^{3}$ is the saturation mass density, we adopt a tabulated EOS.
- at high-density, for $\rho>2\rho_0$ the nuclear matter undergoes a QCD phase transition, and we employ a phenomenological EOS (https://arxiv.org/abs/1804.02783).
- in the inner stellar core, for pressure values larger than $(200, \text{MeV})^4$, we assume that a vacuum energy phase transition occurs, and we use a phenomenological EOS.

The general description for the EOS model is provided in arXiv:23....

In the code we employ units in which G=c=M_sun=1, but rescaling factors are provided.

# Content

Here is a list of the files:
- NSCC.ipynb, the main Julia code to solve the Tolman-Oppenheimer-Volkoff and the tidal deformability equations. EOSs are generated within the code. The code also computes the combined tidal deformabilities.
- eos_generator.ipynb, a Julia code to generate tabulated EOSs which can be used in well-known neutron stars integrator, like the RNS code for fast rotating neutron stars(https://github.com/cgca/rns).
- T4.ipynb, a Julia code that compute the inspiral of neutron stars binaries with post-Newtonian description (employing the Taylor T4 formalism with tidal contributions).
- Rescaledap4.dat and Rescaledsly.dat, tabulated EOS (AP4 and SLy) for the low-density description. They are already rescaled in G=c=M_sun=1 units.
- matrixrhoS.dat and matrixcsS.dat, small datasets of random values for mass density and speed of sound used to construct the phenomenological EOS at high-density.
- matrixrhoL.dat and matrixcsL.dat, large datasets of random values for mass density and speed of sound used to construct the phenomenological EOS at high-density.

**Important notes**: 
- when using the NSCC code there is no need to create the 'modified' EOSs beforehand (e.g. through the eos_generator notebook provided in the repository). The 'modified' EOS are created within the code.
- other tabulated EOS for describing the low-density region can be supplied by the user, but have to carefully be rescaled in the correct units (G=c=M_sun=1).
