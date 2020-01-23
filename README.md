# GammaOUBiGammaOU
Matlab codes relative to the paper

The folder "src" contains the codes with the implementation of the approaches discussed in the paper

- GOUBGOU.m Class that implements the logic of a Gamma-OU and a Bilateral-Gamma-OU and the path-simulation.

- GOUIncrementMethod.m Abstract super class of the following IncrementMethods
- GOUIncrementContTankov.m implements the solution of Cont and Tankov for the increment of a Gamma-OU
- GOUIncrementCufaroSabino.m implements the fast solution proposed in the paper for the increment of a Gamma-OU
- @GOUIncrementCufaroSabinoRejection is the package that implements the acceptance-rejection method proposed in the paper for the increment of a Gamma-OU
- GOUIncrementQDZ.m implements the solution of Qu et Al. for the increment of a Gamma-OU

- BGOUSymmetricIncrementContTankov.m implements the solution of Cont and Tankov for the increment of a Symmetric Bilateral-Gamma-OU
- BGOUSymmetricIncrementCufaroSabino.m implements the fast solution proposed in the paper for the increment of a Symmetric Bilateral-Gamma-OU
- @BGOUSymmetricIncrementCufaroSabinoRejection is the package that implements the acceptance-rejection method proposed in the paper for the increment of a Symmetric Bilateral-Gamma-OU
- BGOUSymmetricIncrementQDZ.m implements the solution of Qu et Al. for the increment of a Symmetric Bilateral-Gamma-OU

- ResolveIncrementMethod.m Resolves which method is to be instantiated

The folder "main" contains the main files to reproduce the tables and figures in the paper

- MainGOU.m reproduces the values in the Tables and in the Figures relative to simulation of a Gamma-OU at one time step only.
- MainGOUTrajectory.m reproduces the values in the Tables relative to simulation of the skeleton of a Gamma-OU over a time grid.
- MainSymmetricGOU.m reproduces the values in the Tables and in the Figures relative to simulation of a Symmetric Bilateral-Gamma-OU at one time step only.
- MainSymmetricGOUTrajectory.m reproduces the values in the Tables relative to simulation of the skeleton of a Symmetric Bilateral-Gamma-OU over a time grid.
